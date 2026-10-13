// Voir apex/gpu/vk_context.hpp et apex/gpu/vk_internal.hpp.
#include "apex/gpu/vk_context.hpp"

#include <stdexcept>

#ifndef APEX_WITH_VULKAN

namespace apex::gpu {

bool compiled_with_vulkan() { return false; }

GpuStatus status() {
    GpuStatus s;
    s.reason = "backend Vulkan non construit (APEX_WITH_VULKAN absent a la compilation)";
    return s;
}

double saxpy_f32(float, const float*, float*, std::int64_t) {
    throw std::runtime_error(status().reason);
}

double saxpy_f64(double, const double*, double*, std::int64_t) {
    throw std::runtime_error(status().reason);
}

}  // namespace apex::gpu

#else  // APEX_WITH_VULKAN

#include "apex/gpu/vk_internal.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mutex>

namespace apex::gpu {

namespace {

struct Context {
    VkInstance instance = VK_NULL_HANDLE;
    VkPhysicalDevice phys = VK_NULL_HANDLE;
    VkDevice device = VK_NULL_HANDLE;
    VkQueue queue = VK_NULL_HANDLE;
    std::uint32_t queue_family = 0;
    VmaAllocator vma = nullptr;
    VkCommandPool pool = VK_NULL_HANDLE;
    VkQueryPool timestamps = VK_NULL_HANDLE;
    std::uint32_t max_groups_x = 65535;
    GpuStatus st;
    std::mutex mutex;   // une seule soumission a la fois (file unique)
    std::map<std::string, detail::Kernel> kernels;
};

// Volontairement jamais detruit : le module vit autant que le processus, et
// detruire des objets Vulkan au dechargement de la DLL est une source connue
// de plantages a la fermeture.
Context* g_ctx = nullptr;
std::once_flag g_once;

std::string vk_fail(const char* what, VkResult r) {
    return std::string(what) + " a echoue (VkResult " + std::to_string(static_cast<int>(r)) + ")";
}

void check(VkResult r, const char* what) {
    if (r != VK_SUCCESS) throw std::runtime_error(vk_fail(what, r));
}

bool env_flag(const char* name) {
    char* v = nullptr;
    std::size_t len = 0;
    if (_dupenv_s(&v, &len, name) != 0 || v == nullptr) return false;
    const bool on = std::strcmp(v, "1") == 0;
    std::free(v);
    return on;
}

void init_context() {
    auto* c = new Context();
    g_ctx = c;
    VkResult r = volkInitialize();
    if (r != VK_SUCCESS) {
        c->st.reason = "vulkan-1.dll introuvable (aucun pilote Vulkan installe)";
        return;
    }
    VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    app.pApplicationName = "apex_core";
    app.apiVersion = VK_API_VERSION_1_3;
    std::vector<const char*> layers;
    if (env_flag("APEX_VK_VALIDATION")) {
        std::uint32_t count = 0;
        vkEnumerateInstanceLayerProperties(&count, nullptr);
        std::vector<VkLayerProperties> props(count);
        vkEnumerateInstanceLayerProperties(&count, props.data());
        for (const auto& p : props)
            if (std::strcmp(p.layerName, "VK_LAYER_KHRONOS_validation") == 0) layers.push_back("VK_LAYER_KHRONOS_validation");
    }
    VkInstanceCreateInfo ici{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    ici.pApplicationInfo = &app;
    ici.enabledLayerCount = static_cast<std::uint32_t>(layers.size());
    ici.ppEnabledLayerNames = layers.data();
    r = vkCreateInstance(&ici, nullptr, &c->instance);
    if (r != VK_SUCCESS) {
        c->st.reason = vk_fail("vkCreateInstance", r);
        return;
    }
    volkLoadInstanceOnly(c->instance);

    std::uint32_t n = 0;
    vkEnumeratePhysicalDevices(c->instance, &n, nullptr);
    std::vector<VkPhysicalDevice> devices(n);
    vkEnumeratePhysicalDevices(c->instance, &n, devices.data());
    // GPU discret, NVIDIA d'abord (cible actuelle : RTX 3090).
    int best_rank = -1;
    for (VkPhysicalDevice d : devices) {
        VkPhysicalDeviceProperties p;
        vkGetPhysicalDeviceProperties(d, &p);
        if (p.deviceType != VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU || p.apiVersion < VK_API_VERSION_1_3) continue;
        const int rank = p.vendorID == 0x10DE ? 2 : 1;
        if (rank > best_rank) { best_rank = rank; c->phys = d; }
    }
    if (c->phys == VK_NULL_HANDLE) {
        c->st.reason = "aucun GPU discret compatible Vulkan 1.3";
        return;
    }
    VkPhysicalDeviceProperties props;
    vkGetPhysicalDeviceProperties(c->phys, &props);
    c->st.device = props.deviceName;
    c->max_groups_x = props.limits.maxComputeWorkGroupCount[0];

    std::uint32_t qn = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(c->phys, &qn, nullptr);
    std::vector<VkQueueFamilyProperties> qprops(qn);
    vkGetPhysicalDeviceQueueFamilyProperties(c->phys, &qn, qprops.data());
    bool found = false;
    for (std::uint32_t i = 0; i < qn && !found; ++i)
        if (qprops[i].queueFlags & VK_QUEUE_COMPUTE_BIT) { c->queue_family = i; found = true; }
    if (!found) {
        c->st.reason = "aucune file compute";
        return;
    }
    const bool has_timestamps = qprops[c->queue_family].timestampValidBits > 0 && props.limits.timestampPeriod > 0;

    VkPhysicalDeviceFeatures supported;
    vkGetPhysicalDeviceFeatures(c->phys, &supported);
    VkPhysicalDeviceFeatures enabled{};
    enabled.shaderFloat64 = supported.shaderFloat64;
    enabled.shaderInt64 = supported.shaderInt64;
    const float prio = 1.0f;
    VkDeviceQueueCreateInfo qci{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    qci.queueFamilyIndex = c->queue_family;
    qci.queueCount = 1;
    qci.pQueuePriorities = &prio;
    VkDeviceCreateInfo dci{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    dci.queueCreateInfoCount = 1;
    dci.pQueueCreateInfos = &qci;
    dci.pEnabledFeatures = &enabled;
    r = vkCreateDevice(c->phys, &dci, nullptr, &c->device);
    if (r != VK_SUCCESS) {
        c->st.reason = vk_fail("vkCreateDevice", r);
        return;
    }
    volkLoadDevice(c->device);
    vkGetDeviceQueue(c->device, c->queue_family, 0, &c->queue);

    VmaVulkanFunctions fn{};
    fn.vkGetInstanceProcAddr = vkGetInstanceProcAddr;
    fn.vkGetDeviceProcAddr = vkGetDeviceProcAddr;
    VmaAllocatorCreateInfo aci{};
    aci.vulkanApiVersion = VK_API_VERSION_1_3;
    aci.physicalDevice = c->phys;
    aci.device = c->device;
    aci.instance = c->instance;
    aci.pVulkanFunctions = &fn;
    r = vmaCreateAllocator(&aci, &c->vma);
    if (r != VK_SUCCESS) {
        c->st.reason = vk_fail("vmaCreateAllocator", r);
        return;
    }
    VkCommandPoolCreateInfo pci{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    pci.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    pci.queueFamilyIndex = c->queue_family;
    r = vkCreateCommandPool(c->device, &pci, nullptr, &c->pool);
    if (r != VK_SUCCESS) {
        c->st.reason = vk_fail("vkCreateCommandPool", r);
        return;
    }
    if (has_timestamps) {
        VkQueryPoolCreateInfo qpi{VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
        qpi.queryType = VK_QUERY_TYPE_TIMESTAMP;
        qpi.queryCount = 2;
        if (vkCreateQueryPool(c->device, &qpi, nullptr, &c->timestamps) == VK_SUCCESS)
            c->st.timestamp_period_ns = props.limits.timestampPeriod;
    }
    c->st.fp64 = enabled.shaderFloat64 == VK_TRUE;
    c->st.ok = true;
}

Context& ctx() {
    std::call_once(g_once, init_context);
    if (!g_ctx->st.ok) throw std::runtime_error(g_ctx->st.reason);
    return *g_ctx;
}

// Commande a usage unique : enregistrer, soumettre, attendre (barriere, 60 s au plus).
template <class Record>
void submit_and_wait(Context& c, Record&& record) {
    VkCommandBufferAllocateInfo ai{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    ai.commandPool = c.pool;
    ai.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    ai.commandBufferCount = 1;
    VkCommandBuffer cb = VK_NULL_HANDLE;
    check(vkAllocateCommandBuffers(c.device, &ai, &cb), "vkAllocateCommandBuffers");
    VkCommandBufferBeginInfo bi{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    bi.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    vkBeginCommandBuffer(cb, &bi);
    record(cb);
    vkEndCommandBuffer(cb);
    VkFenceCreateInfo fci{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    VkFence fence = VK_NULL_HANDLE;
    check(vkCreateFence(c.device, &fci, nullptr, &fence), "vkCreateFence");
    VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    si.commandBufferCount = 1;
    si.pCommandBuffers = &cb;
    VkResult r = vkQueueSubmit(c.queue, 1, &si, fence);
    if (r == VK_SUCCESS) r = vkWaitForFences(c.device, 1, &fence, VK_TRUE, 60ull * 1000 * 1000 * 1000);
    vkDestroyFence(c.device, fence, nullptr);
    vkFreeCommandBuffers(c.device, c.pool, 1, &cb);
    check(r, "soumission GPU");
}

void barrier(VkCommandBuffer cb, VkPipelineStageFlags src, VkAccessFlags src_access, VkPipelineStageFlags dst,
             VkAccessFlags dst_access) {
    VkMemoryBarrier mb{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
    mb.srcAccessMask = src_access;
    mb.dstAccessMask = dst_access;
    vkCmdPipelineBarrier(cb, src, dst, 0, 1, &mb, 0, nullptr, 0, nullptr);
}

// Tampon intermediaire visible de l'hote (transferts).
struct Staging {
    VkBuffer buffer = VK_NULL_HANDLE;
    VmaAllocation alloc = nullptr;
    void* mapped = nullptr;
    Staging(Context& c, std::size_t bytes, bool readback) {
        VkBufferCreateInfo bi{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
        bi.size = bytes;
        bi.usage = readback ? VK_BUFFER_USAGE_TRANSFER_DST_BIT : VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
        VmaAllocationCreateInfo ai{};
        ai.usage = VMA_MEMORY_USAGE_AUTO;
        ai.flags = VMA_ALLOCATION_CREATE_MAPPED_BIT |
                   (readback ? VMA_ALLOCATION_CREATE_HOST_ACCESS_RANDOM_BIT
                             : VMA_ALLOCATION_CREATE_HOST_ACCESS_SEQUENTIAL_WRITE_BIT);
        VmaAllocationInfo info{};
        check(vmaCreateBuffer(c.vma, &bi, &ai, &buffer, &alloc, &info), "vmaCreateBuffer (intermediaire)");
        mapped = info.pMappedData;
    }
    ~Staging() {
        if (buffer != VK_NULL_HANDLE) vmaDestroyBuffer(g_ctx->vma, buffer, alloc);
    }
};

}  // namespace

bool compiled_with_vulkan() { return true; }

GpuStatus status() {
    std::call_once(g_once, init_context);
    return g_ctx->st;
}

namespace detail {

Buffer::Buffer(std::size_t bytes) : bytes_(bytes) {
    Context& c = ctx();
    VkBufferCreateInfo bi{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    bi.size = bytes > 0 ? bytes : 4;
    bi.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    VmaAllocationCreateInfo ai{};
    ai.usage = VMA_MEMORY_USAGE_AUTO_PREFER_DEVICE;
    check(vmaCreateBuffer(c.vma, &bi, &ai, &buffer_, &alloc_, nullptr), "vmaCreateBuffer");
}

Buffer::Buffer(Buffer&& o) noexcept : buffer_(o.buffer_), alloc_(o.alloc_), bytes_(o.bytes_) {
    o.buffer_ = VK_NULL_HANDLE;
    o.alloc_ = nullptr;
    o.bytes_ = 0;
}

Buffer& Buffer::operator=(Buffer&& o) noexcept {
    if (this != &o) {
        release();
        buffer_ = o.buffer_; alloc_ = o.alloc_; bytes_ = o.bytes_;
        o.buffer_ = VK_NULL_HANDLE; o.alloc_ = nullptr; o.bytes_ = 0;
    }
    return *this;
}

Buffer::~Buffer() { release(); }

void Buffer::release() {
    if (buffer_ != VK_NULL_HANDLE) vmaDestroyBuffer(g_ctx->vma, buffer_, alloc_);
    buffer_ = VK_NULL_HANDLE;
    alloc_ = nullptr;
}

void Buffer::upload(const void* data, std::size_t bytes) {
    if (bytes == 0) return;
    Context& c = ctx();
    std::lock_guard<std::mutex> lock(c.mutex);
    Staging s(c, bytes, false);
    std::memcpy(s.mapped, data, bytes);
    vmaFlushAllocation(c.vma, s.alloc, 0, VK_WHOLE_SIZE);
    submit_and_wait(c, [&](VkCommandBuffer cb) {
        VkBufferCopy region{0, 0, bytes};
        vkCmdCopyBuffer(cb, s.buffer, buffer_, 1, &region);
    });
}

void Buffer::download(void* data, std::size_t bytes) const {
    if (bytes == 0) return;
    Context& c = ctx();
    std::lock_guard<std::mutex> lock(c.mutex);
    Staging s(c, bytes, true);
    submit_and_wait(c, [&](VkCommandBuffer cb) {
        VkBufferCopy region{0, 0, bytes};
        vkCmdCopyBuffer(cb, buffer_, s.buffer, 1, &region);
    });
    vmaInvalidateAllocation(c.vma, s.alloc, 0, VK_WHOLE_SIZE);
    std::memcpy(data, s.mapped, bytes);
}

const Kernel& kernel(const std::string& name, const std::uint32_t* spirv, std::size_t words, std::uint32_t buffers,
                     std::uint32_t push_bytes) {
    Context& c = ctx();
    std::lock_guard<std::mutex> lock(c.mutex);
    auto it = c.kernels.find(name);
    if (it != c.kernels.end()) return it->second;
    Kernel k;
    k.buffers = buffers;
    k.push_bytes = push_bytes;
    VkShaderModuleCreateInfo smi{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    smi.codeSize = words * 4;
    smi.pCode = spirv;
    check(vkCreateShaderModule(c.device, &smi, nullptr, &k.module), "vkCreateShaderModule");
    std::vector<VkDescriptorSetLayoutBinding> bindings(buffers);
    for (std::uint32_t i = 0; i < buffers; ++i) {
        bindings[i] = {};
        bindings[i].binding = i;
        bindings[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        bindings[i].descriptorCount = 1;
        bindings[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    }
    VkDescriptorSetLayoutCreateInfo dli{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    dli.bindingCount = buffers;
    dli.pBindings = bindings.data();
    check(vkCreateDescriptorSetLayout(c.device, &dli, nullptr, &k.set_layout), "vkCreateDescriptorSetLayout");
    VkPushConstantRange range{VK_SHADER_STAGE_COMPUTE_BIT, 0, push_bytes};
    VkPipelineLayoutCreateInfo pli{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    pli.setLayoutCount = 1;
    pli.pSetLayouts = &k.set_layout;
    pli.pushConstantRangeCount = push_bytes > 0 ? 1 : 0;
    pli.pPushConstantRanges = &range;
    check(vkCreatePipelineLayout(c.device, &pli, nullptr, &k.layout), "vkCreatePipelineLayout");
    VkComputePipelineCreateInfo cpi{VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
    cpi.stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    cpi.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
    cpi.stage.module = k.module;
    cpi.stage.pName = "main";
    cpi.layout = k.layout;
    check(vkCreateComputePipelines(c.device, VK_NULL_HANDLE, 1, &cpi, nullptr, &k.pipeline),
          "vkCreateComputePipelines");
    return c.kernels.emplace(name, k).first->second;
}

std::uint32_t max_groups() { return ctx().max_groups_x; }

double dispatch(const Kernel& k, const std::vector<const Buffer*>& buffers, const void* push, std::uint32_t groups) {
    Context& c = ctx();
    if (buffers.size() != k.buffers) throw std::invalid_argument("dispatch : nombre de buffers incoherent");
    std::lock_guard<std::mutex> lock(c.mutex);
    VkDescriptorPoolSize ps{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, k.buffers > 0 ? k.buffers : 1};
    VkDescriptorPoolCreateInfo dpi{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    dpi.maxSets = 1;
    dpi.poolSizeCount = 1;
    dpi.pPoolSizes = &ps;
    VkDescriptorPool dpool = VK_NULL_HANDLE;
    check(vkCreateDescriptorPool(c.device, &dpi, nullptr, &dpool), "vkCreateDescriptorPool");
    VkDescriptorSetAllocateInfo dai{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    dai.descriptorPool = dpool;
    dai.descriptorSetCount = 1;
    dai.pSetLayouts = &k.set_layout;
    VkDescriptorSet set = VK_NULL_HANDLE;
    VkResult r = vkAllocateDescriptorSets(c.device, &dai, &set);
    if (r != VK_SUCCESS) {
        vkDestroyDescriptorPool(c.device, dpool, nullptr);
        throw std::runtime_error(vk_fail("vkAllocateDescriptorSets", r));
    }
    std::vector<VkDescriptorBufferInfo> infos(k.buffers);
    std::vector<VkWriteDescriptorSet> writes(k.buffers);
    for (std::uint32_t i = 0; i < k.buffers; ++i) {
        infos[i] = {buffers[i]->handle(), 0, VK_WHOLE_SIZE};
        writes[i] = {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        writes[i].dstSet = set;
        writes[i].dstBinding = i;
        writes[i].descriptorCount = 1;
        writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[i].pBufferInfo = &infos[i];
    }
    vkUpdateDescriptorSets(c.device, k.buffers, writes.data(), 0, nullptr);
    const bool timed = c.timestamps != VK_NULL_HANDLE;
    try {
        submit_and_wait(c, [&](VkCommandBuffer cb) {
            barrier(cb, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_ACCESS_TRANSFER_WRITE_BIT,
                    VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT);
            if (timed) {
                vkCmdResetQueryPool(cb, c.timestamps, 0, 2);
                vkCmdWriteTimestamp(cb, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, c.timestamps, 0);
            }
            vkCmdBindPipeline(cb, VK_PIPELINE_BIND_POINT_COMPUTE, k.pipeline);
            vkCmdBindDescriptorSets(cb, VK_PIPELINE_BIND_POINT_COMPUTE, k.layout, 0, 1, &set, 0, nullptr);
            if (k.push_bytes > 0)
                vkCmdPushConstants(cb, k.layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, k.push_bytes, push);
            vkCmdDispatch(cb, groups, 1, 1);
            if (timed) vkCmdWriteTimestamp(cb, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, c.timestamps, 1);
            barrier(cb, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_WRITE_BIT,
                    VK_PIPELINE_STAGE_TRANSFER_BIT, VK_ACCESS_TRANSFER_READ_BIT);
        });
    } catch (...) {
        vkDestroyDescriptorPool(c.device, dpool, nullptr);
        throw;
    }
    vkDestroyDescriptorPool(c.device, dpool, nullptr);
    if (!timed) return -1.0;
    std::uint64_t ts[2] = {0, 0};
    if (vkGetQueryPoolResults(c.device, c.timestamps, 0, 2, sizeof(ts), ts, sizeof(std::uint64_t),
                              VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT) != VK_SUCCESS)
        return -1.0;
    return static_cast<double>(ts[1] - ts[0]) * c.st.timestamp_period_ns * 1e-9;
}

}  // namespace detail

namespace {

const std::uint32_t kSaxpyF32[] = {
#include "saxpy_f32.spv.inc"
};
const std::uint32_t kSaxpyF64[] = {
#include "saxpy_f64.spv.inc"
};

template <class T, class Push>
double saxpy(const char* name, const std::uint32_t* spirv, std::size_t words, T a, const T* x, T* y,
             std::int64_t n) {
    if (n <= 0) return 0.0;
    const std::size_t bytes = static_cast<std::size_t>(n) * sizeof(T);
    detail::Buffer bx(bytes), by(bytes);
    bx.upload(x, bytes);
    by.upload(y, bytes);
    const detail::Kernel& k = detail::kernel(name, spirv, words, 2, sizeof(Push));
    // Decoupage en dispatchs (limite materielle de groupes et seuil de coupure Windows).
    const std::uint32_t per = std::min<std::uint32_t>(detail::max_groups(), 4096u);
    double gpu_s = 0.0;
    for (std::int64_t start = 0; start < n; start += static_cast<std::int64_t>(per) * 256) {
        const std::int64_t count = std::min<std::int64_t>(n - start, static_cast<std::int64_t>(per) * 256);
        Push p{a, static_cast<std::uint32_t>(n), static_cast<std::uint32_t>(start)};
        const double t = detail::dispatch(k, {&bx, &by}, &p, static_cast<std::uint32_t>((count + 255) / 256));
        gpu_s = (t < 0 || gpu_s < 0) ? -1.0 : gpu_s + t;
    }
    by.download(y, bytes);
    return gpu_s;
}

struct PushF32 { float a; std::uint32_t n, offset; };
struct PushF64 { double a; std::uint32_t n, offset; };

}  // namespace

double saxpy_f32(float a, const float* x, float* y, std::int64_t n) {
    return saxpy<float, PushF32>("saxpy_f32", kSaxpyF32, sizeof(kSaxpyF32) / 4, a, x, y, n);
}

double saxpy_f64(double a, const double* x, double* y, std::int64_t n) {
    if (!status().fp64) throw std::runtime_error("shaderFloat64 indisponible sur ce GPU");
    return saxpy<double, PushF64>("saxpy_f64", kSaxpyF64, sizeof(kSaxpyF64) / 4, a, x, y, n);
}

}  // namespace apex::gpu

#endif  // APEX_WITH_VULKAN
