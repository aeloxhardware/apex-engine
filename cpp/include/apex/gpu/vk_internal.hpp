// Primitives internes du backend Vulkan (compilees seulement avec APEX_WITH_VULKAN).
// Reutilisees par chaque noyau : un noyau = un shader .comp (cpp/shaders/), des
// buffers de stockage (binding 0..n-1, dans l'ordre), des constantes poussees.
#pragma once

#ifdef APEX_WITH_VULKAN

#include <volk.h>
#include <vk_mem_alloc.h>

#include <cstdint>
#include <string>
#include <vector>

namespace apex::gpu::detail {

// Buffer de stockage en memoire du GPU ; transferts par tampon intermediaire.
class Buffer {
public:
    Buffer() = default;
    explicit Buffer(std::size_t bytes);
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
    Buffer(Buffer&& o) noexcept;
    Buffer& operator=(Buffer&& o) noexcept;
    ~Buffer();
    void upload(const void* data, std::size_t bytes);
    void download(void* data, std::size_t bytes) const;
    VkBuffer handle() const { return buffer_; }
    std::size_t size() const { return bytes_; }

private:
    void release();
    VkBuffer buffer_ = VK_NULL_HANDLE;
    VmaAllocation alloc_ = nullptr;
    std::size_t bytes_ = 0;
};

// Pipeline de calcul construit une fois par shader (cache par nom).
struct Kernel {
    VkShaderModule module = VK_NULL_HANDLE;
    VkDescriptorSetLayout set_layout = VK_NULL_HANDLE;
    VkPipelineLayout layout = VK_NULL_HANDLE;
    VkPipeline pipeline = VK_NULL_HANDLE;
    std::uint32_t buffers = 0;
    std::uint32_t push_bytes = 0;
};

const Kernel& kernel(const std::string& name, const std::uint32_t* spirv, std::size_t words, std::uint32_t buffers,
                     std::uint32_t push_bytes);

// Un dispatch de `groups` groupes, les buffers lies dans l'ordre. `push` : les
// constantes (taille push_bytes du noyau). Rend le temps GPU (timestamps, s) ou -1.
// Bloquant (attente de la barriere). Leve std::runtime_error en cas d'echec Vulkan.
double dispatch(const Kernel& k, const std::vector<const Buffer*>& buffers, const void* push,
                std::uint32_t groups);

// Nombre maximal de groupes par dispatch (limite materielle sur x).
std::uint32_t max_groups();

}  // namespace apex::gpu::detail

#endif  // APEX_WITH_VULKAN
