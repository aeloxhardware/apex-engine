// Implementation de VulkanMemoryAllocator (en-tete seul, MIT), compilee une fois.
// Fonctions Vulkan fournies a l'execution par volk (VMA_DYNAMIC_VULKAN_FUNCTIONS) :
// aucun lien statique avec vulkan-1.lib.
#include <volk.h>

#define VMA_IMPLEMENTATION
#define VMA_STATIC_VULKAN_FUNCTIONS 0
#define VMA_DYNAMIC_VULKAN_FUNCTIONS 1
#include <vk_mem_alloc.h>
