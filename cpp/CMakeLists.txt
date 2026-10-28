cmake_minimum_required(VERSION 3.20)
project(apex_core LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(apex_core
  src/core/printed_index.cpp
  src/core/kinematics.cpp
  src/core/support.cpp
  src/core/tet_topology.cpp
  src/core/kdtree.cpp
  src/core/bvh.cpp
  src/core/field_solver.cpp
  src/core/linear_solver.cpp
  src/core/heat_method.cpp
  src/core/reorientation.cpp
  src/core/layer_order.cpp
  src/core/step8_filters.cpp
  src/core/step8_assembly.cpp
  src/core/mesh_repair.cpp
  src/core/mesh_intersect.cpp
  src/core/interpolate.cpp
  src/core/mesh_distance.cpp
  src/core/isosurface.cpp
  src/core/geometry2d.cpp
  src/core/past_material.cpp
  src/core/gcode_export.cpp
  src/core/toolpath2d.cpp
  src/core/step6_layer.cpp
  src/core/overhang.cpp
  src/core/segmentation.cpp
  src/core/scheduling.cpp
  src/core/skin_check.cpp
  src/gpu/vk_context.cpp
  third_party/clipper2/src/clipper.engine.cpp
  third_party/clipper2/src/clipper.offset.cpp
  third_party/clipper2/src/clipper.rectclip.cpp
  bindings/module.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../tetranet/src/tetranet.cpp
)

target_include_directories(apex_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../tetranet/include)
find_package(Threads REQUIRED)
target_link_libraries(apex_core PRIVATE Threads::Threads)

# Eigen (MPL-2.0), header-only, vendore dans third_party/ - voir CLAUDE.md
# limite #56 pour le choix de licence. Introduit au pas 7 (step4) seulement :
# aucun pas anterieur n'en avait besoin, et une dependance de build ne se fige
# pas avant que le code qui l'utilise existe.
target_include_directories(apex_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen)

# Clipper2 2.0.1 (Boost Software License 1.0), vendore : offsets polygonaux de
# step6 (pas 15), a la place de shapely/GEOS.
target_include_directories(apex_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/clipper2/include)

# Backend Vulkan compute (phase 1 du plan « portage C++, puis Vulkan »,
# 2026-10-13). OPTIONNEL : sans lui, src/gpu/vk_context.cpp rend un statut
# « non construit » et tout reste sur CPU. Aucun lien avec vulkan-1.lib : volk
# (MIT) charge vulkan-1.dll a l'execution, le module s'importe donc sur une
# machine sans pilote Vulkan. VMA (MIT) pour la memoire. Shaders compiles en
# SPIR-V par glslc (SDK) au build, inclus comme listes de nombres (-mfmt=num).
option(APEX_WITH_VULKAN "Backend Vulkan compute" OFF)
if(APEX_WITH_VULKAN)
  enable_language(C)
  find_package(Vulkan REQUIRED COMPONENTS glslc)
  set(APEX_SHADER_OUT ${CMAKE_CURRENT_BINARY_DIR}/shaders)
  file(MAKE_DIRECTORY ${APEX_SHADER_OUT})
  file(GLOB APEX_SHADER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.comp)
  set(APEX_SPIRV_INCS)
  foreach(src ${APEX_SHADER_SOURCES})
    get_filename_component(name ${src} NAME_WE)
    set(out ${APEX_SHADER_OUT}/${name}.spv.inc)
    add_custom_command(OUTPUT ${out}
      COMMAND Vulkan::glslc --target-env=vulkan1.3 -O -mfmt=num -o ${out} ${src}
      DEPENDS ${src} VERBATIM)
    list(APPEND APEX_SPIRV_INCS ${out})
  endforeach()
  add_custom_target(apex_shaders DEPENDS ${APEX_SPIRV_INCS})
  add_dependencies(apex_core apex_shaders)
  target_sources(apex_core PRIVATE third_party/volk/volk.c src/gpu/vma_impl.cpp)
  target_include_directories(apex_core PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/volk
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/vma/include
    ${APEX_SHADER_OUT})
  target_link_libraries(apex_core PRIVATE Vulkan::Headers)
  target_compile_definitions(apex_core PRIVATE APEX_WITH_VULKAN=1 VK_NO_PROTOTYPES)
endif()

if(MSVC)
  # /fp:precise est OBLIGATOIRE et non negociable : tout le portage repose sur
  # une comparaison d'equivalence avec la version Python. /fp:fast reordonne
  # les operations flottantes et rendrait cet oracle inutilisable.
  target_compile_options(apex_core PRIVATE /O2 /fp:precise /Zc:__cplusplus /permissive-)
else()
  target_compile_options(apex_core PRIVATE -O2 -ffp-contract=off)
endif()

# Le module est depose dans cpp/lib/, ajoute au sys.path par core/_cpp.py.
# Volontairement hors de la racine du depot : pas d'artefact binaire melange
# aux sources Python.
# APEX_LIB_DIR permet de construire ailleurs (ex. cpp/build_dev/lib) pendant
# qu'un calcul long tient le .pyd de cpp/lib verrouille (Windows) ; le
# chargeur le suit via la variable d'environnement APEX_CPP_LIB_DIR.
set(APEX_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib" CACHE PATH "Dossier de sortie du module apex_core")
set_target_properties(apex_core PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY         ${APEX_LIB_DIR}
  LIBRARY_OUTPUT_DIRECTORY_RELEASE ${APEX_LIB_DIR}
  LIBRARY_OUTPUT_DIRECTORY_DEBUG   ${APEX_LIB_DIR}
  RUNTIME_OUTPUT_DIRECTORY         ${APEX_LIB_DIR}
  RUNTIME_OUTPUT_DIRECTORY_RELEASE ${APEX_LIB_DIR}
  RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${APEX_LIB_DIR}
)

# Worker TetGen : executable SEPARE (AGPL-3.0), jamais lie dans apex_core.
# Les crashs natifs de TetGen restent isoles du processus appelant (limites #35
# et #56). Compile avec TETLIBRARY comme le binding Python.
add_executable(apex_tetgen_worker
  tools/tetgen_worker.cpp
  third_party/tetgen/tetgen.cxx
  third_party/tetgen/predicates.cxx
)
target_compile_definitions(apex_tetgen_worker PRIVATE TETLIBRARY)
target_include_directories(apex_tetgen_worker PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/tetgen)
set_target_properties(apex_tetgen_worker PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${APEX_LIB_DIR}
  RUNTIME_OUTPUT_DIRECTORY_RELEASE ${APEX_LIB_DIR})
