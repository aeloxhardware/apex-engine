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
  src/core/mesh_repair.cpp
  src/core/mesh_intersect.cpp
  src/core/interpolate.cpp
  src/core/mesh_distance.cpp
  bindings/module.cpp
)

target_include_directories(apex_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Eigen (MPL-2.0), header-only, vendore dans third_party/ - voir CLAUDE.md
# limite #56 pour le choix de licence. Introduit au pas 7 (step4) seulement :
# aucun pas anterieur n'en avait besoin, et une dependance de build ne se fige
# pas avant que le code qui l'utilise existe.
target_include_directories(apex_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen)

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
set_target_properties(apex_core PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY         ${CMAKE_CURRENT_SOURCE_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_CURRENT_SOURCE_DIR}/lib
  RUNTIME_OUTPUT_DIRECTORY         ${CMAKE_CURRENT_SOURCE_DIR}/lib
  RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/lib
  RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_CURRENT_SOURCE_DIR}/lib
)
