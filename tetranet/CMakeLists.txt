cmake_minimum_required(VERSION 3.20)
project(tetranet VERSION 0.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(MSVC)
  add_compile_options(/O2 /fp:precise)
else()
  add_compile_options(-O2 -ffp-contract=off)
endif()

add_library(tetranet STATIC src/tetranet.cpp)
target_include_directories(tetranet PUBLIC include)

enable_testing()
add_executable(tetranet_smoke tests/smoke.cpp)
target_link_libraries(tetranet_smoke PRIVATE tetranet)
add_test(NAME smoke COMMAND tetranet_smoke)
