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

find_package(Threads REQUIRED)
target_link_libraries(tetranet PUBLIC Threads::Threads)

add_executable(tetranet_cli tools/tetranet_cli.cpp)
target_link_libraries(tetranet_cli PRIVATE tetranet)

enable_testing()
add_executable(tetranet_smoke tests/smoke.cpp)
target_link_libraries(tetranet_smoke PRIVATE tetranet)
add_test(NAME smoke COMMAND tetranet_smoke)
