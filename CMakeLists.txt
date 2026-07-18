cmake_minimum_required(VERSION 3.24)

project(
  MitoArchitect
  VERSION 0.4.1
  DESCRIPTION "Long-read mitochondrial DNA analysis toolkit"
  LANGUAGES C CXX
)

option(MITO_BUILD_TESTS "Build C++ tests" ON)
option(MITO_BUILD_BENCHMARKS "Build C++ benchmarks" OFF)
option(MITO_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

if(MITO_ENABLE_SANITIZERS)
  if(MSVC)
    message(FATAL_ERROR "MITO_ENABLE_SANITIZERS currently requires GCC or Clang")
  endif()
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
endif()

include(CTest)
include(cmake/CompilerWarnings.cmake)

add_subdirectory(core)
add_subdirectory(ffi)
