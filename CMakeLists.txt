cmake_minimum_required(VERSION 3.24)

project(
  MitoArchitect
  VERSION 0.2.0
  DESCRIPTION "Long-read mitochondrial DNA analysis toolkit"
  LANGUAGES C CXX
)

option(MITO_BUILD_TESTS "Build C++ tests" ON)
option(MITO_BUILD_BENCHMARKS "Build C++ benchmarks" OFF)

include(CTest)
include(cmake/CompilerWarnings.cmake)

add_subdirectory(core)
add_subdirectory(ffi)
