add_library(mito_ffi STATIC
  src/mito_c_api.cpp
)

target_include_directories(mito_ffi
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(mito_ffi PUBLIC mito_core)
target_compile_features(mito_ffi PUBLIC cxx_std_20)
mito_enable_warnings(mito_ffi)

install(TARGETS mito_ffi ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

