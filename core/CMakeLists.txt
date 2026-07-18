add_library(mito_core STATIC
  src/analysis_engine.cpp
)

set(MITO_HDBSCAN_CPP_ROOT "" CACHE PATH "Optional path to rohanmohapatra/hdbscan-cpp HDBSCAN-CPP directory")

find_package(Threads REQUIRED)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(HTSLIB QUIET htslib)
endif()

target_include_directories(mito_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(mito_core PUBLIC cxx_std_20)
target_compile_definitions(mito_core
  PRIVATE
    MITO_DEFAULT_RCRS_FASTA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/rcrs.fasta"
    MITO_CLINICAL_ANNOTATIONS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/clinical_annotations.tsv"
    MITO_PHYLOTREE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/phylotree-rcrs-17.3.xml"
    MITO_PHYLOTREE_WEIGHTS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/phylotree-rcrs-17.3-weights.txt"
    MITO_PHYLOTREE_ALIGNMENT_RULES_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/phylotree-rcrs-17.3-rules.csv"
    MITO_RESOURCE_MANIFEST_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/resource_manifest.tsv"
)
target_link_libraries(mito_core PUBLIC Threads::Threads)
if(HTSLIB_FOUND)
  target_compile_definitions(mito_core PUBLIC MITO_HAS_HTSLIB=1)
  target_include_directories(mito_core SYSTEM PRIVATE ${HTSLIB_INCLUDE_DIRS})
  target_link_libraries(mito_core PUBLIC ${HTSLIB_LIBRARIES})
  target_link_directories(mito_core PUBLIC ${HTSLIB_LIBRARY_DIRS})
endif()

if(MITO_HDBSCAN_CPP_ROOT AND EXISTS "${MITO_HDBSCAN_CPP_ROOT}/Hdbscan/hdbscan.hpp")
  add_library(mito_hdbscan_cpp STATIC
    "${MITO_HDBSCAN_CPP_ROOT}/Distance/EuclideanDistance.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/Distance/IDistanceCalculator.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/Distance/ManhattanDistance.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/Hdbscan/hdbscan.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/HdbscanStar/cluster.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/HdbscanStar/hdbscanAlgorithm.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/HdbscanStar/hdbscanConstraint.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/HdbscanStar/outlierScore.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/HdbscanStar/undirectedGraph.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/Runner/hdbscanParameters.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/Runner/hdbscanResult.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/Runner/hdbscanRunner.cpp"
    "${MITO_HDBSCAN_CPP_ROOT}/Utils/bitSet.cpp"
  )
  target_include_directories(mito_hdbscan_cpp PUBLIC "${MITO_HDBSCAN_CPP_ROOT}")
  target_compile_features(mito_hdbscan_cpp PUBLIC cxx_std_14)
  target_compile_options(mito_hdbscan_cpp PRIVATE -w)
  target_compile_definitions(mito_core PRIVATE MITO_HAS_HDBSCAN_CPP=1)
  target_include_directories(mito_core SYSTEM PRIVATE "${MITO_HDBSCAN_CPP_ROOT}")
  target_link_libraries(mito_core PRIVATE mito_hdbscan_cpp)
  install(TARGETS mito_hdbscan_cpp ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
endif()
mito_enable_warnings(mito_core)

install(TARGETS mito_core ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
install(FILES
  data/rcrs.fasta
  data/clinical_annotations.tsv
  data/phylotree-rcrs-17.3.xml
  data/phylotree-rcrs-17.3-weights.txt
  data/phylotree-rcrs-17.3-rules.csv
  data/phylotree-rcrs-17.3.LICENSE
  data/resource_manifest.tsv
  DESTINATION share/mito-architect/reference
)

if(MITO_BUILD_TESTS AND BUILD_TESTING)
  add_executable(mito_core_smoke_test tests/core_smoke_test.cpp)
  target_link_libraries(mito_core_smoke_test PRIVATE mito_core)
  target_compile_features(mito_core_smoke_test PRIVATE cxx_std_20)
  mito_enable_warnings(mito_core_smoke_test)
  add_test(NAME mito_core_smoke_test COMMAND mito_core_smoke_test)
endif()

if(MITO_BUILD_BENCHMARKS)
  add_executable(mito_core_benchmark benchmarks/core_benchmark.cpp)
  target_link_libraries(mito_core_benchmark PRIVATE mito_core)
  target_compile_features(mito_core_benchmark PRIVATE cxx_std_20)
  mito_enable_warnings(mito_core_benchmark)
endif()
