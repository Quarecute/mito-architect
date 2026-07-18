#include "mito/analysis_engine.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

namespace {

std::size_t read_count_from_env() {
  if (const char *value = std::getenv("MITO_BENCH_READS")) {
    try {
      return std::max<std::size_t>(1, std::stoull(value));
    } catch (...) {
    }
  }
  return 10000;
}

mito::ResultSchema result_schema_from_env() {
  if (const char *value = std::getenv("MITO_BENCH_SCHEMA")) {
    const std::string schema(value);
    if (schema == "0.6") {
      return mito::ResultSchema::v0_6;
    }
    if (schema != "0.5") {
      std::cerr << "MITO_BENCH_SCHEMA must be 0.5 or 0.6\n";
      std::exit(2);
    }
  }
  return mito::ResultSchema::v0_5;
}

std::size_t peak_rss_kib() {
#if defined(__unix__) || defined(__APPLE__)
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0 || usage.ru_maxrss < 0) {
    std::cerr << "getrusage(RUSAGE_SELF) failed\n";
    std::exit(2);
  }
#if defined(__APPLE__)
  return static_cast<std::size_t>(usage.ru_maxrss) / 1024U;
#else
  return static_cast<std::size_t>(usage.ru_maxrss);
#endif
#else
  std::cerr << "peak RSS measurement is unsupported on this platform\n";
  std::exit(2);
#endif
}

std::filesystem::path write_fixture(std::size_t reads) {
  const auto path =
      std::filesystem::temp_directory_path() / "mito_core_benchmark.sam";
  std::ofstream out(path);
  out << "@SQ\tSN:MT\tLN:16569\n";
  for (std::size_t i = 0; i < reads; ++i) {
    const std::size_t start = 1 + (i % 16000);
    out << "bench-" << i << "\t0\tMT\t" << start
        << "\t60\t48M\t*"
           "\t0\t0\tGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATAC\t"
        << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
  }
  return path;
}

} // namespace

int main() {
  const auto reads = read_count_from_env();
  const auto path = write_fixture(reads);

  mito::AnalysisConfig config;
  config.sample_name = "benchmark";
  config.threads = std::thread::hardware_concurrency();
  config.result_schema = result_schema_from_env();

  const auto analysis =
      mito::AnalysisEngine().analyze_profiled(path.string(), "", config);
  const auto &timings = analysis.timings;
  std::cout << "reads=" << reads << "\n";
  std::cout << "schema="
            << (config.result_schema == mito::ResultSchema::v0_6 ? "0.6"
                                                                 : "0.5")
            << "\n";
  std::cout << "elapsed_ms=" << timings.total_us / 1000U << "\n";
  std::cout << "json_bytes=" << analysis.json.size() << "\n";
  std::cout << "peak_rss_kib=" << peak_rss_kib() << "\n";
  std::cout << "reference_load_us=" << timings.reference_load_us << "\n";
  std::cout << "input_ingest_us=" << timings.input_ingest_us << "\n";
  std::cout << "resource_load_us=" << timings.resource_load_us << "\n";
  std::cout << "feature_extraction_us=" << timings.feature_extraction_us
            << "\n";
  std::cout << "event_merge_us=" << timings.event_merge_us << "\n";
  std::cout << "clustering_us=" << timings.clustering_us << "\n";
  std::cout << "evidence_aggregation_us=" << timings.evidence_aggregation_us
            << "\n";
  std::cout << "haplogroup_assignment_us=" << timings.haplogroup_assignment_us
            << "\n";
  std::cout << "serialization_us=" << timings.serialization_us << "\n";
  std::cout << "total_us=" << timings.total_us << "\n";

  std::filesystem::remove(path);
  return 0;
}
