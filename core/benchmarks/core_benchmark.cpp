#include "mito/analysis_engine.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>

namespace {

std::size_t read_count_from_env() {
  if (const char* value = std::getenv("MITO_BENCH_READS")) {
    try {
      return std::max<std::size_t>(1, std::stoull(value));
    } catch (...) {
    }
  }
  return 10000;
}

std::filesystem::path write_fixture(std::size_t reads) {
  const auto path = std::filesystem::temp_directory_path() / "mito_core_benchmark.sam";
  std::ofstream out(path);
  out << "@SQ\tSN:MT\tLN:16569\n";
  for (std::size_t i = 0; i < reads; ++i) {
    const std::size_t start = 1 + (i % 16000);
    out << "bench-" << i << "\t0\tMT\t" << start
        << "\t60\t48M\t*\t0\t0\tGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATAC\t"
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

  const auto start = std::chrono::steady_clock::now();
  const auto json = mito::AnalysisEngine().analyze(path.string(), "", config);
  const auto end = std::chrono::steady_clock::now();

  const auto elapsed_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "reads=" << reads << "\n";
  std::cout << "elapsed_ms=" << elapsed_ms << "\n";
  std::cout << "json_bytes=" << json.size() << "\n";

  std::filesystem::remove(path);
  return 0;
}
