#include "mito_c_api.h"

#include "mito/analysis_engine.hpp"

#include <cstdlib>
#include <cstring>
#include <exception>
#include <new>
#include <string>

namespace {

thread_local std::string g_last_error;

[[nodiscard]] const char* copy_to_c_string(const std::string& value) {
  auto* result = static_cast<char*>(std::malloc(value.size() + 1U));
  if (result == nullptr) {
    g_last_error = "malloc failed while returning analysis JSON";
    return nullptr;
  }
  std::memcpy(result, value.c_str(), value.size() + 1U);
  return result;
}

void set_error(const char* message) {
  g_last_error = message == nullptr ? "unknown error" : message;
}

} // namespace

extern "C" void* mito_engine_new(void) {
  try {
    g_last_error.clear();
    return new mito::AnalysisEngine();
  } catch (const std::bad_alloc&) {
    set_error("allocation failed while creating mito engine");
  } catch (const std::exception& error) {
    g_last_error = error.what();
  } catch (...) {
    set_error("unknown error while creating mito engine");
  }
  return nullptr;
}

extern "C" void mito_engine_delete(void* engine) {
  delete static_cast<mito::AnalysisEngine*>(engine);
}

extern "C" const char* mito_engine_analyze(void* engine,
                                           const char* input_path,
                                           const char* ref_path) {
  return mito_engine_analyze_with_options(engine, input_path, ref_path, true, 1);
}

extern "C" const char* mito_engine_analyze_with_options(void* engine,
                                                        const char* input_path,
                                                        const char* ref_path,
                                                        bool filter_numt,
                                                        std::size_t threads) {
  return mito_engine_analyze_with_cancel(engine, input_path, ref_path, filter_numt, threads,
                                         nullptr, nullptr);
}

extern "C" const char* mito_engine_analyze_with_cancel(void* engine,
                                                       const char* input_path,
                                                       const char* ref_path,
                                                       bool filter_numt,
                                                       std::size_t threads,
                                                       bool (*should_cancel)(void*),
                                                       void* cancel_user_data) {
  if (engine == nullptr) {
    set_error("mito_engine_analyze received a null engine");
    return nullptr;
  }
  if (input_path == nullptr) {
    set_error("mito_engine_analyze received a null input path");
    return nullptr;
  }

  try {
    g_last_error.clear();
    mito::AnalysisConfig config;
    config.filter_numt = filter_numt;
    config.threads = threads == 0 ? 1 : threads;
    if (should_cancel != nullptr) {
      config.should_cancel = [should_cancel, cancel_user_data] {
        return should_cancel(cancel_user_data);
      };
    }
    const std::string reference_path = ref_path == nullptr ? std::string{} : std::string(ref_path);
    const auto json =
        static_cast<mito::AnalysisEngine*>(engine)->analyze(input_path, reference_path, config);
    return copy_to_c_string(json);
  } catch (const std::exception& error) {
    g_last_error = error.what();
  } catch (...) {
    set_error("unknown error during mtDNA analysis");
  }
  return nullptr;
}

extern "C" const char* mito_engine_get_last_error(void) {
  return g_last_error.c_str();
}

extern "C" void mito_engine_free_string(const char* value) {
  std::free(const_cast<char*>(value));
}
