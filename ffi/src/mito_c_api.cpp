#include "mito_c_api.h"

#include "mito/analysis_engine.hpp"
#include "mito/version.hpp"

#include <cstdlib>
#include <cstring>
#include <exception>
#include <new>
#include <string>

namespace {

thread_local std::string g_last_error;
thread_local std::string g_last_error_code;

void clear_error() {
  g_last_error.clear();
  g_last_error_code.clear();
}

void set_error(const char* code, const char* message) {
  g_last_error_code = code == nullptr ? "MITO-E9001" : code;
  g_last_error = message == nullptr ? "unknown error" : message;
}

[[nodiscard]] const char* copy_to_c_string(const std::string& value) {
  auto* result = static_cast<char*>(std::malloc(value.size() + 1U));
  if (result == nullptr) {
    set_error("MITO-E1601", "malloc failed while returning analysis JSON");
    return nullptr;
  }
  std::memcpy(result, value.c_str(), value.size() + 1U);
  return result;
}

} // namespace

extern "C" void* mito_engine_new(void) {
  try {
    clear_error();
    return new mito::AnalysisEngine();
  } catch (const std::bad_alloc&) {
    set_error("MITO-E1601", "allocation failed while creating mito engine");
  } catch (const std::exception& error) {
    set_error("MITO-E9001", error.what());
  } catch (...) {
    set_error("MITO-E9001", "unknown error while creating mito engine");
  }
  return nullptr;
}

extern "C" void mito_engine_delete(void* engine) {
  delete static_cast<mito::AnalysisEngine*>(engine);
}

extern "C" bool mito_engine_has_htslib(void) {
#ifdef MITO_HAS_HTSLIB
  return true;
#else
  return false;
#endif
}

extern "C" const char* mito_engine_version(void) {
  return mito::kEngineVersion;
}

extern "C" const char* mito_engine_schema_version(void) {
  return mito::kResultSchemaVersion;
}

extern "C" const char* mito_engine_error_schema_version(void) {
  return mito::kErrorSchemaVersion;
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
  return mito_engine_analyze_with_config(engine, input_path, ref_path, filter_numt, threads, 20, 10,
                                         0xF00, should_cancel, cancel_user_data);
}

extern "C" const char* mito_engine_analyze_with_config(void* engine,
                                                        const char* input_path,
                                                        const char* ref_path,
                                                        bool filter_numt,
                                                        std::size_t threads,
                                                        unsigned char min_mapping_quality,
                                                        unsigned char min_base_quality,
                                                        unsigned short excluded_snp_flags,
                                                        bool (*should_cancel)(void*),
                                                        void* cancel_user_data) {
  return mito_engine_analyze_with_config_v2(
      engine, input_path, ref_path, filter_numt, threads, min_mapping_quality, min_base_quality,
      excluded_snp_flags, 0.30, false, should_cancel, cancel_user_data);
}

extern "C" const char* mito_engine_analyze_with_config_v2(
    void* engine,
    const char* input_path,
    const char* ref_path,
    bool filter_numt,
    std::size_t threads,
    unsigned char min_mapping_quality,
    unsigned char min_base_quality,
    unsigned short excluded_snp_flags,
    double numt_threshold,
    bool allow_development_tags,
    bool (*should_cancel)(void*),
    void* cancel_user_data) {
  if (engine == nullptr) {
    set_error("MITO-E1001", "mito_engine_analyze received a null engine");
    return nullptr;
  }
  if (input_path == nullptr) {
    set_error("MITO-E1001", "mito_engine_analyze received a null input path");
    return nullptr;
  }

  try {
    clear_error();
    mito::AnalysisConfig config;
    config.filter_numt = filter_numt;
    config.threads = threads == 0 ? 1 : threads;
    config.min_mapping_quality = min_mapping_quality;
    config.min_base_quality = min_base_quality;
    config.excluded_snp_flags = excluded_snp_flags;
    config.numt_threshold = numt_threshold;
    config.allow_development_tags = allow_development_tags;
    if (should_cancel != nullptr) {
      config.should_cancel = [should_cancel, cancel_user_data] {
        return should_cancel(cancel_user_data);
      };
    }
    const std::string reference_path = ref_path == nullptr ? std::string{} : std::string(ref_path);
    const auto json =
        static_cast<mito::AnalysisEngine*>(engine)->analyze(input_path, reference_path, config);
    return copy_to_c_string(json);
  } catch (const mito::AnalysisError& error) {
    const auto code = mito::analysis_error_code_name(error.code());
    g_last_error_code.assign(code.data(), code.size());
    g_last_error = error.what();
  } catch (const std::bad_alloc&) {
    set_error("MITO-E1601", "allocation failed during mtDNA analysis");
  } catch (const std::exception& error) {
    set_error("MITO-E9001", error.what());
  } catch (...) {
    set_error("MITO-E9001", "unknown error during mtDNA analysis");
  }
  return nullptr;
}

extern "C" const char* mito_engine_get_last_error(void) {
  return g_last_error.c_str();
}

extern "C" const char* mito_engine_get_last_error_code(void) {
  return g_last_error_code.empty() ? "MITO-E9001" : g_last_error_code.c_str();
}

extern "C" void mito_engine_free_string(const char* value) {
  std::free(const_cast<char*>(value));
}
