#pragma once

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void* mito_engine_new(void);
void mito_engine_delete(void* engine);
const char* mito_engine_analyze(void* engine, const char* input_path, const char* ref_path);
const char* mito_engine_analyze_with_options(void* engine,
                                             const char* input_path,
                                             const char* ref_path,
                                             bool filter_numt,
                                             size_t threads);
const char* mito_engine_analyze_with_cancel(void* engine,
                                            const char* input_path,
                                            const char* ref_path,
                                            bool filter_numt,
                                            size_t threads,
                                            bool (*should_cancel)(void*),
                                            void* cancel_user_data);
const char* mito_engine_get_last_error(void);
void mito_engine_free_string(const char* value);

#ifdef __cplusplus
}
#endif
