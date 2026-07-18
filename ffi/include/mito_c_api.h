#pragma once

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void *mito_engine_new(void);
void mito_engine_delete(void *engine);
bool mito_engine_has_htslib(void);
const char *mito_engine_version(void);
const char *mito_engine_schema_version(void);
const char *mito_engine_error_schema_version(void);
const char *mito_engine_analyze(void *engine, const char *input_path,
                                const char *ref_path);
const char *mito_engine_analyze_with_options(void *engine,
                                             const char *input_path,
                                             const char *ref_path,
                                             bool filter_numt, size_t threads);
const char *mito_engine_analyze_with_cancel(void *engine,
                                            const char *input_path,
                                            const char *ref_path,
                                            bool filter_numt, size_t threads,
                                            bool (*should_cancel)(void *),
                                            void *cancel_user_data);
const char *mito_engine_analyze_with_config(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    bool (*should_cancel)(void *), void *cancel_user_data);
const char *mito_engine_analyze_with_config_v2(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool (*should_cancel)(void *), void *cancel_user_data);
const char *mito_engine_analyze_with_config_v3(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool emit_evidence_graph, size_t max_evidence_observations,
    bool (*should_cancel)(void *), void *cancel_user_data);
const char *mito_engine_analyze_with_config_v4(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool emit_evidence_graph, size_t max_evidence_observations,
    size_t max_phase_links, bool (*should_cancel)(void *),
    void *cancel_user_data);
const char *mito_engine_analyze_with_config_v5(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool emit_evidence_graph, size_t max_evidence_observations,
    size_t max_phase_links, size_t evidence_page_size,
    const char *molecule_id_tag, const char *umi_tag, const char *duplex_tag,
    bool (*should_cancel)(void *), void *cancel_user_data);
const char *mito_engine_get_last_error(void);
const char *mito_engine_get_last_error_code(void);
void mito_engine_free_string(const char *value);

#ifdef __cplusplus
}
#endif
