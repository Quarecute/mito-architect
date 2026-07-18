#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace mito {

/** Version of the scientific result contract requested by a caller. */
enum class ResultSchema : std::uint8_t {
  /** Frozen RC1 aggregate/read projection. */
  v0_5,
  /** Opt-in RC2 evidence graph with explicit fragments, molecules, and
     observations. */
  v0_6,
};

/** Stable failure categories exposed across the C and Rust boundaries. */
enum class AnalysisErrorCode : std::uint8_t {
  invalid_configuration,
  input_open_failed,
  input_format_unsupported,
  input_parse_failed,
  input_empty,
  reference_open_failed,
  reference_invalid,
  resource_open_failed,
  resource_invalid,
  dependency_unavailable,
  analysis_cancelled,
  resource_exhausted,
  internal_error,
};

/**
 * Versioned, machine-readable analysis failure.
 *
 * Messages are diagnostic and may become more specific. Callers must branch on
 * code(), whose string representation is stable within error schema 1.0.
 */
class AnalysisError final : public std::runtime_error {
public:
  AnalysisError(AnalysisErrorCode code, std::string message)
      : std::runtime_error(std::move(message)), code_(code) {}

  [[nodiscard]] AnalysisErrorCode code() const noexcept { return code_; }

private:
  AnalysisErrorCode code_;
};

[[nodiscard]] std::string_view
analysis_error_code_name(AnalysisErrorCode code) noexcept;

/**
 * Runtime configuration for one mitochondrial analysis job.
 *
 * The current implementation supports FASTQ/SAM inputs directly and uses
 * htslib for BAM/CRAM/SAM when the core is built with htslib available.
 * The bundled rCRS reference is used when no explicit reference is supplied.
 */
struct AnalysisConfig {
  /** Result contract to emit. Schema 0.6 is opt-in while the RC2 contract is
   * reviewed. */
  ResultSchema result_schema = ResultSchema::v0_5;

  /** Discard molecules whose NUMT evidence score exceeds numt_threshold. */
  bool filter_numt = true;

  /** NUMT evidence threshold in [0, 1]. Competitive nuclear evidence scores >=
   * 0.90. */
  double numt_threshold = 0.30;

  /** Enable synthetic snp=/sv= read-name controls. Must remain false for real
   * analysis. */
  bool allow_development_tags = false;

  /** Worker count requested by the caller; the core clamps this to at least
   * one. */
  std::size_t threads = 1;

  /** Minimum insertion/deletion length, in bp, to emit as a structural variant.
   */
  std::size_t sv_min_length = 10;

  /** DBSCAN/HDBSCAN-style epsilon over read SNP/SV token distance. */
  double cluster_epsilon = 0.72;

  /** Minimum neighboring reads required to start a non-outlier cluster. */
  std::size_t min_cluster_size = 1;

  /** Minimum alignment MAPQ for SNP calling and locus-callable allele depth. */
  std::uint8_t min_mapping_quality = 20;

  /** Minimum Phred base quality for SNP calling and locus-callable allele
   * depth. */
  std::uint8_t min_base_quality = 10;

  /** Hard in-memory cap for sparse schema-0.6 observations; prevents unbounded
   * JSON jobs. */
  std::size_t max_evidence_observations = 5'000'000;

  /** Hard cap for callable-aware pairwise phase projections in schema 0.6. */
  std::size_t max_phase_links = 1'000'000;

  /** Maximum observations per deterministic schema-0.6 columnar page. */
  std::size_t evidence_page_size = 4096;

  /** Optional SAM tag whose value is an explicit physical-molecule ID (for
   * example MI). Empty preserves the conservative QNAME proxy policy. */
  std::string molecule_id_tag;

  /** Optional SAM tag carrying UMI metadata. It is audited but never used as a
   * molecule key unless molecule_id_tag explicitly names the same tag. */
  std::string umi_tag;

  /** Optional SAM tag carrying a protocol-declared duplex state or ID. */
  std::string duplex_tag;

  /** SAM flags excluded from SNP evidence: secondary, QC-fail, duplicate,
   * supplementary. */
  std::uint16_t excluded_snp_flags = 0xF00;

  /** Optional sample label to put into result metadata. Defaults to input stem.
   */
  std::string sample_name;

  /** Optional versioned nuclear reference identity supplied by an external
   * mapping stage. */
  std::string nuclear_reference_path;

  /**
   * Optional cooperative cancellation hook checked between expensive steps.
   * It may be invoked concurrently by feature-extraction workers and must be
   * thread-safe.
   */
  std::function<bool()> should_cancel;
};

/** Monotonic wall-clock durations for one analysis invocation, in microseconds.
 */
struct AnalysisPhaseTimings {
  std::uint64_t reference_load_us = 0;
  std::uint64_t input_ingest_us = 0;
  std::uint64_t resource_load_us = 0;
  std::uint64_t feature_extraction_us = 0;
  std::uint64_t event_merge_us = 0;
  std::uint64_t clustering_us = 0;
  std::uint64_t evidence_aggregation_us = 0;
  std::uint64_t haplogroup_assignment_us = 0;
  std::uint64_t serialization_us = 0;
  std::uint64_t total_us = 0;
};

/** Analysis payload plus operational timings kept outside the scientific JSON.
 */
struct ProfiledAnalysis {
  std::string json;
  AnalysisPhaseTimings timings;
};

/**
 * mtDNA analysis engine.
 *
 * The engine owns no external state. Each call to analyze() loads the requested
 * input, the rCRS/custom reference, and the local clinical annotation cache,
 * then returns a single JSON document consumed by the CLI, backend, report, and
 * visualization layers.
 */
class AnalysisEngine {
public:
  AnalysisEngine() = default;

  /**
   * Analyze one mtDNA read file.
   *
   * @param input_path FASTQ, SAM, BAM, or CRAM input path. FASTA read input is
   * rejected.
   * @param reference_path Optional FASTA reference path. Empty uses bundled
   * rCRS.
   * @param config Runtime options for filters, SV thresholds, and metadata.
   * @return JSON result containing metadata, filters, coverage, SVs, clusters,
   *         reads, clinical annotations, and protein structure mappings.
   * @throws std::runtime_error if input/reference parsing fails or unsupported
   *         read formats are requested.
   */
  [[nodiscard]] std::string analyze(const std::string &input_path,
                                    const std::string &reference_path,
                                    const AnalysisConfig &config = {}) const;

  /**
   * Analyze one file and return monotonic phase timings for profiling and
   * service observability. Timings are deliberately excluded from the JSON so
   * identical scientific inputs retain byte-identical result payloads.
   */
  [[nodiscard]] ProfiledAnalysis
  analyze_profiled(const std::string &input_path,
                   const std::string &reference_path,
                   const AnalysisConfig &config = {}) const;
};

} // namespace mito
