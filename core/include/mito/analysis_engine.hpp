#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace mito {

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

[[nodiscard]] std::string_view analysis_error_code_name(AnalysisErrorCode code) noexcept;

/**
 * Runtime configuration for one mitochondrial analysis job.
 *
 * The current implementation supports FASTQ/SAM inputs directly and uses
 * htslib for BAM/CRAM/SAM when the core is built with htslib available.
 * The bundled rCRS reference is used when no explicit reference is supplied.
 */
struct AnalysisConfig {
  /** Discard molecules whose NUMT evidence score exceeds numt_threshold. */
  bool filter_numt = true;

  /** NUMT evidence threshold in [0, 1]. Competitive nuclear evidence scores >= 0.90. */
  double numt_threshold = 0.30;

  /** Enable synthetic snp=/sv= read-name controls. Must remain false for real analysis. */
  bool allow_development_tags = false;

  /** Worker count requested by the caller; the core clamps this to at least one. */
  std::size_t threads = 1;

  /** Minimum insertion/deletion length, in bp, to emit as a structural variant. */
  std::size_t sv_min_length = 10;

  /** DBSCAN/HDBSCAN-style epsilon over read SNP/SV token distance. */
  double cluster_epsilon = 0.72;

  /** Minimum neighboring reads required to start a non-outlier cluster. */
  std::size_t min_cluster_size = 1;

  /** Minimum alignment MAPQ for SNP calling and locus-callable allele depth. */
  std::uint8_t min_mapping_quality = 20;

  /** Minimum Phred base quality for SNP calling and locus-callable allele depth. */
  std::uint8_t min_base_quality = 10;

  /** SAM flags excluded from SNP evidence: secondary, QC-fail, duplicate, supplementary. */
  std::uint16_t excluded_snp_flags = 0xF00;

  /** Optional sample label to put into result metadata. Defaults to input stem. */
  std::string sample_name;

  /** Optional versioned nuclear reference identity supplied by an external mapping stage. */
  std::string nuclear_reference_path;

  /**
   * Optional cooperative cancellation hook checked between expensive steps.
   * It may be invoked concurrently by feature-extraction workers and must be thread-safe.
   */
  std::function<bool()> should_cancel;
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
   * @param input_path FASTQ, SAM, BAM, or CRAM input path. FASTA read input is rejected.
   * @param reference_path Optional FASTA reference path. Empty uses bundled rCRS.
   * @param config Runtime options for filters, SV thresholds, and metadata.
   * @return JSON result containing metadata, filters, coverage, SVs, clusters,
   *         reads, clinical annotations, and protein structure mappings.
   * @throws std::runtime_error if input/reference parsing fails or unsupported
   *         read formats are requested.
   */
  [[nodiscard]] std::string analyze(const std::string& input_path,
                                    const std::string& reference_path,
                                    const AnalysisConfig& config = {}) const;
};

} // namespace mito
