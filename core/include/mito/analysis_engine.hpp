#pragma once

#include <cstddef>
#include <functional>
#include <string>

namespace mito {

/**
 * Runtime configuration for one mitochondrial analysis job.
 *
 * The current implementation supports FASTQ/SAM inputs directly and uses
 * htslib for BAM/CRAM/SAM when the core is built with htslib available.
 * The bundled rCRS reference is used when no explicit reference is supplied.
 */
struct AnalysisConfig {
  /** Discard reads whose NUMT heuristic score exceeds the configured threshold. */
  bool filter_numt = true;

  /** Worker count requested by the caller; the core clamps this to at least one. */
  std::size_t threads = 1;

  /** Minimum insertion/deletion length, in bp, to emit as a structural variant. */
  std::size_t sv_min_length = 10;

  /** DBSCAN/HDBSCAN-style epsilon over read SNP/SV token distance. */
  double cluster_epsilon = 0.72;

  /** Minimum neighboring reads required to start a non-outlier cluster. */
  std::size_t min_cluster_size = 1;

  /** Optional sample label to put into result metadata. Defaults to input stem. */
  std::string sample_name;

  /** Reserved path for future nuclear remapping against hg38 or another reference. */
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
