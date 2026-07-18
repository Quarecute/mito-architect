#include "mito/analysis_engine.hpp"

#include "mito/version.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <charconv>
#include <chrono>
#include <cmath>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <new>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef MITO_HAS_HTSLIB
#include <htslib/sam.h>
#endif

#ifdef MITO_HAS_HDBSCAN_CPP
#include <Hdbscan/hdbscan.hpp>
#endif

namespace mito {

std::string_view
analysis_error_code_name(const AnalysisErrorCode code) noexcept {
  switch (code) {
  case AnalysisErrorCode::invalid_configuration:
    return "MITO-E1001";
  case AnalysisErrorCode::input_open_failed:
    return "MITO-E1101";
  case AnalysisErrorCode::input_format_unsupported:
    return "MITO-E1102";
  case AnalysisErrorCode::input_parse_failed:
    return "MITO-E1103";
  case AnalysisErrorCode::input_empty:
    return "MITO-E1104";
  case AnalysisErrorCode::reference_open_failed:
    return "MITO-E1201";
  case AnalysisErrorCode::reference_invalid:
    return "MITO-E1202";
  case AnalysisErrorCode::resource_open_failed:
    return "MITO-E1301";
  case AnalysisErrorCode::resource_invalid:
    return "MITO-E1302";
  case AnalysisErrorCode::dependency_unavailable:
    return "MITO-E1401";
  case AnalysisErrorCode::analysis_cancelled:
    return "MITO-E1501";
  case AnalysisErrorCode::resource_exhausted:
    return "MITO-E1601";
  case AnalysisErrorCode::internal_error:
    return "MITO-E9001";
  }
  return "MITO-E9001";
}

namespace {

struct GeneAnnotation {
  std::string name;
  std::size_t start;
  std::size_t end;
  std::string strand;
  std::string biotype;
};

struct CigarOperation {
  std::size_t length = 0;
  char code = 'M';
};

struct SupplementaryAlignment {
  std::string reference_name;
  std::size_t reference_start = 0;
  char strand = '+';
  std::string cigar;
  std::vector<CigarOperation> cigar_operations;
  std::uint8_t mapping_quality = 0;
  std::size_t edit_distance = 0;
};

struct ReadRecord {
  std::string id;
  std::string sequence;
  std::string qualities;
  std::size_t reference_start = 0;
  std::string reference_name;
  std::string cigar;
  std::vector<CigarOperation> cigar_operations;
  std::uint16_t flags = 0;
  std::uint8_t mapping_quality = 0;
  std::map<std::string, std::string> aux_tags;
  std::vector<SupplementaryAlignment> supplementary_alignments;
};

struct InputRecords {
  std::vector<ReadRecord> reads;
  bool alignment_input = false;
  bool has_nuclear_contigs = false;
};

struct AlignmentFragmentId {
  std::size_t value = 0;

  auto operator<=>(const AlignmentFragmentId &) const = default;
};

struct MoleculeIndex {
  std::size_t value = 0;

  auto operator<=>(const MoleculeIndex &) const = default;
};

struct AlignmentFragmentEvidence {
  AlignmentFragmentId id;
  MoleculeIndex molecule_index;
  std::size_t source_record_index = 0;
  std::string molecule_id;
  std::string role;
  bool selected_representative = false;
  std::uint16_t flags = 0;
  std::uint8_t mapping_quality = 0;
  std::string reference_name;
  std::size_t reference_start = 0;
  std::string cigar;
};

struct MoleculeAssemblyEvidence {
  MoleculeIndex index;
  std::string id;
  std::string identity_policy;
  std::string assembly_status;
  std::vector<AlignmentFragmentId> fragment_ids;
  AlignmentFragmentId representative_fragment_id;
  std::size_t primary_candidate_count = 0;
  bool ambiguous = false;
  bool analysis_eligible = true;
  std::vector<std::string> source_qnames;
  std::map<std::string, std::string> protocol_metadata;
  std::vector<std::string> protocol_flags;
  std::vector<std::string> exclusion_reasons;
  std::vector<std::string> warnings;
};

struct MoleculeAssemblyResult {
  std::vector<ReadRecord> representatives;
  /** Original alignment records retained for provenance and future
   * multi-fragment evidence. */
  std::vector<ReadRecord> source_alignment_records;
  std::vector<AlignmentFragmentEvidence> fragments;
  std::vector<MoleculeAssemblyEvidence> molecules;
};

struct ClinicalAssertion {
  std::string source;
  std::string assertion_id;
  std::string allele_id;
  std::string disease;
  std::string clinical_significance;
  std::string normalized_significance;
  std::string review_status;
  std::string assertion_date;
  std::string source_url;
  std::vector<std::string> references;
  std::string resource_version;
  std::string retrieved_at;

  auto operator<=>(const ClinicalAssertion &) const = default;
};

struct SnpCall {
  std::size_t position = 0;
  char reference = 'N';
  char alternate = 'N';
  std::string gene;
  std::string consequence;
  std::string protein;
  std::string residue;
  std::string phenotype;
  std::string pathogenicity;
  std::vector<std::string> references;
  std::vector<std::string> sources;
  std::string structure_id;
  std::string structure_chain;
  std::size_t structure_residue = 0;
  std::string structure_complex;
  std::string clinvar_allele_id;
  std::string mitomap_url;
  std::string clinical_conflict_status;
  std::string clinical_consensus_significance;
  std::vector<ClinicalAssertion> clinical_assertions;
};

enum class PhyloMutationKind : std::uint8_t {
  substitution,
  deletion,
  insertion,
};

struct PhyloMutationLocus {
  std::size_t position = 0;
  PhyloMutationKind kind = PhyloMutationKind::substitution;
  std::string insertion_index;

  auto operator<=>(const PhyloMutationLocus &) const = default;
};

struct PhyloMutation {
  PhyloMutationLocus locus;
  PhyloMutationKind kind = PhyloMutationKind::substitution;
  char alternate = 'N';
  std::string inserted_bases;
  std::string encoded;
  double weight = 1.0;
  bool backmutation = false;
};

struct ObservedPhyloMutation {
  std::size_t position = 0;
  PhyloMutationKind kind = PhyloMutationKind::substitution;
  char alternate = 'N';
  std::string insertion_index;
  std::string inserted_bases;
  std::string encoded;
};

struct PhyloAlignmentRule {
  std::vector<std::string> erroneous_markers;
  std::vector<ObservedPhyloMutation> replacement_markers;
};

struct SvCall {
  std::string id;
  std::string type;
  std::size_t start = 0;
  std::size_t end = 0;
  std::size_t length = 0;
  std::vector<std::string> supporting_reads;
  bool known_event = false;
  std::vector<std::string> evidence_sources;
  std::vector<std::string> orientations;
  std::size_t segment_count = 1;
};

struct ComplexSvCall {
  std::string id;
  std::vector<std::string> junction_ids;
  std::vector<std::string> junction_orientations;
  std::vector<std::string> supporting_reads;
  std::size_t segment_count = 0;
};

struct ReadFeature {
  std::string id;
  std::size_t length = 0;
  double mean_quality = 0.0;
  double numt_score = 0.0;
  bool filtered_numt = false;
  std::vector<std::string> numt_evidence;
  int cluster_id = -1;
  bool outlier = false;
  std::uint8_t mapping_quality = 0;
  std::uint16_t flags = 0;
  std::string reference_name;
  std::map<std::string, std::string> aux_tags;
  std::vector<SnpCall> snps;
  std::vector<ObservedPhyloMutation> haplogroup_markers;
  std::vector<std::pair<std::size_t, std::size_t>> haplogroup_ranges;
  bool haplogroup_range_known = false;
  std::vector<std::string> sv_ids;
  std::vector<std::string> complex_event_ids;
};

struct CoverageBin {
  std::size_t start = 0;
  std::size_t end = 0;
  std::size_t depth = 0;
};

struct CoverageResult {
  std::vector<CoverageBin> bins;
  double mean_depth = 0.0;
  double pct_sites_gt20x = 0.0;
  std::size_t max_depth = 0;
};

struct FeatureExtractionResult {
  ReadFeature feature;
  std::vector<SvCall> svs;
  std::vector<ComplexSvCall> complex_events;
};

struct SnpAggregate {
  SnpCall call;
  std::size_t alternate_depth = 0;
  std::size_t reference_depth = 0;
  std::size_t other_depth = 0;
  std::size_t callable_depth = 0;
  double heteroplasmy = 0.0;
  double ci95_low = 0.0;
  double ci95_high = 0.0;
  std::vector<std::string> supporting_reads;
  std::array<std::size_t, 2> alternate_strand_depths{};
  std::array<std::size_t, 2> reference_strand_depths{};
  std::array<std::size_t, 2> other_strand_depths{};
  double alternate_quality_sum = 0.0;
  double reference_quality_sum = 0.0;
  double other_quality_sum = 0.0;
  std::uint8_t alternate_quality_min = std::numeric_limits<std::uint8_t>::max();
  std::uint8_t reference_quality_min = std::numeric_limits<std::uint8_t>::max();
  std::uint8_t other_quality_min = std::numeric_limits<std::uint8_t>::max();
  std::uint8_t alternate_quality_max = 0;
  std::uint8_t reference_quality_max = 0;
  std::uint8_t other_quality_max = 0;
  double alternate_read_position_sum = 0.0;
  double reference_read_position_sum = 0.0;
  double other_read_position_sum = 0.0;
};

struct LocusAlleleAccumulator {
  std::array<std::size_t, 4> depths{};
  std::array<std::array<std::size_t, 2>, 4> strand_depths{};
  std::array<double, 4> quality_sums{};
  std::array<std::uint8_t, 4> quality_mins{};
  std::array<std::uint8_t, 4> quality_maxes{};
  std::array<double, 4> read_position_sums{};

  LocusAlleleAccumulator() {
    quality_mins.fill(std::numeric_limits<std::uint8_t>::max());
  }
};

struct MoleculeAlleleObservation {
  char base = 'N';
  std::uint8_t quality = 0;
  double center_proximity = 0.0;
  std::size_t query_index = 0;
  AlignmentFragmentId alignment_fragment_id;
  std::uint8_t mapping_quality = 0;
  char strand = '+';
  bool has_passing_observation = false;
  bool covered = false;
  bool conflicted = false;
};

enum class ObservationState : std::uint8_t {
  reference,
  alternate,
  event_absent,
  not_callable,
  low_quality,
  conflict,
};

struct EvidenceEvent {
  std::string id;
  std::string type;
  std::size_t start = 0;
  std::size_t end = 0;
  std::size_t length = 0;
  std::string reference;
  std::string alternate;
  std::string normalization;
  std::string source_projection;
  std::string negative_evidence_rule;
  bool absence_assessable = false;
  std::vector<std::string> component_event_ids;
  std::vector<std::string> supporting_molecules;
};

struct EvidenceObservation {
  std::size_t id = 0;
  MoleculeIndex molecule_index;
  std::size_t event_index = 0;
  AlignmentFragmentId alignment_fragment_id;
  ObservationState state = ObservationState::not_callable;
  std::string observed_allele;
  std::optional<std::uint8_t> base_quality;
  std::uint8_t mapping_quality = 0;
  char strand = '+';
  std::optional<double> center_proximity;
  std::string evidence_source;
};

struct AlignmentCallabilityEvidence {
  AlignmentFragmentId alignment_fragment_id;
  bool eligible = false;
  std::string status;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  std::size_t callable_bases = 0;
  std::map<std::string, std::size_t> reference_exclusion_counts;
  std::vector<std::size_t> disrupted_adjacency_anchors;
  std::size_t inserted_query_bases = 0;
  std::size_t soft_clipped_query_bases = 0;
};

struct MoleculeCallabilityEvidence {
  MoleculeIndex molecule_index;
  bool known = false;
  std::string status;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  std::size_t callable_bases = 0;
  std::vector<AlignmentCallabilityEvidence> alignments;
};

struct CallabilityResult {
  std::vector<MoleculeCallabilityEvidence> molecules;
};

struct PhaseLink {
  std::size_t event_a_index = 0;
  std::size_t event_b_index = 0;
  std::size_t jointly_callable = 0;
  std::size_t both_alternate = 0;
  std::size_t a_alternate_b_absent = 0;
  std::size_t a_absent_b_alternate = 0;
  std::size_t neither_alternate = 0;
  std::size_t jointly_uncertain = 0;
  double co_alternate_fraction = 0.0;
  double ci95_low = 0.0;
  double ci95_high = 0.0;
  double expected_co_alternate_fraction = 0.0;
  double linkage_delta = 0.0;
  bool complete_callability = false;
  std::vector<std::size_t> supporting_molecule_indices;
  std::vector<std::size_t> uncertain_molecule_indices;
};

struct NormalizedSmallIndel {
  std::string id;
  std::string type;
  std::size_t start = 0;
  std::size_t end = 0;
  std::size_t length = 0;
  std::string reference;
  std::string alternate;
  std::string normalization;
  std::optional<std::uint8_t> base_quality;
};

struct SparseEvidenceStore {
  std::vector<EvidenceEvent> events;
  std::vector<EvidenceObservation> observations;
};

struct EvidenceAggregationResult {
  std::map<std::string, SnpAggregate> variants;
  SparseEvidenceStore store;
};

struct VariantAlleleEvidenceSummary {
  std::size_t count = 0U;
  std::array<std::size_t, 2> strand_depths{};
  double base_quality_sum = 0.0;
  std::size_t base_quality_count = 0U;
  std::uint8_t base_quality_min = std::numeric_limits<std::uint8_t>::max();
  std::uint8_t base_quality_max = 0U;
  double mapping_quality_sum = 0.0;
  std::uint8_t mapping_quality_min = std::numeric_limits<std::uint8_t>::max();
  std::uint8_t mapping_quality_max = 0U;
  double read_position_sum = 0.0;
  std::size_t read_position_count = 0U;
};

struct UnifiedVariantEvidenceSummary {
  VariantAlleleEvidenceSummary alternate;
  VariantAlleleEvidenceSummary reference;
  VariantAlleleEvidenceSummary other;
  std::size_t event_absent = 0U;
  std::size_t low_quality = 0U;
  std::size_t conflict = 0U;
  bool multi_allelic = false;
};

struct ResourceRecord {
  std::string name;
  std::string version;
  std::string path;
  std::string sha256;
  std::string source;
  std::string license;
  std::string retrieved;
};

struct HaplogroupDefinition {
  std::string name;
  std::vector<PhyloMutation> mutations;
  std::size_t source_order = 0U;
};

struct HaplogroupCandidate {
  std::string name;
  double score = 0.0;
  std::vector<std::string> matched;
  std::vector<std::string> missing;
  std::vector<std::string> extra;
};

struct ClusterHaplogroupAssignment {
  std::string best = "unassigned";
  double quality = 0.0;
  bool contamination_warning = false;
  std::vector<std::string> observed_markers;
  std::vector<std::pair<std::size_t, std::size_t>> callable_ranges;
  std::vector<HaplogroupCandidate> candidates;
};

[[nodiscard]] bool looks_like_nuclear_contig(std::string_view reference_name);

[[nodiscard]] std::string escape_json(std::string_view value) {
  std::string out;
  out.reserve(value.size() + 8);
  for (const char c : value) {
    switch (c) {
    case '"':
      out += "\\\"";
      break;
    case '\\':
      out += "\\\\";
      break;
    case '\b':
      out += "\\b";
      break;
    case '\f':
      out += "\\f";
      break;
    case '\n':
      out += "\\n";
      break;
    case '\r':
      out += "\\r";
      break;
    case '\t':
      out += "\\t";
      break;
    default:
      if (static_cast<unsigned char>(c) < 0x20U) {
        std::ostringstream oss;
        oss << "\\u" << std::hex << std::setw(4) << std::setfill('0')
            << static_cast<int>(c);
        out += oss.str();
      } else {
        out += c;
      }
    }
  }
  return out;
}

[[nodiscard]] std::string quoted(std::string_view value) {
  return "\"" + escape_json(value) + "\"";
}

// Prefer the JSON escaper above over std::quoted when argument-dependent lookup
// sees std::string. The latter only escapes quotes/backslashes and is not a
// complete JSON control-character encoder.
[[nodiscard]] std::string quoted(const std::string &value) {
  return quoted(std::string_view(value));
}

[[nodiscard]] std::string quoted(std::string &value) {
  return quoted(std::string_view(value));
}

[[nodiscard]] std::string quoted(const char *value) {
  return quoted(std::string_view(value));
}

[[nodiscard]] std::vector<std::string> split_tab(std::string_view line) {
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (start <= line.size()) {
    const std::size_t next = line.find('\t', start);
    if (next == std::string_view::npos) {
      fields.emplace_back(line.substr(start));
      break;
    }
    fields.emplace_back(line.substr(start, next - start));
    start = next + 1;
  }
  return fields;
}

[[nodiscard]] std::string trim_copy(std::string_view value) {
  std::size_t first = 0;
  while (first < value.size() &&
         std::isspace(static_cast<unsigned char>(value[first])) != 0) {
    ++first;
  }
  std::size_t last = value.size();
  while (last > first &&
         std::isspace(static_cast<unsigned char>(value[last - 1])) != 0) {
    --last;
  }
  return std::string(value.substr(first, last - first));
}

void strip_trailing_carriage_return(std::string &value) {
  if (!value.empty() && value.back() == '\r') {
    value.pop_back();
  }
}

[[nodiscard]] std::vector<std::string> split_semicolon(std::string_view line) {
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (start <= line.size()) {
    const std::size_t next = line.find(';', start);
    const auto token = trim_copy(
        line.substr(start, next == std::string_view::npos ? line.size() - start
                                                          : next - start));
    if (!token.empty()) {
      fields.push_back(token);
    }
    if (next == std::string_view::npos) {
      break;
    }
    start = next + 1;
  }
  return fields;
}

[[nodiscard]] std::string parse_reference_fasta(std::istream &in,
                                                const std::string &source) {
  std::string reference;
  std::string line;
  std::size_t header_count = 0;
  while (std::getline(in, line)) {
    if (!line.empty() && line.front() == '>') {
      ++header_count;
      if (header_count > 1) {
        throw AnalysisError(
            AnalysisErrorCode::reference_invalid,
            "reference FASTA must contain exactly one sequence: " + source);
      }
      continue;
    }
    for (const char base : line) {
      if (std::isspace(static_cast<unsigned char>(base)) != 0) {
        continue;
      }
      const char upper =
          static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
      if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T' ||
          upper == 'N') {
        reference.push_back(upper);
      } else {
        throw AnalysisError(AnalysisErrorCode::reference_invalid,
                            "reference FASTA contains an unsupported base: " +
                                source);
      }
    }
  }

  if (header_count == 0) {
    throw AnalysisError(AnalysisErrorCode::reference_invalid,
                        "reference FASTA is missing a header: " + source);
  }
  if (reference.empty()) {
    throw AnalysisError(AnalysisErrorCode::reference_invalid,
                        "reference FASTA is empty: " + source);
  }
  return reference;
}

[[nodiscard]] std::string bundled_rcrs_path() {
#ifdef MITO_DEFAULT_RCRS_FASTA_PATH
  return MITO_DEFAULT_RCRS_FASTA_PATH;
#else
  return "core/data/rcrs.fasta";
#endif
}

[[nodiscard]] std::string bundled_clinical_annotations_path() {
#ifdef MITO_CLINICAL_ANNOTATIONS_PATH
  return MITO_CLINICAL_ANNOTATIONS_PATH;
#else
  return "core/data/clinical_annotations.tsv";
#endif
}

[[nodiscard]] std::string clinical_annotations_path() {
  if (const char *override_path = std::getenv("MITO_CLINICAL_ANNOTATIONS")) {
    if (override_path[0] != '\0') {
      return override_path;
    }
  }
  return bundled_clinical_annotations_path();
}

[[nodiscard]] std::string phylotree_path() {
  if (const char *override_path = std::getenv("MITO_PHYLOTREE")) {
    if (override_path[0] != '\0') {
      return override_path;
    }
  }
#ifdef MITO_PHYLOTREE_PATH
  return MITO_PHYLOTREE_PATH;
#else
  return "core/data/phylotree-rcrs-17.3.xml";
#endif
}

[[nodiscard]] std::string phylotree_weights_path() {
  if (const char *override_path = std::getenv("MITO_PHYLOTREE_WEIGHTS")) {
    if (override_path[0] != '\0') {
      return override_path;
    }
  }
#ifdef MITO_PHYLOTREE_WEIGHTS_PATH
  return MITO_PHYLOTREE_WEIGHTS_PATH;
#else
  return "core/data/phylotree-rcrs-17.3-weights.txt";
#endif
}

[[nodiscard]] std::string phylotree_alignment_rules_path() {
  if (const char *override_path =
          std::getenv("MITO_PHYLOTREE_ALIGNMENT_RULES")) {
    if (override_path[0] != '\0') {
      return override_path;
    }
  }
#ifdef MITO_PHYLOTREE_ALIGNMENT_RULES_PATH
  return MITO_PHYLOTREE_ALIGNMENT_RULES_PATH;
#else
  return "core/data/phylotree-rcrs-17.3-rules.csv";
#endif
}

[[nodiscard]] std::string resource_manifest_path() {
#ifdef MITO_RESOURCE_MANIFEST_PATH
  return MITO_RESOURCE_MANIFEST_PATH;
#else
  return "core/data/resource_manifest.tsv";
#endif
}

[[nodiscard]] std::vector<ResourceRecord> load_resource_manifest() {
  std::vector<ResourceRecord> records;
  std::ifstream input(resource_manifest_path());
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open resource manifest: " +
                            resource_manifest_path());
  }
  std::string line;
  bool header = true;
  while (std::getline(input, line)) {
    strip_trailing_carriage_return(line);
    if (line.empty() || header) {
      header = false;
      continue;
    }
    const auto fields = split_tab(line);
    if (fields.size() != 7) {
      throw AnalysisError(AnalysisErrorCode::resource_invalid,
                          "resource manifest row must contain seven fields");
    }
    records.push_back({fields[0], fields[1], fields[2], fields[3], fields[4],
                       fields[5], fields[6]});
  }
  if (records.empty()) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "resource manifest contains no resources");
  }
  return records;
}

void apply_runtime_resource_overrides(std::vector<ResourceRecord> &records) {
  const auto clinical_path = clinical_annotations_path();
  if (clinical_path == bundled_clinical_annotations_path()) {
    return;
  }
  const auto clinical =
      std::find_if(records.begin(), records.end(), [](const auto &resource) {
        return resource.name == "clinical-curated";
      });
  if (clinical == records.end()) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "resource manifest is missing clinical-curated");
  }
  clinical->version = "external-unpinned";
  clinical->path = clinical_path;
  clinical->sha256 = "not-verified";
  clinical->source = "MITO_CLINICAL_ANNOTATIONS override";
  clinical->license = "not-recorded";
  clinical->retrieved = "not-recorded";
}

void throw_if_cancelled(const AnalysisConfig &config) {
  if (config.should_cancel && config.should_cancel()) {
    throw AnalysisError(AnalysisErrorCode::analysis_cancelled,
                        "analysis cancelled");
  }
}

[[nodiscard]] bool valid_sam_tag_name(const std::string &tag) {
  return tag.empty() ||
         (tag.size() == 2U &&
          std::isalpha(static_cast<unsigned char>(tag[0])) != 0 &&
          std::isalnum(static_cast<unsigned char>(tag[1])) != 0);
}

void validate_config(const AnalysisConfig &config) {
  if (!std::isfinite(config.cluster_epsilon) || config.cluster_epsilon < 0.0 ||
      config.cluster_epsilon > 1.0) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "cluster_epsilon must be finite and between 0 and 1");
  }
  if (!std::isfinite(config.numt_threshold) || config.numt_threshold < 0.0 ||
      config.numt_threshold > 1.0) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "numt_threshold must be finite and between 0 and 1");
  }
  if (config.min_cluster_size == 0) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "min_cluster_size must be at least 1");
  }
  if (config.sv_min_length == 0) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "sv_min_length must be at least 1");
  }
  if (config.min_base_quality > 93U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "min_base_quality must not exceed Phred 93");
  }
  if (config.max_evidence_observations == 0U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "max_evidence_observations must be at least 1");
  }
  if (config.max_phase_links == 0U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "max_phase_links must be at least 1");
  }
  if (config.evidence_page_size == 0U ||
      config.evidence_page_size > 1'000'000U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "evidence_page_size must be between 1 and 1000000");
  }
  for (const auto *tag :
       {&config.molecule_id_tag, &config.umi_tag, &config.duplex_tag}) {
    if (!valid_sam_tag_name(*tag)) {
      throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                          "protocol SAM tags must match [A-Za-z][A-Za-z0-9]");
    }
  }
  if (config.result_schema != ResultSchema::v0_6 &&
      (!config.molecule_id_tag.empty() || !config.umi_tag.empty() ||
       !config.duplex_tag.empty())) {
    throw AnalysisError(
        AnalysisErrorCode::invalid_configuration,
        "protocol molecule tags require the schema 0.6 evidence graph");
  }
}

[[nodiscard]] std::string load_reference(const std::string &reference_path) {
  const std::string effective_path =
      reference_path.empty() ? bundled_rcrs_path() : reference_path;
  std::ifstream in(effective_path);
  if (!in) {
    throw AnalysisError(AnalysisErrorCode::reference_open_failed,
                        "could not open reference FASTA: " + effective_path);
  }

  auto reference = parse_reference_fasta(in, effective_path);
  if (reference.size() != static_cast<std::size_t>(kDefaultReferenceLength) &&
      reference_path.empty()) {
    throw AnalysisError(AnalysisErrorCode::reference_invalid,
                        "bundled rCRS reference length is not 16569 bp: " +
                            std::to_string(reference.size()));
  }
  return reference;
}

[[nodiscard]] std::vector<GeneAnnotation> default_genes() {
  return {
      {"MT-RNR1", 648, 1601, "+", "rRNA"},
      {"MT-RNR2", 1671, 3229, "+", "rRNA"},
      {"MT-ND1", 3307, 4262, "+", "protein_coding"},
      {"MT-ND2", 4470, 5511, "+", "protein_coding"},
      {"MT-CO1", 5904, 7445, "+", "protein_coding"},
      {"MT-CO2", 7586, 8269, "+", "protein_coding"},
      {"MT-ATP8", 8366, 8572, "+", "protein_coding"},
      {"MT-ATP6", 8527, 9207, "+", "protein_coding"},
      {"MT-CO3", 9207, 9990, "+", "protein_coding"},
      {"MT-ND3", 10059, 10404, "+", "protein_coding"},
      {"MT-ND4L", 10470, 10766, "+", "protein_coding"},
      {"MT-ND4", 10760, 12137, "+", "protein_coding"},
      {"MT-ND5", 12337, 14148, "+", "protein_coding"},
      {"MT-ND6", 14149, 14673, "-", "protein_coding"},
      {"MT-CYB", 14747, 15887, "+", "protein_coding"},
  };
}

[[nodiscard]] std::optional<std::size_t> parse_size(std::string_view value) {
  std::size_t parsed = 0;
  const auto *begin = value.data();
  const auto *end = value.data() + value.size();
  const auto result = std::from_chars(begin, end, parsed);
  if (result.ec != std::errc{} || result.ptr != end) {
    return std::nullopt;
  }
  return parsed;
}

[[nodiscard]] std::optional<PhyloMutation>
parse_phylo_mutation(std::string_view raw_token) {
  std::string token = trim_copy(raw_token);
  if (token.empty()) {
    return std::nullopt;
  }

  bool backmutation = false;
  if (token.back() == '!') {
    backmutation = true;
    token.pop_back();
  }
  if (token.empty() || token.find('!') != std::string::npos) {
    return std::nullopt;
  }

  std::size_t digit_count = 0;
  while (digit_count < token.size() &&
         std::isdigit(static_cast<unsigned char>(token[digit_count])) != 0) {
    ++digit_count;
  }
  if (digit_count == 0 || digit_count >= token.size()) {
    return std::nullopt;
  }
  const auto position =
      parse_size(std::string_view(token).substr(0, digit_count));
  if (!position || *position == 0 ||
      *position > static_cast<std::size_t>(kDefaultReferenceLength)) {
    return std::nullopt;
  }

  const auto canonical_base = [](char value) -> std::optional<char> {
    const char upper =
        static_cast<char>(std::toupper(static_cast<unsigned char>(value)));
    if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T') {
      return upper;
    }
    return std::nullopt;
  };

  PhyloMutation mutation;
  mutation.locus.position = *position;
  mutation.backmutation = backmutation;
  const auto suffix = std::string_view(token).substr(digit_count);
  if (suffix.size() == 1U && (suffix.front() == 'd' || suffix.front() == 'D')) {
    mutation.kind = PhyloMutationKind::deletion;
    mutation.locus.kind = PhyloMutationKind::substitution;
    mutation.encoded = std::to_string(*position) + "d";
    return mutation;
  }

  if (suffix.front() != '.') {
    if (suffix.size() != 1U) {
      return std::nullopt;
    }
    const auto alternate = canonical_base(suffix.front());
    if (!alternate) {
      return std::nullopt;
    }
    mutation.kind = PhyloMutationKind::substitution;
    mutation.locus.kind = PhyloMutationKind::substitution;
    mutation.alternate = *alternate;
    mutation.encoded = std::to_string(*position) + std::string(1, *alternate);
    return mutation;
  }

  std::size_t index_end = 1U;
  if (index_end < suffix.size() && suffix[index_end] == 'X') {
    ++index_end;
    mutation.locus.insertion_index = "X";
  } else {
    const std::size_t index_start = index_end;
    while (index_end < suffix.size() &&
           std::isdigit(static_cast<unsigned char>(suffix[index_end])) != 0) {
      ++index_end;
    }
    if (index_start == index_end) {
      return std::nullopt;
    }
    const auto insertion_index =
        parse_size(suffix.substr(index_start, index_end - index_start));
    if (!insertion_index || *insertion_index == 0) {
      return std::nullopt;
    }
    mutation.locus.insertion_index = std::to_string(*insertion_index);
  }
  if (index_end >= suffix.size()) {
    return std::nullopt;
  }

  mutation.inserted_bases.reserve(suffix.size() - index_end);
  for (const char value : suffix.substr(index_end)) {
    const auto base = canonical_base(value);
    if (!base) {
      return std::nullopt;
    }
    mutation.inserted_bases.push_back(*base);
  }
  mutation.kind = PhyloMutationKind::insertion;
  mutation.locus.kind = PhyloMutationKind::insertion;
  mutation.encoded = std::to_string(*position) + "." +
                     mutation.locus.insertion_index + mutation.inserted_bases;
  return mutation;
}

[[nodiscard]] std::vector<CigarOperation>
parse_cigar(std::string_view cigar, std::string_view read_id) {
  std::vector<CigarOperation> operations;
  if (cigar.empty() || cigar == "*") {
    return operations;
  }

  std::size_t number_start = 0;
  for (std::size_t i = 0; i < cigar.size(); ++i) {
    if (std::isdigit(static_cast<unsigned char>(cigar[i])) != 0) {
      continue;
    }

    const auto length =
        parse_size(cigar.substr(number_start, i - number_start));
    constexpr std::string_view valid_operations = "MIDNSHP=X";
    if (!length || *length == 0 ||
        valid_operations.find(cigar[i]) == std::string_view::npos) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid CIGAR for read '" + std::string(read_id) +
                              "': " + std::string(cigar));
    }
    operations.push_back({*length, cigar[i]});
    number_start = i + 1;
  }

  if (number_start != cigar.size() || operations.empty()) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "invalid CIGAR for read '" + std::string(read_id) +
                            "': " + std::string(cigar));
  }
  return operations;
}

[[nodiscard]] std::vector<SupplementaryAlignment>
parse_supplementary_alignments(std::string_view encoded_tag,
                               std::string_view read_id) {
  if (encoded_tag.empty()) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "empty SA tag for read '" + std::string(read_id) + "'");
  }

  std::vector<SupplementaryAlignment> alignments;
  std::size_t entry_start = 0;
  std::size_t entry_index = 1;
  while (entry_start < encoded_tag.size()) {
    const auto entry_end = encoded_tag.find(';', entry_start);
    if (entry_end == std::string_view::npos) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "unterminated SA entry " +
                              std::to_string(entry_index) + " for read '" +
                              std::string(read_id) + "'");
    }
    const auto entry = encoded_tag.substr(entry_start, entry_end - entry_start);
    if (entry.empty()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "empty SA entry " + std::to_string(entry_index) +
                              " for read '" + std::string(read_id) + "'");
    }

    std::array<std::string_view, 6> fields{};
    std::size_t field_start = 0;
    std::size_t field_count = 0;
    bool has_extra_field = false;
    while (field_count < fields.size()) {
      const auto comma = entry.find(',', field_start);
      fields[field_count++] =
          entry.substr(field_start, comma == std::string_view::npos
                                        ? entry.size() - field_start
                                        : comma - field_start);
      if (comma == std::string_view::npos) {
        field_start = entry.size();
        break;
      }
      field_start = comma + 1U;
      if (field_count == fields.size()) {
        has_extra_field = true;
      }
    }
    if (field_count != fields.size() || has_extra_field ||
        field_start != entry.size()) {
      throw AnalysisError(
          AnalysisErrorCode::input_parse_failed,
          "SA entry " + std::to_string(entry_index) +
              " must contain exactly 6 comma-separated fields for read '" +
              std::string(read_id) + "'");
    }

    const auto position = parse_size(fields[1]);
    const auto mapq = parse_size(fields[4]);
    const auto edit_distance = parse_size(fields[5]);
    if (fields[0].empty() || fields[0] == "*" || !position || *position == 0 ||
        (fields[2] != "+" && fields[2] != "-") || !mapq ||
        *mapq > std::numeric_limits<std::uint8_t>::max() || !edit_distance) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid SA entry " + std::to_string(entry_index) +
                              " for read '" + std::string(read_id) + "'");
    }

    SupplementaryAlignment alignment;
    alignment.reference_name = std::string(fields[0]);
    alignment.reference_start = *position;
    alignment.strand = fields[2].front();
    alignment.cigar = std::string(fields[3]);
    alignment.cigar_operations =
        parse_cigar(fields[3], std::string(read_id) + " SA");
    if (alignment.cigar_operations.empty()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "SA entry " + std::to_string(entry_index) +
                              " has no alignment CIGAR for read '" +
                              std::string(read_id) + "'");
    }
    alignment.mapping_quality = static_cast<std::uint8_t>(*mapq);
    alignment.edit_distance = *edit_distance;
    alignments.push_back(std::move(alignment));

    entry_start = entry_end + 1U;
    ++entry_index;
  }
  return alignments;
}

[[nodiscard]] std::map<std::string, std::string>
parse_sam_aux_tags(const std::vector<std::string> &fields,
                   std::size_t line_number) {
  std::map<std::string, std::string> tags;
  for (std::size_t i = 11; i < fields.size(); ++i) {
    if (fields[i].size() < 5 || fields[i][2] != ':' || fields[i][4] != ':') {
      if (fields[i].rfind("SA", 0) == 0) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "malformed SA auxiliary field at SAM line " +
                                std::to_string(line_number));
      }
      continue;
    }
    if (fields[i].compare(0, 2, "SA") == 0 && fields[i][3] != 'Z') {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "SA auxiliary field must use SAM type Z at line " +
                              std::to_string(line_number));
    }
    tags.emplace(fields[i].substr(0, 2), fields[i].substr(5));
  }
  return tags;
}

[[nodiscard]] std::string snp_key(std::size_t position, char reference,
                                  char alternate) {
  std::string key = std::to_string(position);
  key.push_back(':');
  key.push_back(
      static_cast<char>(std::toupper(static_cast<unsigned char>(reference))));
  key.push_back(':');
  key.push_back(
      static_cast<char>(std::toupper(static_cast<unsigned char>(alternate))));
  return key;
}

[[nodiscard]] std::string lowercase_token(std::string_view value) {
  std::string normalized;
  normalized.reserve(value.size());
  for (const char c : value) {
    const auto byte = static_cast<unsigned char>(c);
    if (std::isspace(byte) != 0 || c == '-') {
      if (normalized.empty() || normalized.back() != '_') {
        normalized.push_back('_');
      }
    } else {
      normalized.push_back(static_cast<char>(std::tolower(byte)));
    }
  }
  while (!normalized.empty() && normalized.back() == '_') {
    normalized.pop_back();
  }
  return normalized;
}

[[nodiscard]] std::string
normalize_clinical_significance(std::string_view significance) {
  const auto value = lowercase_token(trim_copy(significance));
  if (value.empty() || value == "not_provided" || value == "not_specified") {
    return "not_provided";
  }
  if (value.find("conflict") != std::string::npos ||
      (value.find("pathogenic") != std::string::npos &&
       value.find("benign") != std::string::npos)) {
    return "conflicting";
  }
  if (value.find("pathogenic/likely_pathogenic") != std::string::npos ||
      value.find("pathogenic,_likely_pathogenic") != std::string::npos) {
    return "pathogenic_or_likely_pathogenic";
  }
  if (value.find("likely_pathogenic") != std::string::npos) {
    return "likely_pathogenic";
  }
  if (value.find("pathogenic") != std::string::npos) {
    return "pathogenic";
  }
  if (value.find("benign/likely_benign") != std::string::npos ||
      value.find("benign,_likely_benign") != std::string::npos) {
    return "benign_or_likely_benign";
  }
  if (value.find("likely_benign") != std::string::npos) {
    return "likely_benign";
  }
  if (value.find("benign") != std::string::npos) {
    return "benign";
  }
  if (value.find("uncertain") != std::string::npos || value == "vus") {
    return "uncertain_significance";
  }
  for (const auto category : {"risk_factor", "association", "drug_response",
                              "protective", "affects"}) {
    if (value.find(category) != std::string::npos) {
      return category;
    }
  }
  return "other";
}

[[nodiscard]] std::string
clinical_significance_group(std::string_view normalized) {
  if (normalized == "pathogenic" || normalized == "likely_pathogenic" ||
      normalized == "pathogenic_or_likely_pathogenic") {
    return "pathogenic";
  }
  if (normalized == "benign" || normalized == "likely_benign" ||
      normalized == "benign_or_likely_benign") {
    return "benign";
  }
  if (normalized == "uncertain_significance") {
    return "uncertain";
  }
  if (normalized == "conflicting") {
    return "conflicting";
  }
  if (normalized == "not_provided" || normalized == "other") {
    return "unclassified";
  }
  return "other_assertion";
}

[[nodiscard]] std::string join_values(const std::vector<std::string> &values,
                                      std::string_view delimiter) {
  std::string joined;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0U) {
      joined += delimiter;
    }
    joined += values[i];
  }
  return joined;
}

void finalize_clinical_annotation(SnpCall &call) {
  auto &assertions = call.clinical_assertions;
  std::sort(assertions.begin(), assertions.end());
  assertions.erase(std::unique(assertions.begin(), assertions.end()),
                   assertions.end());

  std::vector<std::string> diseases;
  std::vector<std::string> normalized_significances;
  std::vector<std::string> significance_groups;
  call.references.clear();
  call.sources.clear();
  call.clinvar_allele_id.clear();
  call.mitomap_url.clear();
  for (const auto &assertion : assertions) {
    if (!assertion.disease.empty()) {
      diseases.push_back(assertion.disease);
    }
    if (assertion.normalized_significance != "not_provided") {
      normalized_significances.push_back(assertion.normalized_significance);
    }
    const auto group =
        clinical_significance_group(assertion.normalized_significance);
    if (group != "unclassified") {
      significance_groups.push_back(group);
    }
    call.sources.push_back(assertion.source);
    call.references.insert(call.references.end(), assertion.references.begin(),
                           assertion.references.end());
    if (assertion.source == "ClinVar" && call.clinvar_allele_id.empty()) {
      call.clinvar_allele_id = assertion.allele_id;
    }
    if (assertion.source == "MITOMAP" && call.mitomap_url.empty()) {
      call.mitomap_url = assertion.source_url;
    }
  }
  for (auto *values : {&diseases, &normalized_significances,
                       &significance_groups, &call.references, &call.sources}) {
    std::sort(values->begin(), values->end());
    values->erase(std::unique(values->begin(), values->end()), values->end());
  }

  const bool explicit_conflict =
      std::find(significance_groups.begin(), significance_groups.end(),
                "conflicting") != significance_groups.end();
  const std::size_t classified_group_count =
      static_cast<std::size_t>(std::count_if(
          significance_groups.begin(), significance_groups.end(),
          [](const std::string &group) { return group != "conflicting"; }));
  const bool incompatible_groups = classified_group_count > 1U;
  if (explicit_conflict || incompatible_groups) {
    call.clinical_conflict_status = "conflicting";
    call.clinical_consensus_significance = "conflicting";
  } else if (assertions.size() == 1U) {
    call.clinical_conflict_status = "single_assertion";
  } else {
    call.clinical_conflict_status = "consistent";
  }

  if (call.clinical_consensus_significance.empty()) {
    if (normalized_significances.empty()) {
      call.clinical_consensus_significance = "not_provided";
    } else if (normalized_significances.size() == 1U) {
      call.clinical_consensus_significance = normalized_significances.front();
    } else if (significance_groups.size() == 1U &&
               significance_groups.front() == "pathogenic") {
      call.clinical_consensus_significance = "pathogenic_or_likely_pathogenic";
    } else if (significance_groups.size() == 1U &&
               significance_groups.front() == "benign") {
      call.clinical_consensus_significance = "benign_or_likely_benign";
    } else {
      call.clinical_consensus_significance = normalized_significances.front();
    }
  }
  call.phenotype = join_values(diseases, "; ");
  call.pathogenicity = call.clinical_consensus_significance;
}

[[nodiscard]] std::map<std::string, SnpCall> load_clinical_annotations() {
  std::map<std::string, SnpCall> annotations;
  std::ifstream input(clinical_annotations_path());
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open clinical annotations: " +
                            clinical_annotations_path());
  }

  constexpr std::string_view expected_header =
      "position\tref\talt\tgene\tconsequence\tprotein\tresidue\tstructure_id\t"
      "structure_chain\tstructure_residue\tstructure_"
      "complex\tsource\tassertion_id\t"
      "allele_id\tdisease\tclinical_significance\treview_status\tassertion_"
      "date\t"
      "source_url\treferences\tresource_version\tretrieved_at";
  std::string line;
  std::size_t line_number = 0U;
  bool header_seen = false;
  while (std::getline(input, line)) {
    ++line_number;
    strip_trailing_carriage_return(line);
    if (line.empty()) {
      continue;
    }
    if (!header_seen) {
      header_seen = true;
      if (line != expected_header) {
        throw AnalysisError(
            AnalysisErrorCode::resource_invalid,
            "clinical annotation header does not match schema 1.0");
      }
      continue;
    }

    auto fields = split_tab(line);
    if (fields.size() != 22U) {
      throw AnalysisError(
          AnalysisErrorCode::resource_invalid,
          "clinical annotation row must contain 22 fields at line " +
              std::to_string(line_number));
    }
    for (auto &field : fields) {
      field = trim_copy(field);
    }
    const auto position = parse_size(fields[0]);
    const auto valid_allele = [](const std::string &allele) {
      if (allele.size() != 1U) {
        return false;
      }
      const char base = static_cast<char>(
          std::toupper(static_cast<unsigned char>(allele.front())));
      return base == 'A' || base == 'C' || base == 'G' || base == 'T';
    };
    if (!position || *position == 0U ||
        *position > static_cast<std::size_t>(kDefaultReferenceLength) ||
        !valid_allele(fields[1]) || !valid_allele(fields[2]) ||
        std::toupper(static_cast<unsigned char>(fields[1].front())) ==
            std::toupper(static_cast<unsigned char>(fields[2].front())) ||
        fields[11].empty()) {
      throw AnalysisError(AnalysisErrorCode::resource_invalid,
                          "invalid clinical variant/assertion key at line " +
                              std::to_string(line_number));
    }
    if (!fields[18].empty() && !fields[18].starts_with("https://") &&
        !fields[18].starts_with("http://")) {
      throw AnalysisError(AnalysisErrorCode::resource_invalid,
                          "clinical source_url must be HTTP(S) at line " +
                              std::to_string(line_number));
    }

    SnpCall call;
    call.position = *position;
    call.reference = static_cast<char>(
        std::toupper(static_cast<unsigned char>(fields[1].front())));
    call.alternate = static_cast<char>(
        std::toupper(static_cast<unsigned char>(fields[2].front())));
    call.gene = fields[3];
    call.consequence = fields[4];
    call.protein = fields[5];
    call.residue = fields[6];
    call.structure_id = fields[7];
    call.structure_chain = fields[8];
    if (!fields[9].empty()) {
      const auto residue = parse_size(fields[9]);
      if (!residue || *residue == 0U) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "invalid clinical structure residue at line " +
                                std::to_string(line_number));
      }
      call.structure_residue = *residue;
    }
    call.structure_complex = fields[10];
    ClinicalAssertion assertion;
    assertion.source = fields[11];
    assertion.assertion_id = fields[12];
    assertion.allele_id = fields[13];
    assertion.disease = fields[14];
    assertion.clinical_significance = fields[15];
    assertion.normalized_significance =
        normalize_clinical_significance(fields[15]);
    assertion.review_status = fields[16];
    assertion.assertion_date = fields[17];
    assertion.source_url = fields[18];
    assertion.references = split_semicolon(fields[19]);
    assertion.resource_version = fields[20];
    assertion.retrieved_at = fields[21];
    call.clinical_assertions.push_back(std::move(assertion));

    const auto key = snp_key(call.position, call.reference, call.alternate);
    auto [it, inserted] = annotations.try_emplace(key, call);
    if (!inserted) {
      auto &existing = it->second;
      const auto merge_scalar = [&](std::string &target,
                                    const std::string &incoming,
                                    std::string_view name) {
        if (target.empty()) {
          target = incoming;
        } else if (!incoming.empty() && target != incoming) {
          throw AnalysisError(AnalysisErrorCode::resource_invalid,
                              "conflicting clinical " + std::string(name) +
                                  " at line " + std::to_string(line_number));
        }
      };
      merge_scalar(existing.gene, call.gene, "gene");
      merge_scalar(existing.consequence, call.consequence, "consequence");
      merge_scalar(existing.protein, call.protein, "protein");
      merge_scalar(existing.residue, call.residue, "residue");
      merge_scalar(existing.structure_id, call.structure_id, "structure_id");
      merge_scalar(existing.structure_chain, call.structure_chain,
                   "structure_chain");
      merge_scalar(existing.structure_complex, call.structure_complex,
                   "structure_complex");
      if (existing.structure_residue == 0U) {
        existing.structure_residue = call.structure_residue;
      } else if (call.structure_residue != 0U &&
                 existing.structure_residue != call.structure_residue) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "conflicting clinical structure_residue at line " +
                                std::to_string(line_number));
      }
      existing.clinical_assertions.insert(existing.clinical_assertions.end(),
                                          call.clinical_assertions.begin(),
                                          call.clinical_assertions.end());
    }
  }

  if (!header_seen || annotations.empty()) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "clinical annotation resource contains no assertions");
  }
  for (auto &[_, call] : annotations) {
    finalize_clinical_annotation(call);
  }
  return annotations;
}

[[nodiscard]] InputRecords parse_fastq(std::istream &input,
                                       const AnalysisConfig &config) {
  InputRecords records;
  std::string header;
  std::string sequence;
  std::string plus;
  std::string quality;

  while (std::getline(input, header)) {
    throw_if_cancelled(config);
    strip_trailing_carriage_return(header);
    if (header.empty()) {
      continue;
    }
    if (header.front() != '@') {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "FASTQ parser expected a header starting with '@'");
    }
    if (!std::getline(input, sequence) || !std::getline(input, plus) ||
        !std::getline(input, quality)) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "FASTQ record is truncated");
    }
    strip_trailing_carriage_return(sequence);
    strip_trailing_carriage_return(plus);
    strip_trailing_carriage_return(quality);
    if (plus.empty() || plus.front() != '+') {
      throw AnalysisError(
          AnalysisErrorCode::input_parse_failed,
          "FASTQ parser expected a separator starting with '+' for read '" +
              header.substr(1) + "'");
    }
    if (sequence.empty()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "FASTQ sequence is empty for read '" +
                              header.substr(1) + "'");
    }
    if (quality.size() != sequence.size()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "FASTQ sequence/quality length mismatch for read '" +
                              header.substr(1) + "'");
    }
    ReadRecord record;
    record.id = header.substr(1);
    record.sequence = sequence;
    record.qualities = quality;
    records.reads.push_back(std::move(record));
  }

  return records;
}

[[nodiscard]] InputRecords parse_sam(std::istream &input,
                                     const AnalysisConfig &config) {
  InputRecords records;
  records.alignment_input = true;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    throw_if_cancelled(config);
    ++line_number;
    strip_trailing_carriage_return(line);
    if (line.rfind("@SQ\t", 0) == 0) {
      const auto header_fields = split_tab(line);
      for (const auto &field : header_fields) {
        if (field.rfind("SN:", 0) == 0 &&
            looks_like_nuclear_contig(field.substr(3))) {
          records.has_nuclear_contigs = true;
        }
      }
      continue;
    }
    if (line.empty() || line.front() == '@') {
      continue;
    }

    const auto fields = split_tab(line);
    if (fields.size() < 11) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "SAM record has fewer than 11 fields at line " +
                              std::to_string(line_number));
    }

    const auto flags = parse_size(fields[1]);
    const auto reference_start = parse_size(fields[3]);
    const auto mapping_quality = parse_size(fields[4]);
    if (!flags || *flags > std::numeric_limits<std::uint16_t>::max()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid SAM FLAG at line " +
                              std::to_string(line_number));
    }
    if (!reference_start) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid SAM POS at line " +
                              std::to_string(line_number));
    }
    if (!mapping_quality ||
        *mapping_quality > std::numeric_limits<std::uint8_t>::max()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid SAM MAPQ at line " +
                              std::to_string(line_number));
    }
    if (fields[9] != "*" && fields[10] != "*" &&
        fields[9].size() != fields[10].size()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "SAM sequence/quality length mismatch at line " +
                              std::to_string(line_number));
    }

    ReadRecord record;
    record.id = fields[0];
    record.flags = static_cast<std::uint16_t>(*flags);
    record.reference_name = fields[2];
    record.reference_start = *reference_start;
    record.mapping_quality = static_cast<std::uint8_t>(*mapping_quality);
    record.cigar = fields[5];
    record.cigar_operations = parse_cigar(record.cigar, record.id);
    record.sequence = fields[9] == "*" ? std::string{} : fields[9];
    record.qualities = fields[10] == "*" ? std::string{} : fields[10];
    if (!record.sequence.empty() && !record.cigar_operations.empty()) {
      std::size_t query_length = 0;
      for (const auto &operation : record.cigar_operations) {
        if (operation.code == 'M' || operation.code == 'I' ||
            operation.code == 'S' || operation.code == '=' ||
            operation.code == 'X') {
          if (operation.length >
              std::numeric_limits<std::size_t>::max() - query_length) {
            throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                                "SAM CIGAR query length overflow at line " +
                                    std::to_string(line_number));
          }
          query_length += operation.length;
        }
      }
      if (query_length != record.sequence.size()) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "SAM CIGAR/query length mismatch at line " +
                                std::to_string(line_number));
      }
    }
    record.aux_tags = parse_sam_aux_tags(fields, line_number);
    records.reads.push_back(std::move(record));
  }
  return records;
}

#ifdef MITO_HAS_HTSLIB
class MitoFileReader {
public:
  explicit MitoFileReader(std::string path) : path_(std::move(path)) {}

  [[nodiscard]] InputRecords read_all(const AnalysisConfig &config) const {
    htsFile *raw_file = sam_open(path_.c_str(), "r");
    if (raw_file == nullptr) {
      throw AnalysisError(AnalysisErrorCode::input_open_failed,
                          "htslib could not open alignment file: " + path_);
    }
    HtsFileHandle file(raw_file);

    bam_hdr_t *raw_header = sam_hdr_read(file.get());
    if (raw_header == nullptr) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "htslib could not read alignment header: " + path_);
    }
    BamHeaderHandle header(raw_header);

    bam1_t *raw_record = bam_init1();
    if (raw_record == nullptr) {
      throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                          "htslib could not allocate BAM record");
    }
    BamRecordHandle record(raw_record);

    InputRecords records;
    records.alignment_input = true;
    for (int i = 0; i < header.get()->n_targets; ++i) {
      if (looks_like_nuclear_contig(header.get()->target_name[i])) {
        records.has_nuclear_contigs = true;
        break;
      }
    }
    int read_status = 0;
    while ((read_status = sam_read1(file.get(), header.get(), record.get())) >=
           0) {
      throw_if_cancelled(config);
      records.reads.push_back(convert_record(*header.get(), *record.get()));
    }
    if (read_status < -1) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "htslib failed while reading alignment records: " +
                              path_);
    }
    return records;
  }

private:
  struct HtsFileHandle {
    explicit HtsFileHandle(htsFile *value) : value_(value) {}
    ~HtsFileHandle() { hts_close(value_); }
    HtsFileHandle(const HtsFileHandle &) = delete;
    HtsFileHandle &operator=(const HtsFileHandle &) = delete;
    [[nodiscard]] htsFile *get() const { return value_; }
    htsFile *value_;
  };

  struct BamHeaderHandle {
    explicit BamHeaderHandle(bam_hdr_t *value) : value_(value) {}
    ~BamHeaderHandle() { bam_hdr_destroy(value_); }
    BamHeaderHandle(const BamHeaderHandle &) = delete;
    BamHeaderHandle &operator=(const BamHeaderHandle &) = delete;
    [[nodiscard]] bam_hdr_t *get() const { return value_; }
    bam_hdr_t *value_;
  };

  struct BamRecordHandle {
    explicit BamRecordHandle(bam1_t *value) : value_(value) {}
    ~BamRecordHandle() { bam_destroy1(value_); }
    BamRecordHandle(const BamRecordHandle &) = delete;
    BamRecordHandle &operator=(const BamRecordHandle &) = delete;
    [[nodiscard]] bam1_t *get() const { return value_; }
    bam1_t *value_;
  };

  [[nodiscard]] static ReadRecord convert_record(const bam_hdr_t &header,
                                                 const bam1_t &record) {
    ReadRecord out;
    out.id = bam_get_qname(&record);
    out.flags = record.core.flag;
    out.mapping_quality = record.core.qual;
    out.reference_start =
        record.core.pos < 0 ? 1 : static_cast<std::size_t>(record.core.pos) + 1;
    if (record.core.tid >= 0 && record.core.tid < header.n_targets) {
      out.reference_name = header.target_name[record.core.tid];
    }
    out.cigar = cigar_string(record);
    out.cigar_operations = parse_cigar(out.cigar, out.id);
    out.sequence = sequence_string(record);
    out.qualities = quality_string(record);
    out.aux_tags = aux_tags(record);
    return out;
  }

  [[nodiscard]] static std::string cigar_string(const bam1_t &record) {
    std::string cigar;
    const auto *raw_cigar = bam_get_cigar(&record);
    for (std::uint32_t i = 0; i < record.core.n_cigar; ++i) {
      cigar += std::to_string(bam_cigar_oplen(raw_cigar[i]));
      cigar.push_back(bam_cigar_opchr(raw_cigar[i]));
    }
    return cigar;
  }

  [[nodiscard]] static std::string sequence_string(const bam1_t &record) {
    std::string sequence;
    sequence.reserve(static_cast<std::size_t>(record.core.l_qseq));
    const auto *seq = bam_get_seq(&record);
    for (int i = 0; i < record.core.l_qseq; ++i) {
      sequence.push_back(seq_nt16_str[bam_seqi(seq, i)]);
    }
    return sequence;
  }

  [[nodiscard]] static std::string quality_string(const bam1_t &record) {
    std::string qualities;
    qualities.reserve(static_cast<std::size_t>(record.core.l_qseq));
    const auto *qual = bam_get_qual(&record);
    for (int i = 0; i < record.core.l_qseq; ++i) {
      qualities.push_back(qual[i] == 0xFFU ? '!'
                                           : static_cast<char>(qual[i] + 33U));
    }
    return qualities;
  }

  [[nodiscard]] static std::map<std::string, std::string>
  aux_tags(const bam1_t &record) {
    std::map<std::string, std::string> tags;
    errno = 0;
    const std::uint8_t *value = bam_aux_first(&record);
    if (value == nullptr) {
      if (errno == ENOENT) {
        return tags;
      }
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid auxiliary fields for read '" +
                              std::string(bam_get_qname(&record)) + "'");
    }
    while (value != nullptr) {
      const auto *raw_tag = bam_aux_tag(value);
      const std::string tag(raw_tag, raw_tag + 2);
      const char type = bam_aux_type(value);
      if (tag == "SA" && type != 'Z') {
        throw AnalysisError(
            AnalysisErrorCode::input_parse_failed,
            "SA auxiliary field must use SAM type Z for read '" +
                std::string(bam_get_qname(&record)) + "'");
      }
      std::string encoded;
      if (type == 'Z' || type == 'H') {
        if (const char *text = bam_aux2Z(value)) {
          encoded = text;
        } else {
          throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                              "invalid string auxiliary field " + tag +
                                  " for read '" +
                                  std::string(bam_get_qname(&record)) + "'");
        }
      } else if (type == 'A') {
        encoded.assign(1U, bam_aux2A(value));
      } else if (type == 'c' || type == 'C' || type == 's' || type == 'S' ||
                 type == 'i' || type == 'I') {
        encoded = std::to_string(bam_aux2i(value));
      } else if (type == 'f' || type == 'd') {
        std::ostringstream out;
        out << std::setprecision(std::numeric_limits<double>::max_digits10)
            << bam_aux2f(value);
        encoded = std::move(out).str();
      } else if (type == 'B') {
        const char element_type = static_cast<char>(value[1]);
        const auto length = bam_auxB_len(value);
        std::ostringstream out;
        out << element_type;
        for (std::uint32_t index = 0U; index < length; ++index) {
          out << ',';
          if (element_type == 'f') {
            out << std::setprecision(std::numeric_limits<double>::max_digits10)
                << bam_auxB2f(value, index);
          } else {
            out << bam_auxB2i(value, index);
          }
        }
        encoded = std::move(out).str();
      } else {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "unsupported auxiliary field type for " + tag +
                                " on read '" +
                                std::string(bam_get_qname(&record)) + "'");
      }
      if (!tags.emplace(tag, std::move(encoded)).second) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "duplicate auxiliary field " + tag + " for read '" +
                                std::string(bam_get_qname(&record)) + "'");
      }
      errno = 0;
      value = bam_aux_next(&record, value);
      if (value == nullptr && errno != ENOENT) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "invalid auxiliary fields for read '" +
                                std::string(bam_get_qname(&record)) + "'");
      }
    }
    return tags;
  }

  std::string path_;
};
#endif

[[nodiscard]] bool
has_extension(const std::filesystem::path &path,
              std::initializer_list<std::string_view> extensions) {
  auto ext = path.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char value) {
    return static_cast<char>(std::tolower(value));
  });
  return std::any_of(extensions.begin(), extensions.end(),
                     [&](std::string_view value) { return ext == value; });
}

[[nodiscard]] InputRecords read_input(const std::string &input_path,
                                      const AnalysisConfig &config) {
  const std::filesystem::path path(input_path);
  if (has_extension(path, {".fa", ".fas", ".fasta", ".fna"})) {
    throw AnalysisError(
        AnalysisErrorCode::input_format_unsupported,
        "FASTA input is not supported; use FASTQ, SAM, BAM, or CRAM");
  }

#ifdef MITO_HAS_HTSLIB
  if (has_extension(path, {".sam", ".bam", ".cram"})) {
    return MitoFileReader(input_path).read_all(config);
  }
#endif

  std::ifstream input(input_path, std::ios::binary);
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::input_open_failed,
                        "could not open input file: " + input_path);
  }

  if (has_extension(path, {".sam"})) {
    return parse_sam(input, config);
  }
  if (has_extension(path, {".bam", ".cram"})) {
    throw AnalysisError(AnalysisErrorCode::dependency_unavailable,
                        "BAM/CRAM support requires htslib; install htslib "
                        "development headers and rebuild");
  }
  return parse_fastq(input, config);
}

[[nodiscard]] std::string alignment_as_sa(const ReadRecord &read) {
  const char strand = (read.flags & 0x10U) != 0U ? '-' : '+';
  const auto nm = read.aux_tags.find("NM");
  return read.reference_name + "," + std::to_string(read.reference_start) +
         "," + strand + "," + (read.cigar.empty() ? "*" : read.cigar) + "," +
         std::to_string(static_cast<unsigned int>(read.mapping_quality)) + "," +
         (nm == read.aux_tags.end() ? "0" : nm->second) + ";";
}

[[nodiscard]] std::string alignment_fragment_role(const std::uint16_t flags,
                                                  const bool alignment_input) {
  if (!alignment_input) {
    return "unaligned_read";
  }
  const bool secondary = (flags & 0x100U) != 0U;
  const bool supplementary = (flags & 0x800U) != 0U;
  if (secondary && supplementary) {
    return "secondary_and_supplementary";
  }
  if (secondary) {
    return "secondary";
  }
  if (supplementary) {
    return "supplementary";
  }
  return "primary_candidate";
}

void validate_protocol_tag_value(const std::string &tag,
                                 const std::string &value,
                                 const std::string &read_id) {
  if (value.empty() || value.size() > 256U ||
      std::any_of(value.begin(), value.end(), [](const char value_char) {
        const auto character = static_cast<unsigned char>(value_char);
        return std::iscntrl(character) != 0;
      })) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "invalid protocol tag " + tag + " on read '" + read_id +
                            "'");
  }
}

struct ProtocolTagValues {
  std::set<std::string> values;
  bool missing = false;
};

[[nodiscard]] ProtocolTagValues
protocol_tag_values(const std::vector<ReadRecord> &reads,
                    const std::vector<std::size_t> &indices,
                    const std::string &tag) {
  ProtocolTagValues result;
  if (tag.empty()) {
    return result;
  }
  for (const auto index : indices) {
    const auto found = reads[index].aux_tags.find(tag);
    if (found == reads[index].aux_tags.end()) {
      result.missing = true;
      continue;
    }
    validate_protocol_tag_value(tag, found->second, reads[index].id);
    result.values.insert(found->second);
  }
  return result;
}

[[nodiscard]] std::size_t reference_consuming_length(const ReadRecord &read) {
  std::size_t length = 0U;
  for (const auto &operation : read.cigar_operations) {
    if (operation.code != 'M' && operation.code != '=' &&
        operation.code != 'X' && operation.code != 'D' &&
        operation.code != 'N') {
      continue;
    }
    if (operation.length > std::numeric_limits<std::size_t>::max() - length) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "CIGAR reference span overflows for read '" +
                              read.id + "'");
    }
    length += operation.length;
  }
  return length;
}

void sort_unique_strings(std::vector<std::string> &values) {
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
}

[[nodiscard]] MoleculeAssemblyResult
assemble_molecules(InputRecords &input, const AnalysisConfig &config,
                   const std::size_t reference_length) {
  MoleculeAssemblyResult result;
  if (!input.alignment_input) {
    if (!config.molecule_id_tag.empty() || !config.umi_tag.empty() ||
        !config.duplex_tag.empty()) {
      throw AnalysisError(
          AnalysisErrorCode::invalid_configuration,
          "protocol SAM tags cannot be applied to unaligned FASTQ input");
    }
    result.representatives = std::move(input.reads);
    result.fragments.reserve(result.representatives.size());
    result.molecules.reserve(result.representatives.size());
    std::unordered_map<std::string, std::size_t> name_occurrences;
    name_occurrences.reserve(result.representatives.size());
    for (std::size_t index = 0; index < result.representatives.size();
         ++index) {
      const auto &read = result.representatives[index];
      const auto occurrence = ++name_occurrences[read.id];
      const std::string molecule_id =
          occurrence == 1U ? read.id
                           : read.id + "#record:" + std::to_string(index);
      result.fragments.push_back(
          {AlignmentFragmentId{index}, MoleculeIndex{index}, index, molecule_id,
           alignment_fragment_role(read.flags, false), true, read.flags,
           read.mapping_quality, read.reference_name, read.reference_start,
           read.cigar});
      MoleculeAssemblyEvidence molecule;
      molecule.index = MoleculeIndex{index};
      molecule.id = molecule_id;
      molecule.identity_policy = "fastq_record_proxy";
      molecule.assembly_status = "single_unaligned_fragment";
      molecule.fragment_ids.push_back(AlignmentFragmentId{index});
      molecule.representative_fragment_id = AlignmentFragmentId{index};
      molecule.primary_candidate_count = 1U;
      molecule.source_qnames.push_back(read.id);
      molecule.protocol_flags.emplace_back("READ_PROXY_IDENTITY");
      if (occurrence > 1U) {
        molecule.ambiguous = true;
        molecule.warnings.emplace_back("DUPLICATE_SOURCE_NAME_DISAMBIGUATED");
      }
      result.molecules.push_back(std::move(molecule));
    }
    return result;
  }

  for (auto &record : input.reads) {
    const auto sa = record.aux_tags.find("SA");
    if (sa != record.aux_tags.end()) {
      record.supplementary_alignments =
          parse_supplementary_alignments(sa->second, record.id);
    }
  }

  std::map<std::string, std::vector<std::size_t>> qname_groups;
  for (std::size_t index = 0; index < input.reads.size(); ++index) {
    qname_groups[input.reads[index].id].push_back(index);
  }

  struct AssemblyGroup {
    std::string molecule_id;
    std::vector<std::size_t> indices;
    bool missing_identity_tag = false;
    bool conflicting_identity_tag = false;
    bool partial_identity_tag = false;
  };
  std::map<std::string, AssemblyGroup> groups;
  for (const auto &[qname, indices] : qname_groups) {
    if (config.molecule_id_tag.empty()) {
      auto &group = groups["qname:" + qname];
      group.molecule_id = qname;
      group.indices.insert(group.indices.end(), indices.begin(), indices.end());
      continue;
    }
    const auto tag_values =
        protocol_tag_values(input.reads, indices, config.molecule_id_tag);
    if (tag_values.values.empty()) {
      auto &group = groups["missing:" + qname];
      group.molecule_id = "UNASSIGNED:" + config.molecule_id_tag + ":" + qname;
      group.missing_identity_tag = true;
      group.indices.insert(group.indices.end(), indices.begin(), indices.end());
      continue;
    }
    if (tag_values.values.size() != 1U) {
      auto &group = groups["conflict:" + qname];
      group.molecule_id = "CONFLICT:" + config.molecule_id_tag + ":" + qname;
      group.conflicting_identity_tag = true;
      group.indices.insert(group.indices.end(), indices.begin(), indices.end());
      continue;
    }
    const auto &value = *tag_values.values.begin();
    auto &group = groups["tag:" + value];
    group.molecule_id = config.molecule_id_tag + ":" + value;
    group.partial_identity_tag =
        group.partial_identity_tag || tag_values.missing;
    group.indices.insert(group.indices.end(), indices.begin(), indices.end());
  }

  result.representatives.reserve(groups.size());
  result.fragments.resize(input.reads.size());
  result.molecules.reserve(groups.size());
  for (const auto &[_, group] : groups) {
    const auto &molecule_id = group.molecule_id;
    const auto &indices = group.indices;
    const MoleculeIndex molecule_index{result.molecules.size()};
    std::vector<std::size_t> primary_candidates;
    primary_candidates.reserve(indices.size());
    std::map<std::string, std::size_t> qname_primary_counts;
    std::set<std::string> source_qnames;
    for (const auto index : indices) {
      source_qnames.insert(input.reads[index].id);
      if ((input.reads[index].flags & (0x100U | 0x800U)) == 0U) {
        primary_candidates.push_back(index);
        ++qname_primary_counts[input.reads[index].id];
      } else {
        qname_primary_counts.try_emplace(input.reads[index].id, 0U);
      }
    }
    std::size_t representative_index = primary_candidates.empty()
                                           ? indices.front()
                                           : primary_candidates.front();
    if (!config.molecule_id_tag.empty() && !primary_candidates.empty()) {
      representative_index = *std::max_element(
          primary_candidates.begin(), primary_candidates.end(),
          [&](const auto lhs, const auto rhs) {
            if (input.reads[lhs].mapping_quality !=
                input.reads[rhs].mapping_quality) {
              return input.reads[lhs].mapping_quality <
                     input.reads[rhs].mapping_quality;
            }
            return lhs > rhs;
          });
    }

    MoleculeAssemblyEvidence molecule;
    molecule.index = molecule_index;
    molecule.id = molecule_id;
    molecule.identity_policy = config.molecule_id_tag.empty()
                                   ? "sam_qname"
                                   : "sam_tag:" + config.molecule_id_tag;
    molecule.representative_fragment_id =
        AlignmentFragmentId{representative_index};
    molecule.primary_candidate_count = primary_candidates.size();
    molecule.source_qnames.assign(source_qnames.begin(), source_qnames.end());

    const bool malformed_qname_primary =
        std::any_of(qname_primary_counts.begin(), qname_primary_counts.end(),
                    [](const auto &entry) { return entry.second != 1U; });
    if (config.molecule_id_tag.empty()) {
      molecule.ambiguous = primary_candidates.size() != 1U;
      molecule.analysis_eligible = primary_candidates.size() == 1U;
      molecule.protocol_flags.emplace_back("READ_PROXY_IDENTITY");
      if (primary_candidates.empty()) {
        molecule.assembly_status = "fallback_without_primary";
        molecule.warnings.emplace_back("NO_PRIMARY_ALIGNMENT");
        molecule.exclusion_reasons.emplace_back("NO_PRIMARY_ALIGNMENT");
      } else if (primary_candidates.size() > 1U) {
        molecule.assembly_status = "first_of_multiple_primaries";
        molecule.warnings.emplace_back("MULTIPLE_PRIMARY_ALIGNMENTS");
        molecule.exclusion_reasons.emplace_back("MULTIPLE_PRIMARY_ALIGNMENTS");
      } else {
        molecule.assembly_status = "unique_primary";
      }
    } else {
      molecule.protocol_flags.emplace_back("EXPLICIT_MOLECULE_ID");
      molecule.protocol_metadata["molecule_id_tag"] = config.molecule_id_tag;
      if (!group.missing_identity_tag && !group.conflicting_identity_tag) {
        molecule.protocol_metadata["molecule_id_value"] =
            molecule_id.substr(config.molecule_id_tag.size() + 1U);
      }
      molecule.ambiguous = group.missing_identity_tag ||
                           group.conflicting_identity_tag ||
                           malformed_qname_primary;
      molecule.analysis_eligible = !molecule.ambiguous;
      if (group.missing_identity_tag) {
        molecule.assembly_status = "missing_explicit_molecule_id";
        molecule.warnings.emplace_back("MISSING_MOLECULE_ID_TAG");
        molecule.exclusion_reasons.emplace_back("MISSING_MOLECULE_ID_TAG");
      } else if (group.conflicting_identity_tag) {
        molecule.assembly_status = "conflicting_explicit_molecule_id";
        molecule.warnings.emplace_back("CONFLICTING_MOLECULE_ID_TAG");
        molecule.exclusion_reasons.emplace_back("CONFLICTING_MOLECULE_ID_TAG");
      } else if (malformed_qname_primary) {
        molecule.assembly_status = "invalid_tagged_fragment_group";
        molecule.warnings.emplace_back("INVALID_PRIMARY_COUNT_PER_QNAME");
        molecule.exclusion_reasons.emplace_back(
            "INVALID_PRIMARY_COUNT_PER_QNAME");
      } else if (source_qnames.size() > 1U) {
        molecule.assembly_status = "explicit_tag_multi_qname_group";
      } else {
        molecule.assembly_status = "explicit_tag_single_qname_group";
      }
      if (group.partial_identity_tag) {
        molecule.warnings.emplace_back("PARTIAL_MOLECULE_ID_TAG_INHERITED");
      }
    }

    const auto apply_optional_protocol_tag =
        [&](const std::string &tag, const std::string &name,
            const std::string &recorded_flag) {
          if (tag.empty()) {
            return;
          }
          const auto values = protocol_tag_values(input.reads, indices, tag);
          molecule.protocol_metadata[name + "_tag"] = tag;
          if (values.values.size() == 1U) {
            molecule.protocol_metadata[name + "_value"] =
                *values.values.begin();
            molecule.protocol_flags.push_back(recorded_flag);
            if (values.missing) {
              molecule.warnings.push_back("PARTIAL_" + name + "_TAG_INHERITED");
            }
          } else if (values.values.empty()) {
            molecule.warnings.push_back("MISSING_" + name + "_TAG");
          } else {
            molecule.ambiguous = true;
            molecule.analysis_eligible = false;
            molecule.warnings.push_back("CONFLICTING_" + name + "_TAG_VALUES");
            molecule.exclusion_reasons.push_back("CONFLICTING_" + name +
                                                 "_TAG_VALUES");
          }
        };
    apply_optional_protocol_tag(config.umi_tag, "UMI", "UMI_RECORDED");
    apply_optional_protocol_tag(config.duplex_tag, "DUPLEX",
                                "DUPLEX_METADATA_RECORDED");

    const bool has_eligible_primary = std::any_of(
        primary_candidates.begin(), primary_candidates.end(),
        [&](const auto index) {
          return (input.reads[index].flags & (0x4U | 0x200U | 0x400U)) == 0U;
        });
    if (!primary_candidates.empty() && !has_eligible_primary) {
      molecule.analysis_eligible = false;
      molecule.warnings.emplace_back("NO_ELIGIBLE_PRIMARY_EVIDENCE");
      molecule.exclusion_reasons.emplace_back("NO_ELIGIBLE_PRIMARY_EVIDENCE");
    }

    if (source_qnames.size() > 1U) {
      molecule.protocol_flags.emplace_back("MULTI_QNAME_MOLECULE");
    }
    std::set<std::string> reference_names;
    bool concatemer_candidate = false;
    bool origin_spanning = false;
    for (const auto index : indices) {
      const auto &record = input.reads[index];
      if ((record.flags & 0x100U) != 0U) {
        molecule.protocol_flags.emplace_back("SECONDARY_EVIDENCE_PRESENT");
      }
      if ((record.flags & 0x800U) != 0U) {
        molecule.protocol_flags.emplace_back("SUPPLEMENTARY_EVIDENCE_PRESENT");
      }
      if ((record.flags & 0x400U) != 0U) {
        molecule.protocol_flags.emplace_back("SAM_DUPLICATE_FLAG_PRESENT");
      }
      if (!record.reference_name.empty() && record.reference_name != "*") {
        reference_names.insert(record.reference_name);
      }
      if (reference_length != 0U &&
          record.sequence.size() / 2U >= reference_length) {
        concatemer_candidate = true;
      }
      const auto span = reference_consuming_length(record);
      if (reference_length != 0U && record.reference_start != 0U &&
          span != 0U &&
          span - 1U > reference_length -
                          std::min(reference_length, record.reference_start)) {
        origin_spanning = true;
      }
    }
    if (reference_names.size() > 1U) {
      molecule.protocol_flags.emplace_back("MULTI_REFERENCE_ALIGNMENT");
    }
    if (concatemer_candidate) {
      molecule.protocol_flags.emplace_back("CONCATEMER_LENGTH_CANDIDATE");
    }
    if (origin_spanning) {
      molecule.protocol_flags.emplace_back("ORIGIN_SPANNING_ALIGNMENT");
    }

    molecule.fragment_ids.reserve(indices.size());
    for (const auto index : indices) {
      const auto &record = input.reads[index];
      molecule.fragment_ids.push_back(AlignmentFragmentId{index});
      result.fragments[index] = {
          AlignmentFragmentId{index},
          molecule_index,
          index,
          molecule_id,
          alignment_fragment_role(record.flags, true),
          index == representative_index,
          record.flags,
          record.mapping_quality,
          record.reference_name,
          record.reference_start,
          record.cigar,
      };
    }

    ReadRecord representative = input.reads[representative_index];
    bool has_sa_tag = representative.aux_tags.contains("SA");
    std::string sa =
        has_sa_tag ? representative.aux_tags.at("SA") : std::string{};
    for (const auto index : indices) {
      const bool explicit_fragment_link =
          !config.molecule_id_tag.empty() &&
          (input.reads[index].flags & 0x100U) == 0U;
      if (index == representative_index ||
          ((input.reads[index].flags & 0x800U) == 0U &&
           !explicit_fragment_link) ||
          input.reads[index].reference_name.empty() ||
          input.reads[index].reference_name == "*") {
        continue;
      }
      const auto encoded = alignment_as_sa(input.reads[index]);
      if (encoded == alignment_as_sa(representative)) {
        continue;
      }
      if (sa.find(encoded) == std::string::npos) {
        sa += encoded;
      }
      has_sa_tag = true;
    }
    if (has_sa_tag) {
      representative.aux_tags["SA"] = std::move(sa);
      representative.supplementary_alignments = parse_supplementary_alignments(
          representative.aux_tags.at("SA"), representative.id);
    }
    representative.aux_tags["MA"] = std::to_string(indices.size());
    representative.id = molecule.id;
    sort_unique_strings(molecule.protocol_flags);
    sort_unique_strings(molecule.exclusion_reasons);
    sort_unique_strings(molecule.warnings);
    result.representatives.push_back(std::move(representative));
    result.molecules.push_back(std::move(molecule));
  }
  result.source_alignment_records = std::move(input.reads);
  return result;
}

[[nodiscard]] double mean_quality(const std::string &qualities) {
  if (qualities.empty()) {
    return 0.0;
  }
  const auto total = std::accumulate(
      qualities.begin(), qualities.end(), 0.0, [](double sum, char c) {
        return sum + static_cast<double>(std::max(0, c - 33));
      });
  return total / static_cast<double>(qualities.size());
}

[[nodiscard]] double gc_fraction(std::string_view sequence) {
  if (sequence.empty()) {
    return 0.0;
  }
  std::size_t gc = 0;
  for (const char c : sequence) {
    const char upper =
        static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (upper == 'G' || upper == 'C') {
      ++gc;
    }
  }
  return static_cast<double>(gc) / static_cast<double>(sequence.size());
}

[[nodiscard]] bool looks_like_nuclear_contig(std::string_view reference_name) {
  std::string name(reference_name);
  std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  if (name.rfind("chr", 0) == 0) {
    name.erase(0, 3);
  }
  if (name == "x" || name == "y") {
    return true;
  }
  const auto chromosome = parse_size(name);
  return chromosome && *chromosome >= 1 && *chromosome <= 22;
}

[[nodiscard]] bool has_nuclear_supplementary_alignment(const ReadRecord &read) {
  return std::any_of(
      read.supplementary_alignments.begin(),
      read.supplementary_alignments.end(), [](const auto &alignment) {
        return looks_like_nuclear_contig(alignment.reference_name);
      });
}

[[nodiscard]] std::vector<std::string>
numt_evidence(const ReadRecord &read, std::size_t reference_length) {
  std::vector<std::string> evidence;
  const std::string id_lower = [&] {
    std::string out = read.id;
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
      return static_cast<char>(std::tolower(c));
    });
    return out;
  }();

  if (id_lower.find("numt") != std::string::npos ||
      id_lower.find("nuclear") != std::string::npos) {
    evidence.emplace_back("read_name_heuristic");
  }
  if (looks_like_nuclear_contig(read.reference_name)) {
    evidence.emplace_back("primary_nuclear_alignment");
  }
  if (has_nuclear_supplementary_alignment(read)) {
    evidence.emplace_back("supplementary_nuclear_alignment");
  }

  const double length_ratio = reference_length == 0
                                  ? 0.0
                                  : static_cast<double>(read.sequence.size()) /
                                        static_cast<double>(reference_length);
  const double gc = gc_fraction(read.sequence);
  if (length_ratio > 1.20) {
    evidence.emplace_back("length_exceeds_mtdna");
  }
  if (gc > 0.62 || gc < 0.25) {
    evidence.emplace_back("atypical_gc_fraction");
  }
  return evidence;
}

[[nodiscard]] double numt_score(const std::vector<std::string> &evidence) {
  double score = 0.05;
  for (const auto &item : evidence) {
    if (item == "primary_nuclear_alignment") {
      score = std::max(score, 0.99);
    } else if (item == "supplementary_nuclear_alignment") {
      score = std::max(score, 0.90);
    } else if (item == "read_name_heuristic") {
      score = std::max(score, 0.82);
    } else if (item == "length_exceeds_mtdna") {
      score += 0.30;
    } else if (item == "atypical_gc_fraction") {
      score += 0.18;
    }
  }
  return std::min(0.99, score);
}

[[nodiscard]] bool is_base(char base) {
  base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
  return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

[[nodiscard]] bool passes_snp_alignment_filters(const ReadRecord &read,
                                                const AnalysisConfig &config) {
  return (read.flags & 0x4U) == 0U &&
         (read.flags & config.excluded_snp_flags) == 0U &&
         read.mapping_quality >= config.min_mapping_quality &&
         read.reference_start != 0 && !read.reference_name.empty() &&
         read.reference_name != "*" &&
         !looks_like_nuclear_contig(read.reference_name);
}

[[nodiscard]] std::optional<std::uint8_t>
phred_quality_at(const ReadRecord &read, std::size_t query_index) {
  if (query_index >= read.sequence.size() ||
      query_index >= read.qualities.size()) {
    return std::nullopt;
  }
  const auto encoded = static_cast<unsigned char>(read.qualities[query_index]);
  if (encoded < 33U || encoded > 126U) {
    return std::nullopt;
  }
  return static_cast<std::uint8_t>(encoded - 33U);
}

[[nodiscard]] std::optional<std::size_t> base_index(char base) {
  switch (static_cast<char>(std::toupper(static_cast<unsigned char>(base)))) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return std::nullopt;
  }
}

[[nodiscard]] std::vector<SnpCall>
call_snps_from_alignment(const ReadRecord &read, const std::string &reference,
                         const AnalysisConfig &config) {
  std::vector<SnpCall> snps;
  if (reference.empty() || read.sequence.empty() ||
      read.cigar_operations.empty() ||
      !passes_snp_alignment_filters(read, config)) {
    return snps;
  }

  std::map<std::size_t, char> molecule_alleles;
  std::size_t read_cursor = 0;
  std::size_t reference_cursor =
      read.reference_start == 0 ? 1 : read.reference_start;

  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;

    if (op == 'M' || op == '=' || op == 'X') {
      for (std::size_t offset = 0;
           offset < len && read_cursor + offset < read.sequence.size();
           ++offset) {
        const std::size_t position =
            ((reference_cursor - 1 + offset) % reference.size()) + 1;
        const char reference_base = reference[position - 1];
        const char alternate_base = static_cast<char>(std::toupper(
            static_cast<unsigned char>(read.sequence[read_cursor + offset])));
        const auto quality = phred_quality_at(read, read_cursor + offset);
        if (quality && *quality >= config.min_base_quality &&
            is_base(alternate_base) && is_base(reference_base)) {
          auto [allele, inserted] =
              molecule_alleles.try_emplace(position, alternate_base);
          if (!inserted && allele->second != alternate_base) {
            allele->second = 'N';
          }
        }
      }
      read_cursor += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'I' || op == 'S') {
      read_cursor += len;
    } else if (op == 'D' || op == 'N') {
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'H' || op == 'P') {
      // Does not consume read or reference sequence.
    }
  }

  snps.reserve(molecule_alleles.size());
  for (const auto &[position, alternate] : molecule_alleles) {
    const char reference_base = reference[position - 1U];
    if (is_base(alternate) && alternate != reference_base) {
      SnpCall call;
      call.position = position;
      call.reference = reference_base;
      call.alternate = alternate;
      snps.push_back(std::move(call));
    }
  }

  return snps;
}

[[nodiscard]] std::vector<SnpCall>
call_snps_from_reference_span(const ReadRecord &read,
                              const std::string &reference,
                              const AnalysisConfig &config) {
  std::vector<SnpCall> snps;
  if (reference.empty() || read.sequence.empty() || read.reference_start == 0 ||
      read.reference_name.empty() ||
      !passes_snp_alignment_filters(read, config)) {
    return snps;
  }

  for (std::size_t i = 0; i < read.sequence.size(); ++i) {
    const std::size_t position =
        ((read.reference_start - 1 + i) % reference.size()) + 1;
    const char reference_base = reference[position - 1];
    const char alternate_base = static_cast<char>(
        std::toupper(static_cast<unsigned char>(read.sequence[i])));
    const auto quality = phred_quality_at(read, i);
    if (quality && *quality >= config.min_base_quality &&
        is_base(alternate_base) && is_base(reference_base) &&
        alternate_base != reference_base) {
      SnpCall call;
      call.position = position;
      call.reference = reference_base;
      call.alternate = alternate_base;
      snps.push_back(std::move(call));
    }
  }

  return snps;
}

[[nodiscard]] std::optional<SnpCall> parse_snp_tag(std::string_view token) {
  const auto colon = token.find(':');
  const auto arrow =
      token.find('>', colon == std::string_view::npos ? 0 : colon + 1);
  if (colon == std::string_view::npos || arrow == std::string_view::npos ||
      arrow + 1 >= token.size()) {
    return std::nullopt;
  }

  const auto position = parse_size(token.substr(0, colon));
  if (!position || *position == 0) {
    return std::nullopt;
  }

  const char reference = static_cast<char>(
      std::toupper(static_cast<unsigned char>(token[colon + 1])));
  const char alternate = static_cast<char>(
      std::toupper(static_cast<unsigned char>(token[arrow + 1])));
  const auto valid_base = [](char base) {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
  };
  if (!valid_base(reference) || !valid_base(alternate) ||
      reference == alternate) {
    return std::nullopt;
  }

  SnpCall call;
  call.position = *position;
  call.reference = reference;
  call.alternate = alternate;
  return call;
}

[[nodiscard]] std::vector<SnpCall>
snp_tags_from_header(const ReadRecord &read) {
  std::vector<SnpCall> calls;
  std::istringstream iss(read.id);
  std::string token;
  while (iss >> token) {
    constexpr std::string_view prefix = "snp=";
    if (token.rfind(prefix, 0) == 0) {
      auto call = parse_snp_tag(std::string_view(token).substr(prefix.size()));
      if (call) {
        calls.push_back(std::move(*call));
      }
    }
  }
  return calls;
}

void merge_snps(std::vector<SnpCall> &snps, std::vector<SnpCall> tagged_snps) {
  snps.insert(snps.end(), std::make_move_iterator(tagged_snps.begin()),
              std::make_move_iterator(tagged_snps.end()));
  std::sort(snps.begin(), snps.end(), [](const auto &lhs, const auto &rhs) {
    return snp_key(lhs.position, lhs.reference, lhs.alternate) <
           snp_key(rhs.position, rhs.reference, rhs.alternate);
  });
  snps.erase(std::unique(snps.begin(), snps.end(),
                         [](const auto &lhs, const auto &rhs) {
                           return lhs.position == rhs.position &&
                                  lhs.reference == rhs.reference &&
                                  lhs.alternate == rhs.alternate;
                         }),
             snps.end());
}

[[nodiscard]] bool
observed_phylo_mutation_less(const ObservedPhyloMutation &lhs,
                             const ObservedPhyloMutation &rhs) {
  return std::tie(lhs.position, lhs.kind, lhs.insertion_index, lhs.encoded) <
         std::tie(rhs.position, rhs.kind, rhs.insertion_index, rhs.encoded);
}

[[nodiscard]] ObservedPhyloMutation
observed_phylo_mutation(const PhyloMutation &mutation) {
  ObservedPhyloMutation observed;
  observed.position = mutation.locus.position;
  observed.kind = mutation.kind;
  observed.alternate = mutation.alternate;
  observed.insertion_index = mutation.locus.insertion_index;
  observed.inserted_bases = mutation.inserted_bases;
  observed.encoded = mutation.encoded;
  return observed;
}

[[nodiscard]] std::vector<std::string>
split_alignment_rule_tokens(std::string_view value) {
  std::vector<std::string> tokens;
  std::size_t cursor = 0U;
  while (cursor < value.size()) {
    while (cursor < value.size() &&
           (std::isspace(static_cast<unsigned char>(value[cursor])) != 0 ||
            value[cursor] == ',')) {
      ++cursor;
    }
    const std::size_t start = cursor;
    while (cursor < value.size() &&
           std::isspace(static_cast<unsigned char>(value[cursor])) == 0 &&
           value[cursor] != ',') {
      ++cursor;
    }
    if (start != cursor) {
      tokens.emplace_back(value.substr(start, cursor - start));
    }
  }
  return tokens;
}

[[nodiscard]] bool is_reference_restore_token(std::string_view token) {
  return !token.empty() &&
         std::all_of(token.begin(), token.end(), [](const char value) {
           return std::isdigit(static_cast<unsigned char>(value)) != 0;
         });
}

[[nodiscard]] std::vector<PhyloAlignmentRule> load_phylo_alignment_rules() {
  const auto path = phylotree_alignment_rules_path();
  std::ifstream input(path);
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open PhyloTree alignment rules: " + path);
  }

  std::vector<PhyloAlignmentRule> rules;
  std::string line;
  std::size_t line_number = 0U;
  while (std::getline(input, line)) {
    ++line_number;
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line_number == 1U) {
      if (line != "error,expected") {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "PhyloTree alignment-rule header is invalid");
      }
      continue;
    }
    if (trim_copy(line).empty()) {
      continue;
    }
    const auto delimiter = line.find(',');
    if (delimiter == std::string::npos) {
      throw AnalysisError(AnalysisErrorCode::resource_invalid,
                          "PhyloTree alignment rule has no delimiter at line " +
                              std::to_string(line_number));
    }

    PhyloAlignmentRule rule;
    bool callable_error_side = true;
    for (const auto &token : split_alignment_rule_tokens(
             std::string_view(line).substr(0U, delimiter))) {
      const auto mutation = parse_phylo_mutation(token);
      if (!mutation || mutation->backmutation) {
        // Ambiguity-code rules such as 3106N cannot be emitted by the
        // fail-closed callable-base path and are intentionally inapplicable.
        callable_error_side = false;
        break;
      }
      rule.erroneous_markers.push_back(mutation->encoded);
    }
    if (!callable_error_side) {
      continue;
    }
    if (rule.erroneous_markers.empty()) {
      throw AnalysisError(
          AnalysisErrorCode::resource_invalid,
          "PhyloTree alignment rule has an empty error side at line " +
              std::to_string(line_number));
    }
    std::sort(rule.erroneous_markers.begin(), rule.erroneous_markers.end());
    rule.erroneous_markers.erase(std::unique(rule.erroneous_markers.begin(),
                                             rule.erroneous_markers.end()),
                                 rule.erroneous_markers.end());

    for (const auto &token : split_alignment_rule_tokens(
             std::string_view(line).substr(delimiter + 1U))) {
      if (is_reference_restore_token(token)) {
        continue;
      }
      const auto mutation = parse_phylo_mutation(token);
      if (!mutation || mutation->backmutation) {
        throw AnalysisError(
            AnalysisErrorCode::resource_invalid,
            "unsupported PhyloTree alignment-rule replacement '" + token +
                "' at line " + std::to_string(line_number));
      }
      rule.replacement_markers.push_back(observed_phylo_mutation(*mutation));
    }
    std::sort(rule.replacement_markers.begin(), rule.replacement_markers.end(),
              observed_phylo_mutation_less);
    rule.replacement_markers.erase(
        std::unique(rule.replacement_markers.begin(),
                    rule.replacement_markers.end(),
                    [](const auto &lhs, const auto &rhs) {
                      return lhs.encoded == rhs.encoded;
                    }),
        rule.replacement_markers.end());
    rules.push_back(std::move(rule));
  }
  if (rules.size() < 100U) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "PhyloTree alignment-rule resource is incomplete");
  }
  return rules;
}

void apply_phylo_alignment_rules(std::vector<ObservedPhyloMutation> &observed,
                                 const std::vector<PhyloAlignmentRule> &rules) {
  std::map<std::string, ObservedPhyloMutation> markers;
  for (auto &mutation : observed) {
    markers.insert_or_assign(mutation.encoded, std::move(mutation));
  }

  for (std::size_t pass = 0U; pass <= rules.size(); ++pass) {
    bool changed = false;
    for (const auto &rule : rules) {
      if (!std::all_of(
              rule.erroneous_markers.begin(), rule.erroneous_markers.end(),
              [&](const auto &marker) { return markers.contains(marker); })) {
        continue;
      }
      std::map<std::string, ObservedPhyloMutation> candidate = markers;
      for (const auto &marker : rule.erroneous_markers) {
        candidate.erase(marker);
      }
      for (const auto &replacement : rule.replacement_markers) {
        candidate.insert_or_assign(replacement.encoded, replacement);
      }
      const bool same_keys =
          candidate.size() == markers.size() &&
          std::equal(candidate.begin(), candidate.end(), markers.begin(),
                     [](const auto &lhs, const auto &rhs) {
                       return lhs.first == rhs.first;
                     });
      if (!same_keys) {
        markers = std::move(candidate);
        changed = true;
      }
    }
    if (!changed) {
      observed.clear();
      observed.reserve(markers.size());
      for (auto &[_, marker] : markers) {
        observed.push_back(std::move(marker));
      }
      std::sort(observed.begin(), observed.end(), observed_phylo_mutation_less);
      return;
    }
  }
  throw AnalysisError(AnalysisErrorCode::resource_invalid,
                      "PhyloTree alignment rules did not converge");
}

[[nodiscard]] std::vector<ObservedPhyloMutation>
phylo_tags_from_header(const ReadRecord &read) {
  std::vector<ObservedPhyloMutation> calls;
  std::istringstream iss(read.id);
  std::string token;
  while (iss >> token) {
    constexpr std::string_view prefix = "phylo=";
    if (token.rfind(prefix, 0) != 0) {
      continue;
    }
    const auto mutation =
        parse_phylo_mutation(std::string_view(token).substr(prefix.size()));
    if (!mutation || mutation->backmutation) {
      continue;
    }
    calls.push_back(observed_phylo_mutation(*mutation));
  }
  return calls;
}

[[nodiscard]] std::vector<ObservedPhyloMutation>
call_phylo_indels_from_alignment(const ReadRecord &read,
                                 const AnalysisConfig &config,
                                 const std::string &reference) {
  std::vector<ObservedPhyloMutation> calls;
  const std::size_t reference_length = reference.size();
  if (reference_length == 0 || read.sequence.empty() ||
      read.reference_start == 0U || read.cigar_operations.empty() ||
      !passes_snp_alignment_filters(read, config)) {
    return calls;
  }

  constexpr std::size_t kMaxPhylogeneticIndelLength = 50U;
  std::size_t read_cursor = 0;
  std::size_t reference_cursor =
      ((read.reference_start - 1U) % reference_length) + 1U;
  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      read_cursor += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference_length)) %
                          reference_length) +
                         1U;
    } else if (op == 'I') {
      if (len <= kMaxPhylogeneticIndelLength &&
          read_cursor <= read.sequence.size() &&
          len <= read.sequence.size() - read_cursor) {
        std::string inserted;
        inserted.reserve(len);
        bool callable = true;
        for (std::size_t offset = 0; offset < len; ++offset) {
          const char base = static_cast<char>(std::toupper(
              static_cast<unsigned char>(read.sequence[read_cursor + offset])));
          const auto quality = phred_quality_at(read, read_cursor + offset);
          if (!is_base(base) || !quality ||
              *quality < config.min_base_quality) {
            callable = false;
            break;
          }
          inserted.push_back(base);
        }
        if (callable) {
          std::size_t anchor =
              reference_cursor > 1U ? reference_cursor - 1U : reference_length;
          // Mitochondrial variant nomenclature places repeat-equivalent indels
          // at the most 3-prime coordinate. Rotate the inserted sequence while
          // shifting so that non-homopolymer repeats remain
          // haplotype-equivalent. Do not cross the rCRS coordinate boundary:
          // circular-origin events have their own canonical representation and
          // no unique 3-prime endpoint.
          std::size_t shift = 0U;
          while (anchor < reference_length && !inserted.empty() &&
                 inserted[shift % inserted.size()] == reference[anchor]) {
            ++shift;
            ++anchor;
          }
          if (!inserted.empty() && shift % inserted.size() != 0U) {
            std::rotate(inserted.begin(),
                        inserted.begin() + static_cast<std::ptrdiff_t>(
                                               shift % inserted.size()),
                        inserted.end());
          }
          ObservedPhyloMutation mutation;
          mutation.position = anchor;
          mutation.kind = PhyloMutationKind::insertion;
          mutation.insertion_index = "1";
          mutation.inserted_bases = std::move(inserted);
          mutation.encoded = std::to_string(mutation.position) + ".1" +
                             mutation.inserted_bases;
          calls.push_back(std::move(mutation));
        }
      }
      read_cursor += len;
    } else if (op == 'D') {
      if (len <= kMaxPhylogeneticIndelLength) {
        std::size_t normalized_start = reference_cursor;
        const bool crosses_origin =
            len > reference_length - reference_cursor + 1U;
        while (!crosses_origin && normalized_start <= reference_length &&
               len <= reference_length - normalized_start &&
               reference[normalized_start - 1U] ==
                   reference[normalized_start + len - 1U]) {
          ++normalized_start;
        }
        for (std::size_t offset = 0; offset < len; ++offset) {
          ObservedPhyloMutation mutation;
          mutation.position =
              ((normalized_start - 1U + offset) % reference_length) + 1U;
          mutation.kind = PhyloMutationKind::deletion;
          mutation.encoded = std::to_string(mutation.position) + "d";
          calls.push_back(std::move(mutation));
        }
      }
      reference_cursor = ((reference_cursor - 1U + (len % reference_length)) %
                          reference_length) +
                         1U;
    } else if (op == 'N') {
      reference_cursor = ((reference_cursor - 1U + (len % reference_length)) %
                          reference_length) +
                         1U;
    } else if (op == 'S') {
      read_cursor += len;
    } else if (op == 'H' || op == 'P') {
      // Does not consume query or reference sequence.
    }
  }
  return calls;
}

[[nodiscard]] std::vector<NormalizedSmallIndel>
call_small_indels_from_alignment(const ReadRecord &read,
                                 const AnalysisConfig &config,
                                 const std::string &reference) {
  std::vector<NormalizedSmallIndel> calls;
  if (reference.empty() || read.reference_start == 0U ||
      read.cigar_operations.empty() ||
      !passes_snp_alignment_filters(read, config) ||
      config.sv_min_length <= 1U) {
    return calls;
  }

  constexpr std::size_t kMaximumSmallIndelLength = 50U;
  const auto maximum_length =
      std::min(kMaximumSmallIndelLength, config.sv_min_length - 1U);
  std::size_t query_cursor = 0U;
  std::size_t reference_cursor =
      ((read.reference_start - 1U) % reference.size()) + 1U;
  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      query_cursor += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
      continue;
    }
    if (op == 'I') {
      if (len <= maximum_length && query_cursor <= read.sequence.size() &&
          len <= read.sequence.size() - query_cursor) {
        std::string inserted;
        inserted.reserve(len);
        std::optional<std::uint8_t> minimum_quality;
        bool callable = true;
        for (std::size_t offset = 0U; offset < len; ++offset) {
          const char base =
              static_cast<char>(std::toupper(static_cast<unsigned char>(
                  read.sequence[query_cursor + offset])));
          const auto quality = phred_quality_at(read, query_cursor + offset);
          if (!is_base(base) || !quality ||
              *quality < config.min_base_quality) {
            callable = false;
            break;
          }
          inserted.push_back(base);
          minimum_quality =
              minimum_quality ? std::min(*minimum_quality, *quality) : *quality;
        }
        if (callable) {
          std::size_t anchor =
              reference_cursor > 1U ? reference_cursor - 1U : reference.size();
          std::size_t shift = 0U;
          while (anchor < reference.size() && !inserted.empty() &&
                 inserted[shift % inserted.size()] == reference[anchor]) {
            ++shift;
            ++anchor;
          }
          if (!inserted.empty() && shift % inserted.size() != 0U) {
            std::rotate(inserted.begin(),
                        inserted.begin() + static_cast<std::ptrdiff_t>(
                                               shift % inserted.size()),
                        inserted.end());
          }
          NormalizedSmallIndel event;
          event.type = "SMALL_INSERTION";
          event.start = anchor;
          event.end = anchor;
          event.length = len;
          event.reference.assign(1U, reference[anchor - 1U]);
          event.alternate = event.reference + inserted;
          event.normalization = "rcrs_3prime_small_indel_v1";
          event.base_quality = minimum_quality;
          event.id =
              "indel:insertion:" + std::to_string(anchor) + ":" + inserted;
          calls.push_back(std::move(event));
        }
      }
      query_cursor += len;
      continue;
    }
    if (op == 'D') {
      if (len <= maximum_length && len < reference.size()) {
        std::size_t normalized_start = reference_cursor;
        const bool crosses_origin =
            reference_cursor == 1U ||
            len > reference.size() - reference_cursor + 1U;
        while (!crosses_origin && normalized_start <= reference.size() &&
               len <= reference.size() - normalized_start &&
               reference[normalized_start - 1U] ==
                   reference[normalized_start + len - 1U]) {
          ++normalized_start;
        }
        const auto normalized_end =
            ((normalized_start - 1U + (len - 1U)) % reference.size()) + 1U;
        const auto anchor =
            normalized_start > 1U ? normalized_start - 1U : reference.size();
        std::string deleted;
        deleted.reserve(len);
        for (std::size_t offset = 0U; offset < len; ++offset) {
          deleted.push_back(
              reference[(normalized_start - 1U + offset) % reference.size()]);
        }
        NormalizedSmallIndel event;
        event.type = "SMALL_DELETION";
        event.start = normalized_start;
        event.end = normalized_end;
        event.length = len;
        event.reference.assign(1U, reference[anchor - 1U]);
        event.reference += deleted;
        event.alternate.assign(1U, reference[anchor - 1U]);
        event.normalization = crosses_origin ? "rcrs_circular_small_indel_v1"
                                             : "rcrs_3prime_small_indel_v1";
        event.id = "indel:deletion:" + std::to_string(normalized_start) + "-" +
                   std::to_string(normalized_end) + ":" + deleted;
        calls.push_back(std::move(event));
      }
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
      continue;
    }
    if (op == 'N') {
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'S') {
      query_cursor += len;
    }
  }
  std::sort(calls.begin(), calls.end(),
            [](const auto &lhs, const auto &rhs) { return lhs.id < rhs.id; });
  calls.erase(std::unique(calls.begin(), calls.end(),
                          [](const auto &lhs, const auto &rhs) {
                            return lhs.id == rhs.id;
                          }),
              calls.end());
  return calls;
}

void merge_haplogroup_markers(std::vector<ObservedPhyloMutation> &markers,
                              std::vector<ObservedPhyloMutation> additional) {
  markers.insert(markers.end(), std::make_move_iterator(additional.begin()),
                 std::make_move_iterator(additional.end()));
  std::sort(markers.begin(), markers.end(), observed_phylo_mutation_less);
  markers.erase(std::unique(markers.begin(), markers.end(),
                            [](const auto &lhs, const auto &rhs) {
                              return lhs.encoded == rhs.encoded;
                            }),
                markers.end());
}

[[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>>
ranges_from_bitmap(const std::vector<bool> &positions) {
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  std::size_t position = 1U;
  while (position < positions.size()) {
    while (position < positions.size() && !positions[position]) {
      ++position;
    }
    if (position == positions.size()) {
      break;
    }
    const std::size_t start = position;
    while (position + 1U < positions.size() && positions[position + 1U]) {
      ++position;
    }
    ranges.emplace_back(start, position);
    ++position;
  }
  return ranges;
}

[[nodiscard]] std::optional<std::vector<std::pair<std::size_t, std::size_t>>>
phylo_ranges_from_header(const ReadRecord &read) {
  std::istringstream iss(read.id);
  std::string token;
  while (iss >> token) {
    constexpr std::string_view prefix = "phylo_range=";
    if (token.rfind(prefix, 0) != 0) {
      continue;
    }
    const std::string_view encoded =
        std::string_view(token).substr(prefix.size());
    if (encoded.empty()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "empty development phylo_range for read '" + read.id +
                              "'");
    }
    std::vector<bool> positions(
        static_cast<std::size_t>(kDefaultReferenceLength) + 1U, false);
    std::size_t item_start = 0U;
    while (item_start < encoded.size()) {
      const auto item_end = encoded.find(',', item_start);
      const auto item =
          encoded.substr(item_start, item_end == std::string_view::npos
                                         ? encoded.size() - item_start
                                         : item_end - item_start);
      const auto dash = item.find('-');
      const auto start = parse_size(item.substr(0, dash));
      const auto end = dash == std::string_view::npos
                           ? start
                           : parse_size(item.substr(dash + 1U));
      if (!start || !end || *start == 0U || *end < *start ||
          *end > static_cast<std::size_t>(kDefaultReferenceLength)) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "invalid development phylo_range for read '" +
                                read.id + "'");
      }
      for (std::size_t position = *start; position <= *end; ++position) {
        positions[position] = true;
      }
      if (item_end == std::string_view::npos) {
        break;
      }
      item_start = item_end + 1U;
    }
    return ranges_from_bitmap(positions);
  }
  return std::nullopt;
}

struct CallablePhyloRanges {
  bool known = false;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
};

void append_callable_position(
    std::vector<std::pair<std::size_t, std::size_t>> &ranges,
    std::optional<std::pair<std::size_t, std::size_t>> &current,
    std::size_t position) {
  if (!current) {
    current = std::pair{position, position};
    return;
  }
  if (position == current->second) {
    return;
  }
  if (current->second < std::numeric_limits<std::size_t>::max() &&
      position == current->second + 1U) {
    current->second = position;
    return;
  }
  ranges.push_back(*current);
  current = std::pair{position, position};
}

[[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>>
normalize_ranges(std::vector<std::pair<std::size_t, std::size_t>> ranges) {
  if (ranges.empty()) {
    return ranges;
  }
  std::sort(ranges.begin(), ranges.end());
  std::vector<std::pair<std::size_t, std::size_t>> normalized;
  normalized.reserve(ranges.size());
  for (const auto &range : ranges) {
    if (normalized.empty() ||
        (normalized.back().second < std::numeric_limits<std::size_t>::max() &&
         range.first > normalized.back().second + 1U)) {
      normalized.push_back(range);
      continue;
    }
    normalized.back().second = std::max(normalized.back().second, range.second);
  }
  return normalized;
}

[[nodiscard]] CallablePhyloRanges
callable_phylo_ranges_from_alignment(const ReadRecord &read,
                                     const std::string &reference,
                                     const AnalysisConfig &config) {
  CallablePhyloRanges result;
  if (reference.empty() || read.reference_start == 0U ||
      read.cigar_operations.empty() ||
      !passes_snp_alignment_filters(read, config)) {
    return result;
  }
  result.known = true;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  std::optional<std::pair<std::size_t, std::size_t>> current;
  std::size_t read_cursor = 0U;
  std::size_t reference_cursor =
      ((read.reference_start - 1U) % reference.size()) + 1U;
  constexpr std::size_t kMaxPhylogeneticIndelLength = 50U;
  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      for (std::size_t offset = 0;
           offset < len && read_cursor + offset < read.sequence.size();
           ++offset) {
        const std::size_t position =
            ((reference_cursor - 1U + offset) % reference.size()) + 1U;
        const char base = static_cast<char>(std::toupper(
            static_cast<unsigned char>(read.sequence[read_cursor + offset])));
        const auto quality = phred_quality_at(read, read_cursor + offset);
        if (is_base(base) && is_base(reference[position - 1U]) && quality &&
            *quality >= config.min_base_quality) {
          append_callable_position(ranges, current, position);
        }
      }
      read_cursor += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'I' || op == 'S') {
      read_cursor += len;
    } else if (op == 'D') {
      if (len <= kMaxPhylogeneticIndelLength) {
        for (std::size_t offset = 0; offset < len; ++offset) {
          append_callable_position(
              ranges, current,
              ((reference_cursor - 1U + offset) % reference.size()) + 1U);
        }
      }
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'N') {
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'H' || op == 'P') {
      // Does not consume query or reference sequence.
    }
  }
  if (current) {
    ranges.push_back(*current);
  }
  result.ranges = normalize_ranges(std::move(ranges));
  return result;
}

[[nodiscard]] const ReadRecord &
source_record_for_fragment(const MoleculeAssemblyResult &assembly,
                           const AlignmentFragmentEvidence &fragment) {
  if (!assembly.source_alignment_records.empty()) {
    if (fragment.source_record_index >=
        assembly.source_alignment_records.size()) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "alignment fragment source record is unresolved");
    }
    return assembly.source_alignment_records[fragment.source_record_index];
  }
  if (fragment.molecule_index.value >= assembly.representatives.size()) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "FASTQ fragment source record is unresolved");
  }
  return assembly.representatives[fragment.molecule_index.value];
}

[[nodiscard]] std::string
alignment_callability_exclusion(const ReadRecord &read,
                                const AnalysisConfig &config) {
  if ((read.flags & 0x4U) != 0U) {
    return "UNMAPPED_ALIGNMENT";
  }
  if ((read.flags & 0x100U) != 0U) {
    return "SECONDARY_ALIGNMENT";
  }
  if ((read.flags & 0x200U) != 0U) {
    return "QC_FAILED_ALIGNMENT";
  }
  if ((read.flags & 0x400U) != 0U) {
    return "DUPLICATE_ALIGNMENT";
  }
  const auto configured_exclusions =
      static_cast<std::uint16_t>(config.excluded_snp_flags & ~0x800U);
  if ((read.flags & configured_exclusions) != 0U) {
    return "EXCLUDED_ALIGNMENT_FLAGS";
  }
  if (read.mapping_quality < config.min_mapping_quality) {
    return "LOW_MAPPING_QUALITY";
  }
  if (read.reference_start == 0U || read.reference_name.empty() ||
      read.reference_name == "*") {
    return "MISSING_REFERENCE_PLACEMENT";
  }
  if (looks_like_nuclear_contig(read.reference_name)) {
    return "NON_MITOCHONDRIAL_ALIGNMENT";
  }
  if (read.cigar_operations.empty()) {
    return "MISSING_CIGAR";
  }
  if (read.sequence.empty() || read.qualities.empty()) {
    return "MISSING_QUERY_OR_QUALITIES";
  }
  return {};
}

[[nodiscard]] AlignmentCallabilityEvidence
alignment_callability(const AlignmentFragmentEvidence &fragment,
                      const ReadRecord &read, const std::string &reference,
                      const AnalysisConfig &config) {
  AlignmentCallabilityEvidence evidence;
  evidence.alignment_fragment_id = fragment.id;
  const auto excluded = alignment_callability_exclusion(read, config);
  if (!excluded.empty()) {
    evidence.status = excluded;
    return evidence;
  }
  if (reference.empty()) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "callability requires a non-empty reference");
  }

  evidence.eligible = true;
  evidence.status = "ASSESSED";
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  std::optional<std::pair<std::size_t, std::size_t>> current;
  std::size_t query_cursor = 0U;
  std::size_t reference_cursor =
      ((read.reference_start - 1U) % reference.size()) + 1U;
  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      for (std::size_t offset = 0U; offset < len; ++offset) {
        const auto position =
            ((reference_cursor - 1U + (offset % reference.size())) %
             reference.size()) +
            1U;
        if (query_cursor > std::numeric_limits<std::size_t>::max() - offset) {
          throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                              "CIGAR query coordinate overflow for read '" +
                                  read.id + "'");
        }
        const auto query_index = query_cursor + offset;
        if (query_index >= read.sequence.size() ||
            query_index >= read.qualities.size()) {
          ++evidence
                .reference_exclusion_counts["MISSING_QUERY_BASE_OR_QUALITY"];
          continue;
        }
        const char observed = static_cast<char>(std::toupper(
            static_cast<unsigned char>(read.sequence[query_index])));
        const auto quality = phred_quality_at(read, query_index);
        if (!is_base(reference[position - 1U])) {
          ++evidence.reference_exclusion_counts["NON_CANONICAL_REFERENCE_BASE"];
        } else if (!is_base(observed)) {
          ++evidence.reference_exclusion_counts["NON_CANONICAL_OBSERVED_BASE"];
        } else if (!quality || *quality < config.min_base_quality) {
          ++evidence.reference_exclusion_counts["LOW_OR_INVALID_BASE_QUALITY"];
        } else {
          append_callable_position(ranges, current, position);
        }
      }
      if (len > std::numeric_limits<std::size_t>::max() - query_cursor) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "CIGAR query coordinate overflow for read '" +
                                read.id + "'");
      }
      query_cursor += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'I' || op == 'S') {
      if (len > std::numeric_limits<std::size_t>::max() - query_cursor) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "CIGAR query coordinate overflow for read '" +
                                read.id + "'");
      }
      query_cursor += len;
      if (op == 'I') {
        if (len > std::numeric_limits<std::size_t>::max() -
                      evidence.inserted_query_bases) {
          throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                              "inserted query-base count overflow for read '" +
                                  read.id + "'");
        }
        evidence.inserted_query_bases += len;
        evidence.disrupted_adjacency_anchors.push_back(
            reference_cursor > 1U ? reference_cursor - 1U : reference.size());
      } else {
        if (len > std::numeric_limits<std::size_t>::max() -
                      evidence.soft_clipped_query_bases) {
          throw AnalysisError(
              AnalysisErrorCode::input_parse_failed,
              "soft-clipped query-base count overflow for read '" + read.id +
                  "'");
        }
        evidence.soft_clipped_query_bases += len;
      }
    } else if (op == 'D' || op == 'N') {
      auto &excluded_count =
          evidence.reference_exclusion_counts
              [op == 'D' ? "DELETION_FROM_REFERENCE_PATH" : "REFERENCE_SKIP"];
      if (len > std::numeric_limits<std::size_t>::max() - excluded_count) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "reference exclusion count overflow for read '" +
                                read.id + "'");
      }
      excluded_count += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'H' || op == 'P') {
      // Neither query nor reference is consumed.
    }
  }
  if (current) {
    ranges.push_back(*current);
  }
  evidence.ranges = normalize_ranges(std::move(ranges));
  std::sort(evidence.disrupted_adjacency_anchors.begin(),
            evidence.disrupted_adjacency_anchors.end());
  evidence.disrupted_adjacency_anchors.erase(
      std::unique(evidence.disrupted_adjacency_anchors.begin(),
                  evidence.disrupted_adjacency_anchors.end()),
      evidence.disrupted_adjacency_anchors.end());
  for (const auto &[start, end] : evidence.ranges) {
    evidence.callable_bases += end - start + 1U;
  }
  return evidence;
}

[[nodiscard]] CallabilityResult
build_callability(const MoleculeAssemblyResult &assembly,
                  const std::vector<ReadFeature> &features,
                  const std::string &reference, const AnalysisConfig &config) {
  CallabilityResult result;
  result.molecules.reserve(assembly.molecules.size());
  for (std::size_t molecule_index = 0U;
       molecule_index < assembly.molecules.size(); ++molecule_index) {
    throw_if_cancelled(config);
    const auto &molecule = assembly.molecules[molecule_index];
    MoleculeCallabilityEvidence summary;
    summary.molecule_index = molecule.index;
    if (molecule_index >= features.size()) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "molecule callability feature reference is unresolved");
    }
    if (features[molecule_index].filtered_numt) {
      summary.status = "EXCLUDED_NUMT";
      result.molecules.push_back(std::move(summary));
      continue;
    }

    std::vector<std::pair<std::size_t, std::size_t>> molecule_ranges;
    summary.alignments.reserve(molecule.fragment_ids.size());
    for (const auto fragment_id : molecule.fragment_ids) {
      if (fragment_id.value >= assembly.fragments.size()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "molecule callability fragment reference is unresolved");
      }
      const auto &fragment = assembly.fragments[fragment_id.value];
      const auto &read = source_record_for_fragment(assembly, fragment);
      auto alignment = alignment_callability(fragment, read, reference, config);
      if (alignment.eligible) {
        summary.known = true;
        molecule_ranges.insert(molecule_ranges.end(), alignment.ranges.begin(),
                               alignment.ranges.end());
      }
      summary.alignments.push_back(std::move(alignment));
    }
    summary.ranges = normalize_ranges(std::move(molecule_ranges));
    for (const auto &[start, end] : summary.ranges) {
      summary.callable_bases += end - start + 1U;
    }
    summary.status =
        summary.known
            ? (summary.callable_bases == 0U ? "ASSESSED_EMPTY" : "ASSESSED")
            : "NOT_ASSESSABLE";
    result.molecules.push_back(std::move(summary));
  }
  return result;
}

[[nodiscard]] bool callable_position_in_ranges(
    const std::vector<std::pair<std::size_t, std::size_t>> &ranges,
    const std::size_t position) {
  const auto it =
      std::upper_bound(ranges.begin(), ranges.end(), position,
                       [](const std::size_t value, const auto &range) {
                         return value < range.first;
                       });
  return it != ranges.begin() && position <= std::prev(it)->second;
}

[[nodiscard]] bool
span_in_ranges(const std::vector<std::pair<std::size_t, std::size_t>> &ranges,
               const std::size_t start, const std::size_t end) {
  if (start == 0U || end < start) {
    return false;
  }
  const auto it =
      std::upper_bound(ranges.begin(), ranges.end(), start,
                       [](const std::size_t value, const auto &range) {
                         return value < range.first;
                       });
  return it != ranges.begin() && std::prev(it)->second >= end;
}

[[nodiscard]] std::size_t
circular_previous_position(std::size_t position, std::size_t reference_length);
[[nodiscard]] std::size_t circular_next_position(std::size_t position,
                                                 std::size_t reference_length);

[[nodiscard]] bool
reference_adjacency_callable(const MoleculeCallabilityEvidence &callability,
                             const std::size_t left, const std::size_t right) {
  for (const auto &alignment : callability.alignments) {
    if (!alignment.eligible ||
        !callable_position_in_ranges(alignment.ranges, left) ||
        !callable_position_in_ranges(alignment.ranges, right)) {
      continue;
    }
    if ((right == left + 1U || right < left) &&
        !std::binary_search(alignment.disrupted_adjacency_anchors.begin(),
                            alignment.disrupted_adjacency_anchors.end(),
                            left)) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] bool reference_span_with_flanks_callable(
    const MoleculeCallabilityEvidence &callability, const std::size_t start,
    const std::size_t end, const std::size_t reference_length) {
  if (start == 0U || end < start || end > reference_length) {
    return false;
  }
  const auto left = circular_previous_position(start, reference_length);
  const auto right = circular_next_position(end, reference_length);
  for (const auto &alignment : callability.alignments) {
    if (alignment.eligible && span_in_ranges(alignment.ranges, start, end) &&
        callable_position_in_ranges(alignment.ranges, left) &&
        callable_position_in_ranges(alignment.ranges, right)) {
      return true;
    }
  }
  return false;
}

void apply_clinical_annotations(
    std::vector<SnpCall> &snps,
    const std::map<std::string, SnpCall> &annotations) {
  for (auto &snp : snps) {
    const auto it =
        annotations.find(snp_key(snp.position, snp.reference, snp.alternate));
    if (it == annotations.end()) {
      continue;
    }
    const auto position = snp.position;
    const auto reference = snp.reference;
    const auto alternate = snp.alternate;
    snp = it->second;
    snp.position = position;
    snp.reference = reference;
    snp.alternate = alternate;
  }
}

class JsonBuffer {
public:
  explicit JsonBuffer(const std::size_t reserve_bytes) {
    buffer_.reserve(reserve_bytes);
  }

  JsonBuffer &operator<<(const char *value) {
    buffer_.append(value);
    return *this;
  }

  JsonBuffer &operator<<(const std::string &value) {
    buffer_.append(value);
    return *this;
  }

  JsonBuffer &operator<<(const std::string_view value) {
    buffer_.append(value);
    return *this;
  }

  JsonBuffer &operator<<(const char value) {
    buffer_.push_back(value);
    return *this;
  }

  JsonBuffer &operator<<(const bool value) {
    buffer_.append(value ? "true" : "false");
    return *this;
  }

  template <typename Integer>
    requires(std::is_integral_v<Integer> && !std::is_same_v<Integer, bool> &&
             !std::is_same_v<Integer, char>)
  JsonBuffer &operator<<(const Integer value) {
    std::array<char, 32> encoded{};
    const auto [end, error] =
        std::to_chars(encoded.data(), encoded.data() + encoded.size(), value);
    if (error != std::errc{}) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "integer JSON serialization failed");
    }
    buffer_.append(encoded.data(), end);
    return *this;
  }

  JsonBuffer &operator<<(const double value) {
    std::array<char, 64> encoded{};
    const auto [end, error] =
        std::to_chars(encoded.data(), encoded.data() + encoded.size(), value,
                      std::chars_format::fixed, 9);
    if (error != std::errc{}) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "floating-point JSON serialization failed");
    }
    buffer_.append(encoded.data(), end);
    return *this;
  }

  [[nodiscard]] std::string take() && { return std::move(buffer_); }

private:
  std::string buffer_;
};

void write_string_array(JsonBuffer &out,
                        const std::vector<std::string> &values) {
  out << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << quoted(values[i]);
  }
  out << "]";
}

template <typename Integer>
  requires(std::is_integral_v<Integer>)
void write_integer_array(JsonBuffer &out, const std::vector<Integer> &values) {
  out << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << values[i];
  }
  out << "]";
}

[[nodiscard]] bool has_clinical_annotation(const SnpCall &snp) {
  return !snp.clinical_assertions.empty();
}

void write_clinical_annotation(JsonBuffer &out, const SnpCall &snp) {
  out << "{"
      << "\"schema_version\":" << quoted(kClinicalAnnotationSchemaVersion)
      << ","
      << "\"conflict_status\":" << quoted(snp.clinical_conflict_status) << ","
      << "\"consensus_significance\":"
      << quoted(snp.clinical_consensus_significance) << ","
      << "\"pathogenicity\":" << quoted(snp.pathogenicity) << ","
      << "\"phenotype\":" << quoted(snp.phenotype) << ","
      << "\"references\":";
  write_string_array(out, snp.references);
  out << ",\"sources\":";
  write_string_array(out, snp.sources);
  out << ",\"source\":" << quoted("local-cache") << ","
      << "\"assertions\":[";
  for (std::size_t i = 0; i < snp.clinical_assertions.size(); ++i) {
    if (i != 0U) {
      out << ",";
    }
    const auto &assertion = snp.clinical_assertions[i];
    out << "{"
        << "\"source\":" << quoted(assertion.source) << ","
        << "\"assertion_id\":" << quoted(assertion.assertion_id) << ","
        << "\"allele_id\":" << quoted(assertion.allele_id) << ","
        << "\"disease\":" << quoted(assertion.disease) << ","
        << "\"clinical_significance\":"
        << quoted(assertion.clinical_significance) << ","
        << "\"normalized_significance\":"
        << quoted(assertion.normalized_significance) << ","
        << "\"review_status\":" << quoted(assertion.review_status) << ","
        << "\"assertion_date\":" << quoted(assertion.assertion_date) << ","
        << "\"source_url\":" << quoted(assertion.source_url) << ","
        << "\"references\":";
    write_string_array(out, assertion.references);
    out << ",\"resource_version\":" << quoted(assertion.resource_version) << ","
        << "\"retrieved_at\":" << quoted(assertion.retrieved_at) << "}";
  }
  out << "]";
  if (!snp.clinvar_allele_id.empty()) {
    out << ",\"clinvar_allele_id\":" << quoted(snp.clinvar_allele_id);
  }
  if (!snp.mitomap_url.empty()) {
    out << ",\"mitomap_url\":" << quoted(snp.mitomap_url);
  }
  out << "}";
}

void write_quality_summary(JsonBuffer &out, std::size_t count,
                           double quality_sum, std::uint8_t minimum,
                           std::uint8_t maximum) {
  out << "{\"count\":" << count << ",\"mean_phred\":";
  if (count == 0U) {
    out << "null,\"min_phred\":null,\"max_phred\":null}";
    return;
  }
  out << (quality_sum / static_cast<double>(count)) << ","
      << "\"min_phred\":" << static_cast<unsigned int>(minimum) << ","
      << "\"max_phred\":" << static_cast<unsigned int>(maximum) << "}";
}

void write_nullable_mean(JsonBuffer &out, double sum, std::size_t count) {
  if (count == 0U) {
    out << "null";
  } else {
    out << (sum / static_cast<double>(count));
  }
}

void write_string_map(JsonBuffer &out,
                      const std::map<std::string, std::string> &values) {
  out << "{";
  std::size_t index = 0;
  for (const auto &[key, value] : values) {
    if (index++ != 0) {
      out << ",";
    }
    out << quoted(key) << ":" << quoted(value);
  }
  out << "}";
}

[[nodiscard]] std::optional<SvCall> parse_sv_tag(std::string_view read_id,
                                                 std::string_view token) {
  const auto colon = token.find(':');
  const auto dash =
      token.find('-', colon == std::string_view::npos ? 0 : colon + 1);
  if (colon == std::string_view::npos || dash == std::string_view::npos) {
    return std::nullopt;
  }

  const auto type = token.substr(0, colon);
  const auto start = parse_size(token.substr(colon + 1, dash - colon - 1));
  const auto end = parse_size(token.substr(dash + 1));
  if (!start || !end || *start == 0 || *end < *start) {
    return std::nullopt;
  }

  SvCall sv;
  sv.type = std::string(type);
  std::transform(
      sv.type.begin(), sv.type.end(), sv.type.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (sv.type == "del") {
    sv.type = "deletion";
  } else if (sv.type == "ins") {
    sv.type = "insertion";
  }
  sv.start = *start;
  sv.end = *end;
  sv.length = sv.end - sv.start + 1;
  sv.supporting_reads.emplace_back(read_id);
  sv.evidence_sources.emplace_back("development_tag");
  sv.known_event = sv.type == "deletion" && sv.start >= 8400 &&
                   sv.start <= 8500 && sv.end >= 13400 && sv.end <= 13550;
  sv.id =
      sv.type + ":" + std::to_string(sv.start) + "-" + std::to_string(sv.end);
  return sv;
}

[[nodiscard]] std::vector<SvCall> sv_tags_from_header(const ReadRecord &read) {
  std::vector<SvCall> calls;
  std::istringstream iss(read.id);
  std::string token;
  while (iss >> token) {
    constexpr std::string_view prefix = "sv=";
    if (token.rfind(prefix, 0) == 0) {
      auto call =
          parse_sv_tag(read.id, std::string_view(token).substr(prefix.size()));
      if (call) {
        calls.push_back(std::move(*call));
      }
    }
  }
  return calls;
}

[[nodiscard]] std::vector<SvCall>
parse_cigar_svs(const ReadRecord &read, std::size_t sv_min_length,
                std::size_t reference_length) {
  std::vector<SvCall> calls;
  if (read.cigar_operations.empty() || (read.flags & 0x4U) != 0U ||
      read.reference_start == 0 || read.reference_name.empty() ||
      read.reference_name == "*" ||
      looks_like_nuclear_contig(read.reference_name)) {
    return calls;
  }

  std::size_t reference_cursor =
      read.reference_start == 0 ? 1 : read.reference_start;
  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;
    if ((op == 'D' || op == 'N') && len >= sv_min_length) {
      if (len - 1U >
          std::numeric_limits<std::size_t>::max() - reference_cursor) {
        throw AnalysisError(
            AnalysisErrorCode::input_parse_failed,
            "CIGAR structural-variant coordinate overflow for read '" +
                read.id + "'");
      }
      SvCall sv;
      sv.type = "deletion";
      sv.start = reference_cursor;
      sv.end = reference_cursor + len - 1;
      sv.length = len;
      sv.supporting_reads.emplace_back(read.id);
      sv.evidence_sources.emplace_back("cigar");
      sv.known_event = sv.start >= 8400 && sv.start <= 8500 &&
                       sv.end >= 13400 && sv.end <= 13550;
      sv.id =
          "deletion:" + std::to_string(sv.start) + "-" + std::to_string(sv.end);
      calls.push_back(std::move(sv));
    } else if (op == 'I' && len >= sv_min_length) {
      SvCall sv;
      sv.type = "insertion";
      sv.start =
          reference_cursor > 1U ? reference_cursor - 1U : reference_length;
      sv.end = sv.start;
      sv.length = len;
      sv.supporting_reads.emplace_back(read.id);
      sv.evidence_sources.emplace_back("cigar");
      sv.id =
          "insertion:" + std::to_string(sv.start) + "+" + std::to_string(len);
      calls.push_back(std::move(sv));
    } else if (op == 'S' && len >= sv_min_length) {
      SvCall sv;
      const bool leading =
          reference_cursor ==
          (read.reference_start == 0 ? 1 : read.reference_start);
      sv.type = leading ? "soft_clip_left" : "soft_clip_right";
      sv.start = leading ? reference_cursor
                         : (reference_cursor > 1U ? reference_cursor - 1U
                                                  : reference_length);
      sv.end = reference_cursor;
      sv.length = len;
      sv.supporting_reads.emplace_back(read.id);
      sv.evidence_sources.emplace_back("cigar");
      sv.id =
          sv.type + ":" + std::to_string(sv.start) + "+" + std::to_string(len);
      calls.push_back(std::move(sv));
    }

    if (op == 'M' || op == '=' || op == 'X' || op == 'D' || op == 'N') {
      if (len > std::numeric_limits<std::size_t>::max() - reference_cursor) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "CIGAR reference coordinate overflow for read '" +
                                read.id + "'");
      }
      reference_cursor += len;
    }
  }

  return calls;
}

struct AlignmentSegment {
  std::string reference_name;
  std::size_t reference_start = 0;
  std::size_t reference_end = 0;
  std::size_t query_start = 0;
  std::size_t query_end = 0;
  char strand = '+';
  std::uint8_t mapping_quality = 0;
};

[[nodiscard]] bool is_mitochondrial_contig(std::string_view name) {
  std::string normalized(name);
  std::transform(
      normalized.begin(), normalized.end(), normalized.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return normalized == "mt" || normalized == "chrm" || normalized == "m" ||
         normalized == "nc_012920.1" || normalized == "nc_012920";
}

[[nodiscard]] std::optional<AlignmentSegment>
make_alignment_segment(std::string reference_name, std::size_t reference_start,
                       char strand, std::uint8_t mapping_quality,
                       const std::vector<CigarOperation> &operations) {
  if (!is_mitochondrial_contig(reference_name) || reference_start == 0 ||
      operations.empty()) {
    return std::nullopt;
  }
  std::size_t leading_clip = 0;
  std::size_t query_extent = 0;
  std::size_t query_aligned = 0;
  std::size_t reference_span = 0;
  bool seen_alignment = false;
  for (const auto &operation : operations) {
    if (operation.code == 'M' || operation.code == 'I' ||
        operation.code == 'S' || operation.code == 'H' ||
        operation.code == '=' || operation.code == 'X') {
      if (operation.length >
          std::numeric_limits<std::size_t>::max() - query_extent) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment query extent overflow");
      }
      query_extent += operation.length;
    }
    if (!seen_alignment && (operation.code == 'S' || operation.code == 'H')) {
      if (operation.length >
          std::numeric_limits<std::size_t>::max() - leading_clip) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment query coordinate overflow");
      }
      leading_clip += operation.length;
      continue;
    }
    if (operation.code == 'M' || operation.code == '=' ||
        operation.code == 'X') {
      seen_alignment = true;
      if (operation.length >
              std::numeric_limits<std::size_t>::max() - query_aligned ||
          operation.length >
              std::numeric_limits<std::size_t>::max() - reference_span) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment segment length overflow");
      }
      query_aligned += operation.length;
      reference_span += operation.length;
    } else if (operation.code == 'I') {
      seen_alignment = true;
      if (operation.length >
          std::numeric_limits<std::size_t>::max() - query_aligned) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment query length overflow");
      }
      query_aligned += operation.length;
    } else if (operation.code == 'D' || operation.code == 'N') {
      seen_alignment = true;
      if (operation.length >
          std::numeric_limits<std::size_t>::max() - reference_span) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment reference length overflow");
      }
      reference_span += operation.length;
    }
  }
  if (query_aligned == 0 || reference_span == 0) {
    return std::nullopt;
  }
  if (leading_clip > query_extent ||
      query_aligned > query_extent - leading_clip) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "alignment clipping exceeds query extent");
  }
  const auto query_start = strand == '+'
                               ? leading_clip
                               : query_extent - leading_clip - query_aligned;
  if (reference_span - 1U >
          std::numeric_limits<std::size_t>::max() - reference_start ||
      query_aligned - 1U >
          std::numeric_limits<std::size_t>::max() - query_start) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "alignment segment coordinate overflow");
  }
  return AlignmentSegment{std::move(reference_name),
                          reference_start,
                          reference_start + reference_span - 1U,
                          query_start,
                          query_start + query_aligned - 1U,
                          strand,
                          mapping_quality};
}

[[nodiscard]] std::vector<AlignmentSegment>
alignment_segments(const ReadRecord &read) {
  std::vector<AlignmentSegment> segments;
  if (const auto primary =
          make_alignment_segment(read.reference_name, read.reference_start,
                                 (read.flags & 0x10U) != 0U ? '-' : '+',
                                 read.mapping_quality, read.cigar_operations)) {
    segments.push_back(*primary);
  }
  for (const auto &alignment : read.supplementary_alignments) {
    if (const auto segment = make_alignment_segment(
            alignment.reference_name, alignment.reference_start,
            alignment.strand, alignment.mapping_quality,
            alignment.cigar_operations)) {
      segments.push_back(*segment);
    }
  }
  std::sort(
      segments.begin(), segments.end(), [](const auto &lhs, const auto &rhs) {
        return std::tie(lhs.query_start, lhs.query_end, lhs.reference_start,
                        lhs.strand) < std::tie(rhs.query_start, rhs.query_end,
                                               rhs.reference_start, rhs.strand);
      });
  segments.erase(std::unique(segments.begin(), segments.end(),
                             [](const auto &lhs, const auto &rhs) {
                               return lhs.query_start == rhs.query_start &&
                                      lhs.query_end == rhs.query_end &&
                                      lhs.reference_start ==
                                          rhs.reference_start &&
                                      lhs.reference_end == rhs.reference_end &&
                                      lhs.strand == rhs.strand;
                             }),
                 segments.end());
  return segments;
}

[[nodiscard]] std::vector<SvCall>
split_alignment_svs(const ReadRecord &read, std::size_t reference_length,
                    std::size_t sv_min_length) {
  std::vector<SvCall> calls;
  const auto segments = alignment_segments(read);
  if (segments.size() < 2 || reference_length == 0) {
    return calls;
  }
  for (std::size_t i = 1; i < segments.size(); ++i) {
    const auto &previous = segments[i - 1];
    const auto &current = segments[i];
    const auto query_gap = static_cast<std::int64_t>(current.query_start) -
                           static_cast<std::int64_t>(previous.query_end) - 1;
    SvCall event;
    event.supporting_reads.push_back(read.id);
    event.evidence_sources.emplace_back("split_alignment");
    event.segment_count = segments.size();
    event.orientations.push_back(std::string(1, previous.strand) + "/" +
                                 current.strand);

    const auto leaving_position = previous.strand == '+'
                                      ? previous.reference_end
                                      : previous.reference_start;
    const auto entering_position =
        current.strand == '+' ? current.reference_start : current.reference_end;
    if (leaving_position > reference_length ||
        entering_position > reference_length) {
      throw AnalysisError(
          AnalysisErrorCode::input_parse_failed,
          "split alignment coordinate exceeds the mitochondrial reference for "
          "read '" +
              read.id + "'");
    }

    if (previous.strand != current.strand) {
      event.type = "inversion";
      event.start = std::min(leaving_position, entering_position);
      event.end = std::max(leaving_position, entering_position);
      event.length = event.end - event.start + 1U;
    } else {
      const bool forward = previous.strand == '+';
      const bool follows_linear_reference =
          forward ? entering_position > leaving_position
                  : entering_position < leaving_position;
      const auto linear_gap =
          follows_linear_reference
              ? (forward ? entering_position - leaving_position - 1U
                         : leaving_position - entering_position - 1U)
              : 0U;
      const auto overlap =
          follows_linear_reference
              ? 0U
              : (forward ? leaving_position - entering_position + 1U
                         : entering_position - leaving_position + 1U);
      const auto circular_gap =
          follows_linear_reference
              ? reference_length
              : (forward ? reference_length - leaving_position +
                               entering_position - 1U
                         : leaving_position - 1U + reference_length -
                               entering_position);
      const bool crosses_origin =
          !follows_linear_reference && circular_gap < overlap;
      const bool ambiguous_topology =
          !follows_linear_reference && circular_gap == overlap;

      if (crosses_origin) {
        event.type = "circular_origin";
        event.start = std::max(leaving_position, entering_position);
        event.end = std::min(leaving_position, entering_position);
        event.length = circular_gap;
      } else if (ambiguous_topology && overlap >= sv_min_length) {
        event.type = "ambiguous_adjacency";
        event.start = std::min(leaving_position, entering_position);
        event.end = std::max(leaving_position, entering_position);
        event.length = overlap;
      } else if (follows_linear_reference && linear_gap >= sv_min_length) {
        event.type = "deletion";
        event.start = std::min(leaving_position, entering_position) + 1U;
        event.end = std::max(leaving_position, entering_position) - 1U;
        event.length = linear_gap;
      } else if (!follows_linear_reference && overlap >= sv_min_length) {
        event.type = "duplication";
        event.start = std::min(leaving_position, entering_position);
        event.end = std::max(leaving_position, entering_position);
        event.length = overlap;
      } else if (query_gap >= static_cast<std::int64_t>(sv_min_length)) {
        event.type = "insertion";
        event.start = std::min(leaving_position, entering_position);
        event.end = event.start;
        event.length = static_cast<std::size_t>(query_gap);
      } else {
        continue;
      }
    }
    event.known_event = event.type == "deletion" && event.start >= 8400 &&
                        event.start <= 8500 && event.end >= 13400 &&
                        event.end <= 13550;
    if (event.type == "insertion") {
      event.id = event.type + ":" + std::to_string(event.start) + "+" +
                 std::to_string(event.length);
    } else {
      event.id = event.type + ":" + std::to_string(event.start) + "-" +
                 std::to_string(event.end);
    }
    calls.push_back(std::move(event));
  }
  return calls;
}

struct CanonicalJunctionPath {
  std::vector<std::string> ids;
  std::vector<std::string> orientations;
};

[[nodiscard]] char flip_strand(char strand) {
  if (strand == '+') {
    return '-';
  }
  if (strand == '-') {
    return '+';
  }
  throw AnalysisError(AnalysisErrorCode::internal_error,
                      "invalid strand in split-junction orientation");
}

[[nodiscard]] std::string
reverse_complement_orientation(std::string_view orientation) {
  if (orientation.size() != 3U || orientation[1] != '/') {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "invalid split-junction orientation '" +
                            std::string(orientation) + "'");
  }
  return std::string{flip_strand(orientation[2]), '/',
                     flip_strand(orientation[0])};
}

[[nodiscard]] CanonicalJunctionPath
canonical_junction_path(const std::vector<SvCall> &junctions) {
  CanonicalJunctionPath forward;
  std::vector<std::string> forward_tokens;
  forward.ids.reserve(junctions.size());
  forward.orientations.reserve(junctions.size());
  forward_tokens.reserve(junctions.size());
  for (const auto &junction : junctions) {
    if (junction.orientations.size() != 1U) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "split junction does not have exactly one orientation");
    }
    forward.ids.push_back(junction.id);
    forward.orientations.push_back(junction.orientations.front());
    forward_tokens.push_back(junction.id + "@" + junction.orientations.front());
  }

  CanonicalJunctionPath reverse;
  std::vector<std::string> reverse_tokens;
  reverse.ids.reserve(junctions.size());
  reverse.orientations.reserve(junctions.size());
  reverse_tokens.reserve(junctions.size());
  for (auto it = junctions.rbegin(); it != junctions.rend(); ++it) {
    const auto orientation =
        reverse_complement_orientation(it->orientations.front());
    reverse.ids.push_back(it->id);
    reverse.orientations.push_back(orientation);
    reverse_tokens.push_back(it->id + "@" + orientation);
  }
  return reverse_tokens < forward_tokens ? std::move(reverse)
                                         : std::move(forward);
}

[[nodiscard]] std::optional<ComplexSvCall>
coalesce_complex_sv(std::string_view read_id,
                    const std::vector<SvCall> &junctions) {
  if (junctions.size() < 2U) {
    return std::nullopt;
  }

  ComplexSvCall event;
  auto canonical_path = canonical_junction_path(junctions);
  event.junction_ids = std::move(canonical_path.ids);
  event.junction_orientations = std::move(canonical_path.orientations);
  event.supporting_reads.emplace_back(read_id);
  for (const auto &junction : junctions) {
    event.segment_count = std::max(event.segment_count, junction.segment_count);
  }

  // The exact canonical path is the identifier. This is deliberately not a
  // truncated hash: an event ID remains collision-free and independently
  // reviewable, while reversal canonicalization merges opposite-strand reads.
  event.id = "complex:";
  for (std::size_t i = 0; i < event.junction_ids.size(); ++i) {
    if (i != 0U) {
      event.id.push_back('|');
    }
    event.id += event.junction_ids[i];
    event.id.push_back('@');
    event.id += event.junction_orientations[i];
  }
  return event;
}

void merge_sv(std::map<std::string, SvCall> &svs, SvCall sv) {
  const auto existing = svs.find(sv.id);
  if (existing == svs.end()) {
    for (auto *values :
         {&sv.supporting_reads, &sv.evidence_sources, &sv.orientations}) {
      std::sort(values->begin(), values->end());
      values->erase(std::unique(values->begin(), values->end()), values->end());
    }
    svs.emplace(sv.id, std::move(sv));
    return;
  }

  auto &merged = existing->second;
  if (merged.type != sv.type || merged.start != sv.start ||
      merged.end != sv.end || merged.length != sv.length) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "canonical SV ID collision for '" + sv.id + "'");
  }
  const auto merge_unique = [](std::vector<std::string> &target,
                               const std::vector<std::string> &incoming) {
    target.insert(target.end(), incoming.begin(), incoming.end());
    std::sort(target.begin(), target.end());
    target.erase(std::unique(target.begin(), target.end()), target.end());
  };
  merge_unique(merged.supporting_reads, sv.supporting_reads);
  merge_unique(merged.evidence_sources, sv.evidence_sources);
  merge_unique(merged.orientations, sv.orientations);
  merged.known_event = merged.known_event || sv.known_event;
  merged.segment_count = std::max(merged.segment_count, sv.segment_count);
}

void merge_complex_sv(std::map<std::string, ComplexSvCall> &events,
                      ComplexSvCall event) {
  const auto existing = events.find(event.id);
  if (existing == events.end()) {
    std::sort(event.supporting_reads.begin(), event.supporting_reads.end());
    event.supporting_reads.erase(std::unique(event.supporting_reads.begin(),
                                             event.supporting_reads.end()),
                                 event.supporting_reads.end());
    events.emplace(event.id, std::move(event));
    return;
  }

  auto &merged = existing->second;
  if (merged.junction_ids != event.junction_ids ||
      merged.junction_orientations != event.junction_orientations) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "canonical complex-SV ID collision for '" + event.id +
                            "'");
  }
  merged.supporting_reads.insert(merged.supporting_reads.end(),
                                 event.supporting_reads.begin(),
                                 event.supporting_reads.end());
  std::sort(merged.supporting_reads.begin(), merged.supporting_reads.end());
  merged.supporting_reads.erase(std::unique(merged.supporting_reads.begin(),
                                            merged.supporting_reads.end()),
                                merged.supporting_reads.end());
  merged.segment_count = std::max(merged.segment_count, event.segment_count);
}

[[nodiscard]] std::size_t effective_thread_count(std::size_t requested_threads,
                                                 std::size_t work_items) {
  if (work_items <= 1) {
    return 1;
  }
  const auto hardware_threads =
      std::max(1U, std::thread::hardware_concurrency());
  const std::size_t requested = requested_threads == 0 ? 1 : requested_threads;
  return std::max<std::size_t>(
      1, std::min<std::size_t>({requested, work_items,
                                static_cast<std::size_t>(hardware_threads)}));
}

[[nodiscard]] FeatureExtractionResult extract_read_feature(
    const ReadRecord &read, const std::string &reference,
    const AnalysisConfig &config,
    const std::map<std::string, SnpCall> &clinical_annotations) {
  throw_if_cancelled(config);
  FeatureExtractionResult result;
  auto &feature = result.feature;
  feature.id = read.id;
  feature.length = read.sequence.size();
  feature.mean_quality = mean_quality(read.qualities);
  feature.numt_evidence = numt_evidence(read, reference.size());
  feature.numt_score = numt_score(feature.numt_evidence);
  feature.filtered_numt =
      config.filter_numt && feature.numt_score > config.numt_threshold;
  feature.mapping_quality = read.mapping_quality;
  feature.flags = read.flags;
  feature.reference_name = read.reference_name;
  feature.aux_tags = read.aux_tags;

  if (!feature.filtered_numt) {
    feature.snps = call_snps_from_alignment(read, reference, config);
    if (feature.snps.empty() && read.cigar.empty()) {
      feature.snps = call_snps_from_reference_span(read, reference, config);
    }
    if (config.allow_development_tags) {
      merge_snps(feature.snps, snp_tags_from_header(read));
    }
    if (std::any_of(
            feature.snps.begin(), feature.snps.end(), [&](const auto &snp) {
              return snp.position == 0 || snp.position > reference.size();
            })) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "SNP coordinate is outside the reference for read '" +
                              read.id + "'");
    }
    apply_clinical_annotations(feature.snps, clinical_annotations);
    feature.haplogroup_markers.reserve(feature.snps.size());
    for (const auto &snp : feature.snps) {
      ObservedPhyloMutation mutation;
      mutation.position = snp.position;
      mutation.kind = PhyloMutationKind::substitution;
      mutation.alternate = snp.alternate;
      mutation.encoded =
          std::to_string(snp.position) + std::string(1, snp.alternate);
      feature.haplogroup_markers.push_back(std::move(mutation));
    }
    auto aligned_phylo_indels =
        call_phylo_indels_from_alignment(read, config, reference);
    std::vector<std::size_t> aligned_phylo_indel_positions;
    aligned_phylo_indel_positions.reserve(aligned_phylo_indels.size());
    for (const auto &mutation : aligned_phylo_indels) {
      aligned_phylo_indel_positions.push_back(mutation.position);
    }
    merge_haplogroup_markers(feature.haplogroup_markers,
                             std::move(aligned_phylo_indels));
    if (config.allow_development_tags) {
      merge_haplogroup_markers(feature.haplogroup_markers,
                               phylo_tags_from_header(read));
    }
    const auto callable_ranges =
        callable_phylo_ranges_from_alignment(read, reference, config);
    feature.haplogroup_range_known = callable_ranges.known;
    feature.haplogroup_ranges = callable_ranges.ranges;
    if (feature.haplogroup_range_known) {
      for (const auto position : aligned_phylo_indel_positions) {
        feature.haplogroup_ranges.emplace_back(position, position);
      }
      feature.haplogroup_ranges =
          normalize_ranges(std::move(feature.haplogroup_ranges));
    }
    if (config.allow_development_tags) {
      if (auto development_ranges = phylo_ranges_from_header(read)) {
        feature.haplogroup_range_known = true;
        feature.haplogroup_ranges = std::move(*development_ranges);
      }
    }
    if (config.allow_development_tags) {
      result.svs = sv_tags_from_header(read);
    }
    auto cigar_svs =
        parse_cigar_svs(read, config.sv_min_length, reference.size());
    auto split_svs =
        split_alignment_svs(read, reference.size(), config.sv_min_length);
    if (auto complex_event = coalesce_complex_sv(read.id, split_svs)) {
      feature.complex_event_ids.push_back(complex_event->id);
      result.complex_events.push_back(std::move(*complex_event));
    }
    if (!split_svs.empty()) {
      cigar_svs.erase(std::remove_if(cigar_svs.begin(), cigar_svs.end(),
                                     [](const auto &sv) {
                                       return sv.type == "soft_clip_left" ||
                                              sv.type == "soft_clip_right";
                                     }),
                      cigar_svs.end());
    }
    result.svs.insert(result.svs.end(),
                      std::make_move_iterator(cigar_svs.begin()),
                      std::make_move_iterator(cigar_svs.end()));
    result.svs.insert(result.svs.end(),
                      std::make_move_iterator(split_svs.begin()),
                      std::make_move_iterator(split_svs.end()));
    if (std::any_of(result.svs.begin(), result.svs.end(), [&](const auto &sv) {
          return sv.start == 0 || sv.end > reference.size();
        })) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "SV coordinate is outside the reference for read '" +
                              read.id + "'");
    }

    feature.sv_ids.reserve(result.svs.size());
    for (const auto &sv : result.svs) {
      feature.sv_ids.push_back(sv.id);
    }
    std::sort(feature.sv_ids.begin(), feature.sv_ids.end());
    feature.sv_ids.erase(
        std::unique(feature.sv_ids.begin(), feature.sv_ids.end()),
        feature.sv_ids.end());
    std::sort(feature.complex_event_ids.begin(),
              feature.complex_event_ids.end());
    feature.complex_event_ids.erase(
        std::unique(feature.complex_event_ids.begin(),
                    feature.complex_event_ids.end()),
        feature.complex_event_ids.end());
  }

  return result;
}

[[nodiscard]] std::vector<FeatureExtractionResult> extract_features_parallel(
    const std::vector<ReadRecord> &reads, const std::string &reference,
    const AnalysisConfig &config,
    const std::map<std::string, SnpCall> &clinical_annotations) {
  std::vector<FeatureExtractionResult> results(reads.size());
  if (reads.empty()) {
    return results;
  }

  const std::size_t worker_count =
      effective_thread_count(config.threads, reads.size());
  if (worker_count == 1) {
    for (std::size_t i = 0; i < reads.size(); ++i) {
      throw_if_cancelled(config);
      results[i] = extract_read_feature(reads[i], reference, config,
                                        clinical_annotations);
    }
    return results;
  }

  constexpr std::size_t kChunkSize = 32;
  std::atomic<std::size_t> next_index{0};
  std::atomic<bool> failed{false};
  std::vector<std::exception_ptr> errors(worker_count);
  std::vector<std::jthread> workers;
  workers.reserve(worker_count);

  for (std::size_t worker_id = 0; worker_id < worker_count; ++worker_id) {
    workers.emplace_back([&, worker_id] {
      try {
        while (!failed.load(std::memory_order_relaxed)) {
          throw_if_cancelled(config);
          const std::size_t begin =
              next_index.fetch_add(kChunkSize, std::memory_order_relaxed);
          if (begin >= reads.size()) {
            break;
          }
          const std::size_t end = std::min(reads.size(), begin + kChunkSize);
          for (std::size_t i = begin; i < end; ++i) {
            throw_if_cancelled(config);
            results[i] = extract_read_feature(reads[i], reference, config,
                                              clinical_annotations);
          }
        }
      } catch (...) {
        errors[worker_id] = std::current_exception();
        failed.store(true, std::memory_order_relaxed);
      }
    });
  }

  for (auto &worker : workers) {
    worker.join();
  }
  for (const auto &error : errors) {
    if (error) {
      std::rethrow_exception(error);
    }
  }

  return results;
}

using FeatureTokens = std::vector<std::string>;

[[nodiscard]] FeatureTokens feature_tokens(const ReadFeature &feature) {
  FeatureTokens tokens;
  tokens.reserve(feature.snps.size() + feature.haplogroup_markers.size() +
                 feature.sv_ids.size() + feature.complex_event_ids.size() + 1U);
  for (const auto &snp : feature.snps) {
    tokens.push_back("snp:" +
                     snp_key(snp.position, snp.reference, snp.alternate));
  }
  for (const auto &marker : feature.haplogroup_markers) {
    if (marker.kind != PhyloMutationKind::substitution) {
      tokens.push_back("phylo_indel:" + marker.encoded);
    }
  }
  for (const auto &sv_id : feature.sv_ids) {
    tokens.push_back("sv:" + sv_id);
  }
  for (const auto &complex_event_id : feature.complex_event_ids) {
    tokens.push_back("complex_sv:" + complex_event_id);
  }
  if (tokens.empty()) {
    tokens.push_back("length:" + std::to_string(feature.length / 100U));
  }
  std::sort(tokens.begin(), tokens.end());
  tokens.erase(std::unique(tokens.begin(), tokens.end()), tokens.end());
  return tokens;
}

[[nodiscard]] double token_distance(const FeatureTokens &lhs,
                                    const FeatureTokens &rhs) {
  std::size_t intersection = 0;
  auto lhs_it = lhs.begin();
  auto rhs_it = rhs.begin();
  while (lhs_it != lhs.end() && rhs_it != rhs.end()) {
    if (*lhs_it == *rhs_it) {
      ++intersection;
      ++lhs_it;
      ++rhs_it;
    } else if (*lhs_it < *rhs_it) {
      ++lhs_it;
    } else {
      ++rhs_it;
    }
  }

  const std::size_t union_size = lhs.size() + rhs.size() - intersection;
  if (union_size == 0) {
    return 0.0;
  }
  return 1.0 -
         (static_cast<double>(intersection) / static_cast<double>(union_size));
}

struct TokenProfile {
  FeatureTokens tokens;
  std::vector<std::size_t> feature_indices;
};

struct Neighborhood {
  std::vector<std::size_t> profiles;
  std::size_t read_count = 0;
};

using InvertedTokenIndex =
    std::unordered_map<std::string_view, std::vector<std::size_t>>;

[[nodiscard]] InvertedTokenIndex
build_inverted_token_index(const std::vector<TokenProfile> &profiles) {
  InvertedTokenIndex index;
  for (std::size_t profile_index = 0; profile_index < profiles.size();
       ++profile_index) {
    for (const auto &token : profiles[profile_index].tokens) {
      index[token].push_back(profile_index);
    }
  }
  return index;
}

[[nodiscard]] Neighborhood
region_query(const std::vector<TokenProfile> &profiles,
             const InvertedTokenIndex &inverted_index,
             std::vector<std::size_t> &candidate_marks, std::size_t &mark_epoch,
             std::size_t index, double epsilon) {
  Neighborhood neighborhood;
  std::vector<std::size_t> candidates;
  if (epsilon >= 1.0) {
    candidates.resize(profiles.size());
    std::iota(candidates.begin(), candidates.end(), 0U);
  } else {
    if (mark_epoch == std::numeric_limits<std::size_t>::max()) {
      std::fill(candidate_marks.begin(), candidate_marks.end(), 0U);
      mark_epoch = 0;
    }
    ++mark_epoch;
    for (const auto &token : profiles[index].tokens) {
      const auto posting = inverted_index.find(token);
      if (posting == inverted_index.end()) {
        continue;
      }
      for (const auto candidate : posting->second) {
        if (candidate_marks[candidate] != mark_epoch) {
          candidate_marks[candidate] = mark_epoch;
          candidates.push_back(candidate);
        }
      }
    }
    std::sort(candidates.begin(), candidates.end());
  }

  for (const auto candidate : candidates) {
    const auto lhs_size = profiles[index].tokens.size();
    const auto rhs_size = profiles[candidate].tokens.size();
    const double maximum_similarity =
        static_cast<double>(std::min(lhs_size, rhs_size)) /
        static_cast<double>(std::max(lhs_size, rhs_size));
    if (maximum_similarity < 1.0 - epsilon) {
      continue;
    }
    if (token_distance(profiles[index].tokens, profiles[candidate].tokens) <=
        epsilon) {
      neighborhood.profiles.push_back(candidate);
      neighborhood.read_count += profiles[candidate].feature_indices.size();
    }
  }
  return neighborhood;
}

void append_new_neighbors(std::vector<std::size_t> &queue,
                          std::vector<unsigned char> &queued,
                          const std::vector<std::size_t> &candidates) {
  for (const std::size_t candidate : candidates) {
    if (candidate >= queued.size() || queued[candidate] != 0U) {
      continue;
    }
    queued[candidate] = 1U;
    queue.push_back(candidate);
  }
}

void assign_dbscan_clusters(std::vector<ReadFeature> &features, double epsilon,
                            std::size_t min_cluster_size,
                            const AnalysisConfig &config) {
  std::vector<std::size_t> active_indices;
  active_indices.reserve(features.size());
  for (std::size_t i = 0; i < features.size(); ++i) {
    if (!features[i].filtered_numt) {
      active_indices.push_back(i);
      features[i].cluster_id = -1;
      features[i].outlier = false;
    }
  }

  std::vector<TokenProfile> profiles;
  profiles.reserve(active_indices.size());
  {
    std::map<FeatureTokens, std::size_t> profile_lookup;
    for (const auto index : active_indices) {
      auto tokens = feature_tokens(features[index]);
      const auto existing = profile_lookup.find(tokens);
      if (existing != profile_lookup.end()) {
        profiles[existing->second].feature_indices.push_back(index);
        continue;
      }
      const auto profile_index = profiles.size();
      profile_lookup.emplace(tokens, profile_index);
      profiles.push_back({std::move(tokens), {index}});
    }
  }

  constexpr int kUnvisited = -2;
  std::vector<int> labels(profiles.size(), kUnvisited);
  const auto inverted_index = build_inverted_token_index(profiles);
  std::vector<std::size_t> candidate_marks(profiles.size(), 0U);
  std::size_t mark_epoch = 0;
  const std::size_t min_points = std::max<std::size_t>(1, min_cluster_size);
  int next_cluster_id = 0;

  for (std::size_t point = 0; point < profiles.size(); ++point) {
    throw_if_cancelled(config);
    if (labels[point] != kUnvisited) {
      continue;
    }

    auto neighborhood = region_query(profiles, inverted_index, candidate_marks,
                                     mark_epoch, point, epsilon);
    if (neighborhood.read_count < min_points) {
      labels[point] = -1;
      continue;
    }

    const int cluster_id = next_cluster_id++;
    labels[point] = cluster_id;
    auto &neighbors = neighborhood.profiles;
    std::vector<unsigned char> queued(profiles.size(), 0U);
    for (const auto neighbor : neighbors) {
      if (neighbor < queued.size()) {
        queued[neighbor] = 1U;
      }
    }
    for (std::size_t cursor = 0; cursor < neighbors.size(); ++cursor) {
      throw_if_cancelled(config);
      const std::size_t neighbor = neighbors[cursor];
      if (labels[neighbor] == -1) {
        labels[neighbor] = cluster_id;
      }
      if (labels[neighbor] != kUnvisited) {
        continue;
      }
      labels[neighbor] = cluster_id;
      auto expanded = region_query(profiles, inverted_index, candidate_marks,
                                   mark_epoch, neighbor, epsilon);
      if (expanded.read_count >= min_points) {
        append_new_neighbors(neighbors, queued, expanded.profiles);
      }
    }
  }

  for (std::size_t i = 0; i < profiles.size(); ++i) {
    for (const auto feature_index : profiles[i].feature_indices) {
      auto &feature = features[feature_index];
      feature.cluster_id = labels[i];
      feature.outlier = labels[i] < 0;
    }
  }
}

#ifdef MITO_HAS_HDBSCAN_CPP
[[nodiscard]] bool assign_hdbscan_clusters(std::vector<ReadFeature> &features,
                                           std::size_t min_cluster_size,
                                           const AnalysisConfig &config) {
  std::vector<std::size_t> active_indices;
  active_indices.reserve(features.size());
  for (std::size_t i = 0; i < features.size(); ++i) {
    if (!features[i].filtered_numt) {
      active_indices.push_back(i);
    }
  }
  if (active_indices.empty()) {
    return true;
  }

  std::vector<FeatureTokens> token_sets;
  token_sets.reserve(active_indices.size());
  std::map<std::string, std::size_t> token_columns;
  for (const auto index : active_indices) {
    token_sets.push_back(feature_tokens(features[index]));
    for (const auto &token : token_sets.back()) {
      if (!token_columns.contains(token)) {
        token_columns.emplace(token, token_columns.size());
      }
    }
  }

  std::vector<std::vector<double>> dataset(
      token_sets.size(),
      std::vector<double>(std::max<std::size_t>(1, token_columns.size()), 0.0));
  for (std::size_t row = 0; row < token_sets.size(); ++row) {
    for (const auto &token : token_sets[row]) {
      dataset[row][token_columns.at(token)] = 1.0;
    }
  }

  try {
    throw_if_cancelled(config);
    Hdbscan hdbscan("");
    hdbscan.dataset = std::move(dataset);
    const auto min_points =
        static_cast<int>(std::max<std::size_t>(1, min_cluster_size));
    hdbscan.execute(min_points, min_points, "Euclidean");
    if (hdbscan.normalizedLabels_.size() != active_indices.size()) {
      return false;
    }
    for (std::size_t i = 0; i < active_indices.size(); ++i) {
      const int label = hdbscan.normalizedLabels_[i];
      auto &feature = features[active_indices[i]];
      feature.cluster_id = label <= 0 ? -1 : label - 1;
      feature.outlier = label <= 0;
    }
    return true;
  } catch (const std::bad_alloc &) {
    throw;
  } catch (...) {
    throw_if_cancelled(config);
    return false;
  }
}
#endif

void add_circular_coverage_span(std::vector<std::int64_t> &difference,
                                std::uint64_t &uniform_depth, std::size_t start,
                                std::size_t length) {
  const std::size_t reference_length = difference.size() - 1U;
  if (reference_length == 0 || length == 0) {
    return;
  }

  const auto complete_cycles =
      static_cast<std::uint64_t>(length / reference_length);
  if (complete_cycles >
      std::numeric_limits<std::uint64_t>::max() - uniform_depth) {
    throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                        "coverage depth overflow");
  }
  uniform_depth += complete_cycles;
  const std::size_t remainder = length % reference_length;
  if (remainder == 0) {
    return;
  }

  const std::size_t zero_based_start = (start - 1U) % reference_length;
  const std::size_t linear_end = zero_based_start + remainder;
  ++difference[zero_based_start];
  if (linear_end <= reference_length) {
    --difference[linear_end];
    return;
  }

  --difference[reference_length];
  ++difference[0];
  --difference[linear_end - reference_length];
}

[[nodiscard]] CoverageResult
compute_coverage(const std::vector<ReadRecord> &reads,
                 const std::vector<ReadFeature> &features,
                 const AnalysisConfig &config, std::size_t reference_length,
                 std::size_t bin_count = 180) {
  CoverageResult result;
  if (reference_length == 0 || bin_count == 0) {
    return result;
  }

  std::vector<std::int64_t> difference(reference_length + 1U, 0);
  std::uint64_t uniform_depth = 0;
  for (std::size_t read_index = 0; read_index < reads.size(); ++read_index) {
    throw_if_cancelled(config);
    const auto &read = reads[read_index];
    if (read_index >= features.size() || features[read_index].filtered_numt ||
        (read.flags & 0x4U) != 0U || read.reference_start == 0 ||
        read.reference_name.empty() || read.reference_name == "*" ||
        read.cigar_operations.empty() ||
        looks_like_nuclear_contig(read.reference_name)) {
      continue;
    }

    std::size_t reference_cursor =
        ((read.reference_start - 1U) % reference_length) + 1U;
    for (const auto &operation : read.cigar_operations) {
      if (operation.code == 'M' || operation.code == '=' ||
          operation.code == 'X') {
        add_circular_coverage_span(difference, uniform_depth, reference_cursor,
                                   operation.length);
      }
      if (operation.code == 'M' || operation.code == '=' ||
          operation.code == 'X' || operation.code == 'D' ||
          operation.code == 'N') {
        reference_cursor =
            ((reference_cursor - 1U + (operation.length % reference_length)) %
             reference_length) +
            1U;
      }
    }
  }

  std::vector<std::uint64_t> site_depth(reference_length, 0);
  std::int64_t running_delta = 0;
  long double total_depth = 0.0L;
  std::size_t sites_gt20 = 0;
  for (std::size_t i = 0; i < reference_length; ++i) {
    running_delta += difference[i];
    if (running_delta < 0) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "internal coverage accumulator underflow");
    }
    const auto variable_depth = static_cast<std::uint64_t>(running_delta);
    if (variable_depth >
        std::numeric_limits<std::uint64_t>::max() - uniform_depth) {
      throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                          "coverage depth overflow");
    }
    const auto depth = uniform_depth + variable_depth;
    if (depth > std::numeric_limits<std::size_t>::max()) {
      throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                          "coverage depth exceeds platform size limit");
    }
    site_depth[i] = depth;
    total_depth += static_cast<long double>(depth);
    result.max_depth =
        std::max(result.max_depth, static_cast<std::size_t>(depth));
    if (depth > 20U) {
      ++sites_gt20;
    }
  }

  result.mean_depth = static_cast<double>(
      total_depth / static_cast<long double>(reference_length));
  result.pct_sites_gt20x = (static_cast<double>(sites_gt20) * 100.0) /
                           static_cast<double>(reference_length);

  auto &bins = result.bins;
  bins.reserve(bin_count);
  const std::size_t bin_width =
      std::max<std::size_t>(1, reference_length / bin_count);
  for (std::size_t i = 0; i < bin_count; ++i) {
    const std::size_t start = i * bin_width + 1;
    if (start > reference_length) {
      break;
    }
    const std::size_t end =
        i == bin_count - 1 ? reference_length
                           : std::min(reference_length, start + bin_width - 1U);
    const auto first =
        site_depth.begin() + static_cast<std::ptrdiff_t>(start - 1U);
    const auto last = site_depth.begin() + static_cast<std::ptrdiff_t>(end);
    const long double bin_total = std::accumulate(
        first, last, 0.0L, [](long double sum, std::uint64_t depth) {
          return sum + static_cast<long double>(depth);
        });
    const auto site_count = end - start + 1U;
    const auto mean_depth = static_cast<std::size_t>(
        (bin_total / static_cast<long double>(site_count)) + 0.5L);
    bins.push_back({start, end, mean_depth});
  }

  return result;
}

[[nodiscard]] EvidenceAggregationResult
aggregate_snps(const MoleculeAssemblyResult &assembly,
               const std::vector<ReadFeature> &features,
               const AnalysisConfig &config, std::size_t reference_length,
               const std::string &reference,
               const std::map<std::string, SnpCall> &clinical_annotations) {
  EvidenceAggregationResult result;
  auto &aggregates = result.variants;
  const auto &reads = assembly.representatives;
  const bool capture_sparse_store = config.result_schema == ResultSchema::v0_6;
  std::unordered_set<std::size_t> target_positions;
  for (std::size_t feature_index = 0U; feature_index < features.size();
       ++feature_index) {
    const auto &feature = features[feature_index];
    if (feature.filtered_numt ||
        (capture_sparse_store &&
         (feature_index >= assembly.molecules.size() ||
          !assembly.molecules[feature_index].analysis_eligible))) {
      continue;
    }
    for (const auto &snp : feature.snps) {
      const auto key = snp_key(snp.position, snp.reference, snp.alternate);
      auto [it, inserted] = aggregates.try_emplace(key);
      if (inserted) {
        it->second.call = snp;
      }
      target_positions.insert(snp.position);
    }
  }
  AnalysisConfig fragment_config = config;
  fragment_config.excluded_snp_flags =
      static_cast<std::uint16_t>(fragment_config.excluded_snp_flags & ~0x800U);
  if (capture_sparse_store) {
    for (std::size_t molecule_index = 0U;
         molecule_index < assembly.molecules.size(); ++molecule_index) {
      if (molecule_index >= features.size() ||
          features[molecule_index].filtered_numt ||
          !assembly.molecules[molecule_index].analysis_eligible) {
        continue;
      }
      for (const auto fragment_id :
           assembly.molecules[molecule_index].fragment_ids) {
        const auto &fragment = assembly.fragments[fragment_id.value];
        const auto &read = source_record_for_fragment(assembly, fragment);
        auto fragment_snps =
            call_snps_from_alignment(read, reference, fragment_config);
        apply_clinical_annotations(fragment_snps, clinical_annotations);
        for (const auto &snp : fragment_snps) {
          const auto key = snp_key(snp.position, snp.reference, snp.alternate);
          auto [it, inserted] = aggregates.try_emplace(key);
          if (inserted) {
            it->second.call = snp;
          }
          target_positions.insert(snp.position);
        }
      }
    }
  }

  struct EventAggregateBinding {
    SnpAggregate *aggregate = nullptr;
    std::size_t event_index = 0;
  };
  std::unordered_map<std::size_t, std::vector<EventAggregateBinding>>
      aggregates_by_position;
  aggregates_by_position.reserve(target_positions.size());
  if (capture_sparse_store) {
    result.store.events.reserve(aggregates.size());
  }
  for (auto &[key, aggregate] : aggregates) {
    const std::size_t event_index =
        capture_sparse_store ? result.store.events.size() : 0U;
    if (capture_sparse_store) {
      EvidenceEvent event;
      event.id = "snv:" + key;
      event.type = "SNV";
      event.start = aggregate.call.position;
      event.end = aggregate.call.position;
      event.length = 1U;
      event.reference.assign(1U, aggregate.call.reference);
      event.alternate.assign(1U, aggregate.call.alternate);
      event.normalization = "rcrs_circular_snv_v1";
      event.source_projection = "variants";
      event.negative_evidence_rule = "callable_base_allele";
      event.absence_assessable = true;
      result.store.events.push_back(std::move(event));
    }
    aggregates_by_position[aggregate.call.position].push_back(
        {&aggregate, event_index});
  }

  std::unordered_map<std::size_t, LocusAlleleAccumulator> allele_evidence;
  allele_evidence.reserve(target_positions.size());
  for (std::size_t read_index = 0; read_index < reads.size(); ++read_index) {
    throw_if_cancelled(config);
    if (read_index >= features.size() || features[read_index].filtered_numt ||
        (capture_sparse_store &&
         !assembly.molecules[read_index].analysis_eligible)) {
      continue;
    }
    std::map<std::size_t, MoleculeAlleleObservation> molecule_alleles;
    std::vector<std::pair<const ReadRecord *, AlignmentFragmentId>>
        evidence_records;
    if (capture_sparse_store) {
      const auto &molecule = assembly.molecules[read_index];
      evidence_records.reserve(molecule.fragment_ids.size());
      for (const auto fragment_id : molecule.fragment_ids) {
        const auto &fragment = assembly.fragments[fragment_id.value];
        evidence_records.emplace_back(
            &source_record_for_fragment(assembly, fragment), fragment.id);
      }
    } else {
      evidence_records.emplace_back(
          &reads[read_index],
          assembly.molecules[read_index].representative_fragment_id);
    }

    for (const auto &[read_pointer, fragment_id] : evidence_records) {
      const auto &read = *read_pointer;
      const auto &effective_config =
          capture_sparse_store ? fragment_config : config;
      if (!passes_snp_alignment_filters(read, effective_config) ||
          read.cigar_operations.empty() || reference_length == 0U) {
        continue;
      }
      std::size_t query_cursor = 0U;
      std::size_t reference_cursor =
          ((read.reference_start - 1U) % reference_length) + 1U;
      for (const auto &operation : read.cigar_operations) {
        const auto len = operation.length;
        const char op = operation.code;
        if (op == 'M' || op == '=' || op == 'X') {
          for (std::size_t offset = 0U;
               offset < len && query_cursor + offset < read.sequence.size();
               ++offset) {
            const std::size_t position =
                ((reference_cursor - 1U + offset) % reference_length) + 1U;
            if (!target_positions.contains(position)) {
              continue;
            }
            auto &molecule_observation = molecule_alleles[position];
            if (!molecule_observation.covered ||
                fragment_id < molecule_observation.alignment_fragment_id) {
              molecule_observation.alignment_fragment_id = fragment_id;
              molecule_observation.mapping_quality = read.mapping_quality;
              molecule_observation.strand =
                  (read.flags & 0x10U) == 0U ? '+' : '-';
            }
            molecule_observation.covered = true;
            const auto quality = phred_quality_at(read, query_cursor + offset);
            const char base =
                static_cast<char>(std::toupper(static_cast<unsigned char>(
                    read.sequence[query_cursor + offset])));
            if (!quality || *quality < config.min_base_quality ||
                !base_index(base)) {
              continue;
            }
            const auto query_index = query_cursor + offset;
            const double center_proximity =
                read.sequence.size() <= 1U
                    ? 0.0
                    : (2.0 * static_cast<double>(
                                 std::min(query_index, read.sequence.size() -
                                                           1U - query_index))) /
                          static_cast<double>(read.sequence.size() - 1U);
            if (!molecule_observation.has_passing_observation) {
              molecule_observation.base = base;
              molecule_observation.quality = *quality;
              molecule_observation.center_proximity = center_proximity;
              molecule_observation.query_index = query_index;
              molecule_observation.alignment_fragment_id = fragment_id;
              molecule_observation.mapping_quality = read.mapping_quality;
              molecule_observation.strand =
                  (read.flags & 0x10U) == 0U ? '+' : '-';
              molecule_observation.has_passing_observation = true;
            } else if (!molecule_observation.conflicted) {
              if (molecule_observation.base != base) {
                molecule_observation.conflicted = true;
                molecule_observation.base = 'N';
              } else if (*quality > molecule_observation.quality ||
                         (*quality == molecule_observation.quality &&
                          (fragment_id <
                               molecule_observation.alignment_fragment_id ||
                           (fragment_id ==
                                molecule_observation.alignment_fragment_id &&
                            query_index < molecule_observation.query_index)))) {
                molecule_observation.quality = *quality;
                molecule_observation.center_proximity = center_proximity;
                molecule_observation.query_index = query_index;
                molecule_observation.alignment_fragment_id = fragment_id;
                molecule_observation.mapping_quality = read.mapping_quality;
                molecule_observation.strand =
                    (read.flags & 0x10U) == 0U ? '+' : '-';
              }
            }
          }
          query_cursor += len;
          reference_cursor =
              ((reference_cursor - 1U + (len % reference_length)) %
               reference_length) +
              1U;
        } else if (op == 'I' || op == 'S') {
          query_cursor += len;
        } else if (op == 'D' || op == 'N') {
          reference_cursor =
              ((reference_cursor - 1U + (len % reference_length)) %
               reference_length) +
              1U;
        }
      }
    }

    for (const auto &[position, observation] : molecule_alleles) {
      const auto position_aggregates = aggregates_by_position.find(position);
      if (position_aggregates == aggregates_by_position.end()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "target SNP position has no normalized event");
      }
      if (capture_sparse_store) {
        if (read_index >= assembly.molecules.size()) {
          throw AnalysisError(
              AnalysisErrorCode::internal_error,
              "molecule assembly and representative order diverged");
        }
        const auto &molecule = assembly.molecules[read_index];
        for (const auto &binding : position_aggregates->second) {
          if (result.store.observations.size() >=
              config.max_evidence_observations) {
            throw AnalysisError(
                AnalysisErrorCode::resource_exhausted,
                "schema 0.6 sparse observation limit exceeded "
                "(max_evidence_observations=" +
                    std::to_string(config.max_evidence_observations) + ")");
          }
          EvidenceObservation stored;
          stored.id = result.store.observations.size();
          stored.molecule_index = molecule.index;
          stored.event_index = binding.event_index;
          stored.alignment_fragment_id = observation.alignment_fragment_id;
          stored.mapping_quality = observation.mapping_quality;
          stored.strand = observation.strand;
          stored.evidence_source = "aligned_base";
          if (observation.conflicted) {
            stored.state = ObservationState::conflict;
          } else if (!observation.has_passing_observation) {
            stored.state = ObservationState::low_quality;
          } else {
            stored.observed_allele.assign(1U, observation.base);
            stored.base_quality = observation.quality;
            stored.center_proximity = observation.center_proximity;
            if (observation.base == binding.aggregate->call.alternate) {
              stored.state = ObservationState::alternate;
            } else if (observation.base == binding.aggregate->call.reference) {
              stored.state = ObservationState::reference;
            } else {
              stored.state = ObservationState::event_absent;
            }
          }
          result.store.observations.push_back(std::move(stored));
        }
      }
      if (observation.conflicted || !observation.has_passing_observation) {
        continue;
      }
      if (const auto index = base_index(observation.base)) {
        const std::size_t strand_index = observation.strand == '+' ? 0U : 1U;
        auto &evidence = allele_evidence[position];
        ++evidence.depths[*index];
        ++evidence.strand_depths[*index][strand_index];
        evidence.quality_sums[*index] +=
            static_cast<double>(observation.quality);
        evidence.quality_mins[*index] =
            std::min(evidence.quality_mins[*index], observation.quality);
        evidence.quality_maxes[*index] =
            std::max(evidence.quality_maxes[*index], observation.quality);
        evidence.read_position_sums[*index] += observation.center_proximity;
        for (const auto &binding : position_aggregates->second) {
          if (observation.base == binding.aggregate->call.alternate) {
            binding.aggregate->supporting_reads.push_back(
                features[read_index].id);
          }
        }
      }
    }
  }

  constexpr double z = 1.959963984540054;
  for (auto &[_, aggregate] : aggregates) {
    auto &support = aggregate.supporting_reads;
    std::sort(support.begin(), support.end());
    support.erase(std::unique(support.begin(), support.end()), support.end());
    const auto evidence_it = allele_evidence.find(aggregate.call.position);
    if (evidence_it == allele_evidence.end()) {
      continue;
    }
    const auto &evidence = evidence_it->second;
    const auto &depths = evidence.depths;
    aggregate.callable_depth =
        std::accumulate(depths.begin(), depths.end(), std::size_t{0});
    const auto alternate = base_index(aggregate.call.alternate);
    const auto reference_allele = base_index(aggregate.call.reference);
    if (!alternate || !reference_allele) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "aggregate SNP contains a non-canonical allele");
    }
    aggregate.alternate_depth = depths[*alternate];
    aggregate.reference_depth = depths[*reference_allele];
    aggregate.other_depth = aggregate.callable_depth -
                            aggregate.alternate_depth -
                            aggregate.reference_depth;
    aggregate.alternate_strand_depths = evidence.strand_depths[*alternate];
    aggregate.reference_strand_depths =
        evidence.strand_depths[*reference_allele];
    aggregate.alternate_quality_sum = evidence.quality_sums[*alternate];
    aggregate.reference_quality_sum = evidence.quality_sums[*reference_allele];
    aggregate.alternate_quality_min = evidence.quality_mins[*alternate];
    aggregate.reference_quality_min = evidence.quality_mins[*reference_allele];
    aggregate.alternate_quality_max = evidence.quality_maxes[*alternate];
    aggregate.reference_quality_max = evidence.quality_maxes[*reference_allele];
    aggregate.alternate_read_position_sum =
        evidence.read_position_sums[*alternate];
    aggregate.reference_read_position_sum =
        evidence.read_position_sums[*reference_allele];
    for (std::size_t index = 0U; index < depths.size(); ++index) {
      if (index == *alternate || index == *reference_allele) {
        continue;
      }
      for (std::size_t strand = 0U; strand < 2U; ++strand) {
        aggregate.other_strand_depths[strand] +=
            evidence.strand_depths[index][strand];
      }
      aggregate.other_quality_sum += evidence.quality_sums[index];
      aggregate.other_read_position_sum += evidence.read_position_sums[index];
      if (depths[index] != 0U) {
        aggregate.other_quality_min =
            std::min(aggregate.other_quality_min, evidence.quality_mins[index]);
        aggregate.other_quality_max = std::max(aggregate.other_quality_max,
                                               evidence.quality_maxes[index]);
      }
    }
    if (aggregate.callable_depth == 0) {
      continue;
    }
    const double n = static_cast<double>(aggregate.callable_depth);
    const double proportion =
        static_cast<double>(aggregate.alternate_depth) / n;
    aggregate.heteroplasmy = proportion;
    const double z2 = z * z;
    const double denominator = 1.0 + z2 / n;
    const double center = (proportion + z2 / (2.0 * n)) / denominator;
    const double margin = z *
                          std::sqrt((proportion * (1.0 - proportion) / n) +
                                    (z2 / (4.0 * n * n))) /
                          denominator;
    aggregate.ci95_low = std::max(0.0, center - margin);
    aggregate.ci95_high = std::min(1.0, center + margin);
  }

  if (capture_sparse_store) {
    if (assembly.molecules.size() != reads.size()) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 molecule/representative cardinality mismatch");
    }
    std::vector<std::size_t> fragment_reference_counts(
        assembly.fragments.size(), 0U);
    std::set<std::string> molecule_ids;
    for (std::size_t molecule_index = 0;
         molecule_index < assembly.molecules.size(); ++molecule_index) {
      const auto &molecule = assembly.molecules[molecule_index];
      if (molecule.index.value != molecule_index ||
          !molecule_ids.insert(molecule.id).second) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 molecule IDs are not unique and contiguous");
      }
      if (molecule.representative_fragment_id.value >=
          assembly.fragments.size()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 representative fragment reference is unresolved");
      }
      for (const auto fragment_id : molecule.fragment_ids) {
        if (fragment_id.value >= assembly.fragments.size() ||
            assembly.fragments[fragment_id.value].molecule_index !=
                molecule.index) {
          throw AnalysisError(
              AnalysisErrorCode::internal_error,
              "schema 0.6 molecule fragment reference is unresolved");
        }
        ++fragment_reference_counts[fragment_id.value];
      }
    }
    if (std::any_of(fragment_reference_counts.begin(),
                    fragment_reference_counts.end(),
                    [](const std::size_t count) { return count != 1U; })) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 alignment fragment must resolve to exactly one molecule");
    }

    std::set<std::string> event_ids;
    for (std::size_t event_index = 0; event_index < result.store.events.size();
         ++event_index) {
      if (!event_ids.insert(result.store.events[event_index].id).second) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "schema 0.6 normalized event IDs are not unique");
      }
    }

    struct ProjectionCounts {
      std::size_t alternate = 0;
      std::size_t reference = 0;
      std::size_t other = 0;
      std::vector<std::string> supporting_molecules;
    };
    std::vector<ProjectionCounts> projected(result.store.events.size());
    std::set<std::pair<std::size_t, std::size_t>> molecule_event_pairs;
    for (std::size_t observation_index = 0;
         observation_index < result.store.observations.size();
         ++observation_index) {
      const auto &observation = result.store.observations[observation_index];
      if (observation.id != observation_index ||
          observation.molecule_index.value >= assembly.molecules.size() ||
          observation.event_index >= result.store.events.size() ||
          observation.alignment_fragment_id.value >=
              assembly.fragments.size()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "schema 0.6 observation reference is unresolved");
      }
      if (assembly.fragments[observation.alignment_fragment_id.value]
              .molecule_index != observation.molecule_index) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 observation fragment belongs to another molecule");
      }
      if (!molecule_event_pairs
               .emplace(observation.molecule_index.value,
                        observation.event_index)
               .second) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 molecule/event observation is not unique");
      }
      auto &counts = projected[observation.event_index];
      switch (observation.state) {
      case ObservationState::alternate:
        ++counts.alternate;
        counts.supporting_molecules.push_back(
            assembly.molecules[observation.molecule_index.value].id);
        break;
      case ObservationState::reference:
        ++counts.reference;
        break;
      case ObservationState::event_absent:
        ++counts.other;
        break;
      case ObservationState::low_quality:
      case ObservationState::conflict:
        break;
      case ObservationState::not_callable:
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "sparse schema 0.6 store must encode NOT_CALLABLE implicitly");
      }
    }

    for (std::size_t event_index = 0; event_index < result.store.events.size();
         ++event_index) {
      const auto &event = result.store.events[event_index];
      const auto aggregate_it = aggregates.find(snp_key(
          event.start, event.reference.front(), event.alternate.front()));
      if (aggregate_it == aggregates.end()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "schema 0.6 event has no variant projection");
      }
      const auto &aggregate = aggregate_it->second;
      auto &counts = projected[event_index];
      std::sort(counts.supporting_molecules.begin(),
                counts.supporting_molecules.end());
      result.store.events[event_index].supporting_molecules =
          counts.supporting_molecules;
      if (counts.alternate != aggregate.alternate_depth ||
          counts.reference != aggregate.reference_depth ||
          counts.other != aggregate.other_depth ||
          counts.alternate + counts.reference + counts.other !=
              aggregate.callable_depth ||
          counts.supporting_molecules != aggregate.supporting_reads) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 evidence/variant projection invariant failed for " +
                event.id);
      }
    }
  }
  return result;
}

void append_sparse_observation(SparseEvidenceStore &store,
                               EvidenceObservation observation,
                               const AnalysisConfig &config) {
  if (store.observations.size() >= config.max_evidence_observations) {
    throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                        "schema 0.6 sparse observation limit exceeded "
                        "(max_evidence_observations=" +
                            std::to_string(config.max_evidence_observations) +
                            ")");
  }
  observation.id = store.observations.size();
  store.observations.push_back(std::move(observation));
}

[[nodiscard]] std::unordered_map<std::string, std::size_t>
molecule_index_by_id(const MoleculeAssemblyResult &assembly) {
  std::unordered_map<std::string, std::size_t> indices;
  indices.reserve(assembly.molecules.size());
  for (std::size_t index = 0U; index < assembly.molecules.size(); ++index) {
    if (!indices.emplace(assembly.molecules[index].id, index).second) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 molecule ID is not unique");
    }
  }
  return indices;
}

[[nodiscard]] std::size_t
circular_previous_position(const std::size_t position,
                           const std::size_t reference_length) {
  return position > 1U ? position - 1U : reference_length;
}

[[nodiscard]] std::size_t
circular_next_position(const std::size_t position,
                       const std::size_t reference_length) {
  return position < reference_length ? position + 1U : 1U;
}

void append_small_indel_evidence(SparseEvidenceStore &store,
                                 const MoleculeAssemblyResult &assembly,
                                 const std::vector<ReadFeature> &features,
                                 const CallabilityResult &callability,
                                 const std::string &reference,
                                 const AnalysisConfig &config) {
  struct Support {
    std::size_t molecule_index = 0U;
    AlignmentFragmentId alignment_fragment_id;
    std::uint8_t mapping_quality = 0U;
    char strand = '+';
    NormalizedSmallIndel call;
  };
  struct Aggregate {
    NormalizedSmallIndel representative;
    std::map<std::size_t, Support> supports;
  };
  std::map<std::string, Aggregate> aggregates;
  AnalysisConfig fragment_config = config;
  fragment_config.excluded_snp_flags =
      static_cast<std::uint16_t>(fragment_config.excluded_snp_flags & ~0x800U);
  for (std::size_t molecule_index = 0U;
       molecule_index < assembly.molecules.size(); ++molecule_index) {
    throw_if_cancelled(config);
    if (molecule_index >= features.size() ||
        features[molecule_index].filtered_numt ||
        !assembly.molecules[molecule_index].analysis_eligible) {
      continue;
    }
    for (const auto fragment_id :
         assembly.molecules[molecule_index].fragment_ids) {
      const auto &fragment = assembly.fragments[fragment_id.value];
      const auto &read = source_record_for_fragment(assembly, fragment);
      auto calls =
          call_small_indels_from_alignment(read, fragment_config, reference);
      for (auto &call : calls) {
        auto [it, inserted] = aggregates.try_emplace(call.id);
        if (inserted) {
          it->second.representative = call;
        } else if (it->second.representative.type != call.type ||
                   it->second.representative.start != call.start ||
                   it->second.representative.end != call.end ||
                   it->second.representative.reference != call.reference ||
                   it->second.representative.alternate != call.alternate) {
          throw AnalysisError(AnalysisErrorCode::internal_error,
                              "normalized small-indel ID collision for '" +
                                  call.id + "'");
        }
        Support candidate{molecule_index, fragment_id, read.mapping_quality,
                          (read.flags & 0x10U) == 0U ? '+' : '-',
                          std::move(call)};
        const auto support = it->second.supports.find(molecule_index);
        if (support == it->second.supports.end() ||
            candidate.call.base_quality.value_or(0U) >
                support->second.call.base_quality.value_or(0U) ||
            (candidate.call.base_quality == support->second.call.base_quality &&
             candidate.alignment_fragment_id <
                 support->second.alignment_fragment_id)) {
          it->second.supports.insert_or_assign(molecule_index,
                                               std::move(candidate));
        }
      }
    }
  }

  using IndelLocus = std::pair<std::string, std::size_t>;
  std::map<IndelLocus,
           std::vector<std::pair<const std::string *, const Aggregate *>>>
      aggregates_by_locus;
  for (const auto &[event_id, aggregate] : aggregates) {
    aggregates_by_locus[{aggregate.representative.type,
                         aggregate.representative.start}]
        .emplace_back(&event_id, &aggregate);
  }

  for (auto &[event_id, aggregate] : aggregates) {
    const auto event_index = store.events.size();
    EvidenceEvent event;
    event.id = aggregate.representative.id;
    event.type = aggregate.representative.type;
    event.start = aggregate.representative.start;
    event.end = aggregate.representative.end;
    event.length = aggregate.representative.length;
    event.reference = aggregate.representative.reference;
    event.alternate = aggregate.representative.alternate;
    event.normalization = aggregate.representative.normalization;
    event.source_projection = "small_indels";
    event.absence_assessable =
        event.type == "SMALL_INSERTION" || event.end >= event.start;
    event.negative_evidence_rule =
        event.type == "SMALL_INSERTION"
            ? "same_fragment_callable_reference_adjacency"
            : "same_fragment_callable_deleted_span_with_flanks";
    for (const auto &[molecule_index, support] : aggregate.supports) {
      (void)support;
      event.supporting_molecules.push_back(
          assembly.molecules[molecule_index].id);
    }
    std::sort(event.supporting_molecules.begin(),
              event.supporting_molecules.end());
    event.supporting_molecules.erase(
        std::unique(event.supporting_molecules.begin(),
                    event.supporting_molecules.end()),
        event.supporting_molecules.end());
    store.events.push_back(event);

    std::vector<const Support *> support_by_molecule(assembly.molecules.size(),
                                                     nullptr);
    for (const auto &[molecule_index, support] : aggregate.supports) {
      support_by_molecule[molecule_index] = &support;
    }
    for (std::size_t molecule_index = 0U;
         molecule_index < assembly.molecules.size(); ++molecule_index) {
      if (molecule_index >= features.size() ||
          features[molecule_index].filtered_numt ||
          !assembly.molecules[molecule_index].analysis_eligible) {
        continue;
      }
      const auto &molecule = assembly.molecules[molecule_index];
      EvidenceObservation observation;
      observation.molecule_index = molecule.index;
      observation.event_index = event_index;
      observation.alignment_fragment_id = molecule.representative_fragment_id;
      const auto &read = assembly.representatives[molecule_index];
      observation.mapping_quality = read.mapping_quality;
      observation.strand = (read.flags & 0x10U) == 0U ? '+' : '-';
      if (const auto *support = support_by_molecule[molecule_index]) {
        observation.alignment_fragment_id = support->alignment_fragment_id;
        observation.mapping_quality = support->mapping_quality;
        observation.strand = support->strand;
        observation.state = ObservationState::alternate;
        observation.observed_allele = support->call.alternate;
        observation.base_quality = support->call.base_quality;
        observation.evidence_source = "cigar_small_indel";
        append_sparse_observation(store, std::move(observation), config);
        continue;
      }
      const Support *alternative_support = nullptr;
      const auto locus = aggregates_by_locus.find({event.type, event.start});
      if (locus == aggregates_by_locus.end()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "small-indel locus index is unresolved");
      }
      for (const auto &[alternative_id, alternative] : locus->second) {
        if (*alternative_id == event_id) {
          continue;
        }
        const auto found = alternative->supports.find(molecule_index);
        if (found != alternative->supports.end()) {
          alternative_support = &found->second;
          break;
        }
      }
      if (alternative_support != nullptr) {
        observation.alignment_fragment_id =
            alternative_support->alignment_fragment_id;
        observation.mapping_quality = alternative_support->mapping_quality;
        observation.strand = alternative_support->strand;
        observation.state = ObservationState::event_absent;
        observation.observed_allele = alternative_support->call.alternate;
        observation.base_quality = alternative_support->call.base_quality;
        observation.evidence_source = "cigar_alternative_small_indel";
        append_sparse_observation(store, std::move(observation), config);
        continue;
      }
      bool absent_is_callable = false;
      if (event.type == "SMALL_INSERTION") {
        absent_is_callable = reference_adjacency_callable(
            callability.molecules[molecule_index], event.start,
            circular_next_position(event.start, reference.size()));
      } else if (event.end >= event.start) {
        absent_is_callable = reference_span_with_flanks_callable(
            callability.molecules[molecule_index], event.start, event.end,
            reference.size());
      }
      if (absent_is_callable) {
        observation.state = ObservationState::event_absent;
        observation.observed_allele = event.reference;
        observation.evidence_source = "callable_reference_path";
        append_sparse_observation(store, std::move(observation), config);
      }
    }
  }
}

[[nodiscard]] std::string uppercase_ascii(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](const unsigned char character) {
                   return static_cast<char>(std::toupper(character));
                 });
  return value;
}

void append_structural_evidence(
    SparseEvidenceStore &store, const MoleculeAssemblyResult &assembly,
    const std::vector<ReadFeature> &features,
    const CallabilityResult &callability,
    const std::map<std::string, SvCall> &svs,
    const std::map<std::string, ComplexSvCall> &complex_events,
    const std::size_t reference_length, const AnalysisConfig &config) {
  const auto molecule_indices = molecule_index_by_id(assembly);
  for (const auto &[_, sv] : svs) {
    const auto event_index = store.events.size();
    EvidenceEvent event;
    event.id = "sv:" + sv.id;
    event.type = "SV_" + uppercase_ascii(sv.type);
    event.start = sv.start;
    event.end = sv.end;
    event.length = sv.length;
    event.normalization = "canonical_mtdna_adjacency_v1";
    event.source_projection = "svs";
    event.absence_assessable =
        sv.type == "insertion" || (sv.type == "deletion" && sv.end >= sv.start);
    event.negative_evidence_rule =
        sv.type == "insertion"
            ? "same_fragment_callable_reference_adjacency"
            : (sv.type == "deletion" && sv.end >= sv.start
                   ? "same_fragment_callable_deleted_span_with_flanks"
                   : "support_only_no_negative_inference");
    for (const auto &molecule_id : sv.supporting_reads) {
      const auto found = molecule_indices.find(molecule_id);
      if (found == molecule_indices.end()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "SV support references an unknown molecule '" +
                                molecule_id + "'");
      }
      const auto molecule_index = found->second;
      if (molecule_index < features.size() &&
          !features[molecule_index].filtered_numt &&
          assembly.molecules[molecule_index].analysis_eligible) {
        event.supporting_molecules.push_back(molecule_id);
      }
    }
    std::sort(event.supporting_molecules.begin(),
              event.supporting_molecules.end());
    event.supporting_molecules.erase(
        std::unique(event.supporting_molecules.begin(),
                    event.supporting_molecules.end()),
        event.supporting_molecules.end());
    if (event.supporting_molecules.empty()) {
      continue;
    }
    store.events.push_back(event);

    std::vector<bool> supports(assembly.molecules.size(), false);
    for (const auto &molecule_id : event.supporting_molecules) {
      const auto found = molecule_indices.find(molecule_id);
      if (found == molecule_indices.end()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "SV support references an unknown molecule '" +
                                molecule_id + "'");
      }
      supports[found->second] = true;
    }
    for (std::size_t molecule_index = 0U;
         molecule_index < assembly.molecules.size(); ++molecule_index) {
      if (molecule_index >= features.size() ||
          features[molecule_index].filtered_numt ||
          !assembly.molecules[molecule_index].analysis_eligible) {
        continue;
      }
      const auto &molecule = assembly.molecules[molecule_index];
      EvidenceObservation observation;
      observation.molecule_index = molecule.index;
      observation.event_index = event_index;
      observation.alignment_fragment_id = molecule.representative_fragment_id;
      const auto &read = assembly.representatives[molecule_index];
      observation.mapping_quality = read.mapping_quality;
      observation.strand = (read.flags & 0x10U) == 0U ? '+' : '-';
      if (supports[molecule_index]) {
        observation.state = ObservationState::alternate;
        observation.observed_allele = "<" + uppercase_ascii(sv.type) + ">";
        observation.evidence_source = "normalized_sv_support";
        append_sparse_observation(store, std::move(observation), config);
        continue;
      }

      bool absent_is_callable = false;
      if (sv.type == "insertion") {
        absent_is_callable = reference_adjacency_callable(
            callability.molecules[molecule_index], sv.start,
            circular_next_position(sv.start, reference_length));
      } else if (sv.type == "deletion" && sv.end >= sv.start) {
        absent_is_callable = reference_span_with_flanks_callable(
            callability.molecules[molecule_index], sv.start, sv.end,
            reference_length);
      }
      if (absent_is_callable) {
        observation.state = ObservationState::event_absent;
        observation.observed_allele = "REFERENCE_ADJACENCY";
        observation.evidence_source = "callable_reference_path";
        append_sparse_observation(store, std::move(observation), config);
      }
    }
  }

  for (const auto &[_, complex] : complex_events) {
    const auto event_index = store.events.size();
    EvidenceEvent event;
    event.id = complex.id;
    event.type = "COMPLEX_SV_PATH";
    event.normalization = "strand_invariant_ordered_path_v1";
    event.source_projection = "complex_events";
    event.negative_evidence_rule = "support_only_no_negative_inference";
    event.component_event_ids.reserve(complex.junction_ids.size());
    for (const auto &junction_id : complex.junction_ids) {
      event.component_event_ids.push_back("sv:" + junction_id);
    }
    for (const auto &molecule_id : complex.supporting_reads) {
      const auto found = molecule_indices.find(molecule_id);
      if (found == molecule_indices.end()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "complex-event support references an unknown molecule '" +
                molecule_id + "'");
      }
      const auto molecule_index = found->second;
      if (molecule_index < features.size() &&
          !features[molecule_index].filtered_numt &&
          assembly.molecules[molecule_index].analysis_eligible) {
        event.supporting_molecules.push_back(molecule_id);
      }
    }
    std::sort(event.supporting_molecules.begin(),
              event.supporting_molecules.end());
    event.supporting_molecules.erase(
        std::unique(event.supporting_molecules.begin(),
                    event.supporting_molecules.end()),
        event.supporting_molecules.end());
    if (event.supporting_molecules.empty()) {
      continue;
    }
    store.events.push_back(event);
    for (const auto &molecule_id : event.supporting_molecules) {
      const auto found = molecule_indices.find(molecule_id);
      if (found == molecule_indices.end()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "complex-event support references an unknown molecule '" +
                molecule_id + "'");
      }
      const auto molecule_index = found->second;
      const auto &molecule = assembly.molecules[molecule_index];
      EvidenceObservation observation;
      observation.molecule_index = molecule.index;
      observation.event_index = event_index;
      observation.alignment_fragment_id = molecule.representative_fragment_id;
      observation.state = ObservationState::alternate;
      observation.observed_allele = "COMPLEX_PATH";
      observation.evidence_source = "split_alignment_path";
      const auto &read = assembly.representatives[molecule_index];
      observation.mapping_quality = read.mapping_quality;
      observation.strand = (read.flags & 0x10U) == 0U ? '+' : '-';
      append_sparse_observation(store, std::move(observation), config);
    }
  }
}

void validate_unified_evidence_graph(const MoleculeAssemblyResult &assembly,
                                     const CallabilityResult &callability,
                                     const SparseEvidenceStore &store,
                                     const std::size_t reference_length) {
  if (callability.molecules.size() != assembly.molecules.size()) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "schema 0.6 callability/molecule cardinality mismatch");
  }
  for (const auto &molecule : assembly.molecules) {
    if (molecule.source_qnames.empty() || molecule.identity_policy.empty() ||
        (molecule.analysis_eligible && !molecule.exclusion_reasons.empty()) ||
        (!molecule.analysis_eligible && molecule.exclusion_reasons.empty())) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 molecule assembly contract is invalid");
    }
  }
  for (std::size_t index = 0U; index < callability.molecules.size(); ++index) {
    const auto &summary = callability.molecules[index];
    if (summary.molecule_index.value != index) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 callability molecule reference is unresolved");
    }
    std::size_t callable_bases = 0U;
    std::size_t previous_end = 0U;
    for (const auto &[start, end] : summary.ranges) {
      if (start == 0U || end < start || end > reference_length ||
          (previous_end != 0U && start <= previous_end)) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 callable ranges are invalid or overlapping");
      }
      callable_bases += end - start + 1U;
      previous_end = end;
    }
    if (callable_bases != summary.callable_bases ||
        (!summary.known && summary.callable_bases != 0U)) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 callable base count invariant failed");
    }
  }

  std::set<std::string> event_ids;
  for (const auto &event : store.events) {
    if (event.id.empty() || !event_ids.insert(event.id).second) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 normalized event IDs are empty or duplicated");
    }
  }
  std::vector<std::vector<std::string>> projected_support(store.events.size());
  std::set<std::pair<std::size_t, std::size_t>> molecule_event_pairs;
  for (std::size_t index = 0U; index < store.observations.size(); ++index) {
    const auto &observation = store.observations[index];
    if (observation.id != index ||
        observation.molecule_index.value >= assembly.molecules.size() ||
        observation.event_index >= store.events.size() ||
        observation.alignment_fragment_id.value >= assembly.fragments.size()) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 unified observation reference is unresolved");
    }
    if (assembly.fragments[observation.alignment_fragment_id.value]
            .molecule_index != observation.molecule_index) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 unified observation alignment belongs to "
                          "another molecule");
    }
    if (observation.state == ObservationState::not_callable) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 sparse store materialized NOT_CALLABLE");
    }
    if (observation.evidence_source.empty()) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 observation evidence source is empty");
    }
    if (!molecule_event_pairs
             .emplace(observation.molecule_index.value, observation.event_index)
             .second) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 unified molecule/event observation is not unique");
    }
    if (observation.state == ObservationState::alternate) {
      projected_support[observation.event_index].push_back(
          assembly.molecules[observation.molecule_index.value].id);
    }
  }
  for (std::size_t event_index = 0U; event_index < store.events.size();
       ++event_index) {
    const auto &event = store.events[event_index];
    if (event.negative_evidence_rule.empty() ||
        (event.absence_assessable == (event.negative_evidence_rule ==
                                      "support_only_no_negative_inference"))) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 event negative-evidence contract is invalid for " +
              event.id);
    }
    auto &support = projected_support[event_index];
    std::sort(support.begin(), support.end());
    support.erase(std::unique(support.begin(), support.end()), support.end());
    if (support != event.supporting_molecules) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 event/support projection invariant failed for " +
              event.id);
    }
  }
}

[[nodiscard]] bool is_callable_phase_state(const ObservationState state) {
  return state == ObservationState::reference ||
         state == ObservationState::alternate ||
         state == ObservationState::event_absent;
}

[[nodiscard]] std::vector<PhaseLink>
build_phase_links(const SparseEvidenceStore &store,
                  const MoleculeAssemblyResult &assembly,
                  const AnalysisConfig &config) {
  std::vector<std::vector<const EvidenceObservation *>> by_molecule(
      assembly.molecules.size());
  for (const auto &observation : store.observations) {
    if (assembly.molecules[observation.molecule_index.value]
            .analysis_eligible) {
      by_molecule[observation.molecule_index.value].push_back(&observation);
    }
  }
  for (auto &observations : by_molecule) {
    std::sort(observations.begin(), observations.end(),
              [](const auto *lhs, const auto *rhs) {
                return lhs->event_index < rhs->event_index;
              });
  }

  using Pair = std::pair<std::size_t, std::size_t>;
  std::set<Pair> candidate_pairs;
  for (const auto &observations : by_molecule) {
    for (std::size_t first = 0U; first < observations.size(); ++first) {
      if (observations[first]->state != ObservationState::alternate) {
        continue;
      }
      for (std::size_t second = 0U; second < observations.size(); ++second) {
        if (first == second) {
          continue;
        }
        const Pair pair{std::min(observations[first]->event_index,
                                 observations[second]->event_index),
                        std::max(observations[first]->event_index,
                                 observations[second]->event_index)};
        if (!candidate_pairs.contains(pair) &&
            candidate_pairs.size() >= config.max_phase_links) {
          throw AnalysisError(
              AnalysisErrorCode::resource_exhausted,
              "schema 0.6 phase-link limit exceeded (max_phase_links=" +
                  std::to_string(config.max_phase_links) + ")");
        }
        candidate_pairs.insert(pair);
      }
    }
  }

  std::map<Pair, PhaseLink> projected;
  for (const auto &pair : candidate_pairs) {
    PhaseLink link;
    link.event_a_index = pair.first;
    link.event_b_index = pair.second;
    link.complete_callability = store.events[pair.first].absence_assessable &&
                                store.events[pair.second].absence_assessable;
    projected.emplace(pair, link);
  }
  std::vector<std::vector<std::size_t>> candidate_neighbors(
      store.events.size());
  for (const auto &[event_a, event_b] : candidate_pairs) {
    candidate_neighbors[event_a].push_back(event_b);
  }
  for (std::size_t molecule_index = 0U; molecule_index < by_molecule.size();
       ++molecule_index) {
    const auto &observations = by_molecule[molecule_index];
    for (const auto *observation_a : observations) {
      for (const auto event_b_index :
           candidate_neighbors[observation_a->event_index]) {
        const auto observation_b = std::lower_bound(
            observations.begin(), observations.end(), event_b_index,
            [](const auto *observation, const std::size_t event_index) {
              return observation->event_index < event_index;
            });
        if (observation_b == observations.end() ||
            (*observation_b)->event_index != event_b_index) {
          continue;
        }
        const Pair pair{observation_a->event_index, event_b_index};
        const auto found = projected.find(pair);
        if (found == projected.end()) {
          throw AnalysisError(AnalysisErrorCode::internal_error,
                              "phase candidate adjacency is unresolved");
        }
        auto &link = found->second;
        const auto state_a = observation_a->state;
        const auto state_b = (*observation_b)->state;
        if (!is_callable_phase_state(state_a) ||
            !is_callable_phase_state(state_b)) {
          ++link.jointly_uncertain;
          link.uncertain_molecule_indices.push_back(molecule_index);
          continue;
        }
        ++link.jointly_callable;
        const bool alternate_a = state_a == ObservationState::alternate;
        const bool alternate_b = state_b == ObservationState::alternate;
        if (alternate_a && alternate_b) {
          ++link.both_alternate;
          link.supporting_molecule_indices.push_back(molecule_index);
        } else if (alternate_a) {
          ++link.a_alternate_b_absent;
        } else if (alternate_b) {
          ++link.a_absent_b_alternate;
        } else {
          ++link.neither_alternate;
        }
      }
    }
  }

  constexpr double z = 1.959963984540054;
  std::vector<PhaseLink> links;
  links.reserve(projected.size());
  for (auto &[_, link] : projected) {
    std::sort(link.supporting_molecule_indices.begin(),
              link.supporting_molecule_indices.end());
    std::sort(link.uncertain_molecule_indices.begin(),
              link.uncertain_molecule_indices.end());
    if (link.jointly_callable != 0U) {
      const double n = static_cast<double>(link.jointly_callable);
      const double both = static_cast<double>(link.both_alternate);
      link.co_alternate_fraction = both / n;
      const double z2 = z * z;
      const double denominator = 1.0 + z2 / n;
      const double center =
          (link.co_alternate_fraction + z2 / (2.0 * n)) / denominator;
      const double margin = z *
                            std::sqrt((link.co_alternate_fraction *
                                       (1.0 - link.co_alternate_fraction) / n) +
                                      (z2 / (4.0 * n * n))) /
                            denominator;
      link.ci95_low = std::max(0.0, center - margin);
      link.ci95_high = std::min(1.0, center + margin);
      const double p_a =
          static_cast<double>(link.both_alternate + link.a_alternate_b_absent) /
          n;
      const double p_b =
          static_cast<double>(link.both_alternate + link.a_absent_b_alternate) /
          n;
      link.expected_co_alternate_fraction = p_a * p_b;
      link.linkage_delta =
          link.co_alternate_fraction - link.expected_co_alternate_fraction;
    }
    links.push_back(link);
  }
  return links;
}

[[nodiscard]] std::vector<HaplogroupDefinition> load_haplogroups() {
  std::ifstream input(phylotree_path());
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open PhyloTree resource: " +
                            phylotree_path());
  }
  std::vector<HaplogroupDefinition> definitions;
  std::vector<HaplogroupDefinition> stack;
  std::string line;
  std::size_t mutation_count = 0;
  std::size_t source_order = 0;
  while (std::getline(input, line)) {
    const auto opening = line.find("<haplogroup name=\"");
    if (opening != std::string::npos) {
      const auto name_start =
          opening + std::string_view("<haplogroup name=\"").size();
      const auto name_end = line.find('"', name_start);
      if (name_end == std::string::npos) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "malformed haplogroup name in PhyloTree resource");
      }
      HaplogroupDefinition node;
      node.name = line.substr(name_start, name_end - name_start);
      node.source_order = source_order++;
      if (!stack.empty()) {
        node.mutations = stack.back().mutations;
      }
      stack.push_back(std::move(node));
    }

    const auto poly_start = line.find("<poly>");
    if (poly_start != std::string::npos) {
      if (stack.empty()) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "PhyloTree mutation is outside a haplogroup");
      }
      const auto value_start = poly_start + std::string_view("<poly>").size();
      const auto value_end = line.find("</poly>", value_start);
      if (value_end == std::string::npos) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "unterminated PhyloTree mutation element");
      }
      const auto encoded = trim_copy(
          std::string_view(line).substr(value_start, value_end - value_start));
      const auto mutation = parse_phylo_mutation(encoded);
      if (!mutation) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "unsupported PhyloTree mutation token: " + encoded);
      }
      ++mutation_count;
      auto mutation_it = std::lower_bound(
          stack.back().mutations.begin(), stack.back().mutations.end(),
          mutation->locus, [](const auto &item, const auto &locus) {
            return item.locus < locus;
          });
      const bool same_locus = mutation_it != stack.back().mutations.end() &&
                              mutation_it->locus == mutation->locus;
      if (mutation->backmutation) {
        if (same_locus) {
          stack.back().mutations.erase(mutation_it);
        }
      } else if (same_locus) {
        *mutation_it = *mutation;
      } else {
        stack.back().mutations.insert(mutation_it, *mutation);
      }
    }

    if (line.find("</haplogroup>") != std::string::npos) {
      if (stack.empty()) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "unbalanced PhyloTree haplogroup closing tag");
      }
      definitions.push_back(stack.back());
      stack.pop_back();
    }
  }
  if (!stack.empty() || definitions.size() < 5000U || mutation_count == 0U) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "PhyloTree resource is incomplete or unbalanced");
  }
  return definitions;
}

[[nodiscard]] std::unordered_map<std::string, double> load_phylo_weights() {
  std::ifstream input(phylotree_weights_path());
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open PhyloTree weights: " +
                            phylotree_weights_path());
  }
  std::unordered_map<std::string, double> weights;
  std::string line;
  while (std::getline(input, line)) {
    const auto fields = split_tab(line);
    if (fields.size() < 2) {
      continue;
    }
    try {
      const double weight = std::stod(fields[1]);
      if (std::isfinite(weight) && weight > 0.0) {
        weights[fields[0]] = weight;
      }
    } catch (const std::exception &) {
      // Ignore a malformed optional weight and use the deterministic default.
    }
  }
  return weights;
}

[[nodiscard]] double
mutation_weight(const std::unordered_map<std::string, double> &weights,
                const std::string &mutation) {
  const auto it = weights.find(mutation);
  return it == weights.end() ? 1.0 : it->second;
}

[[nodiscard]] std::string macro_haplogroup(std::string_view name) {
  if (name.empty()) {
    return {};
  }
  std::string macro(1, name.front());
  if (name.front() == 'L' && name.size() > 1 &&
      std::isdigit(static_cast<unsigned char>(name[1])) != 0) {
    macro.push_back(name[1]);
  }
  return macro;
}

[[nodiscard]] bool
phylo_mutation_matches(const PhyloMutation &expected,
                       const ObservedPhyloMutation &observed) {
  if (expected.locus.position != observed.position ||
      expected.kind != observed.kind) {
    return false;
  }
  if (expected.kind == PhyloMutationKind::substitution) {
    return expected.alternate == observed.alternate;
  }
  if (expected.kind == PhyloMutationKind::deletion) {
    return true;
  }
  if (expected.inserted_bases.empty() || observed.inserted_bases.empty()) {
    return false;
  }
  if (expected.locus.insertion_index == "X") {
    for (std::size_t i = 0; i < observed.inserted_bases.size(); ++i) {
      if (observed.inserted_bases[i] !=
          expected.inserted_bases[i % expected.inserted_bases.size()]) {
        return false;
      }
    }
    return true;
  }

  const auto expected_index = parse_size(expected.locus.insertion_index);
  const auto observed_index = parse_size(observed.insertion_index);
  if (!expected_index || !observed_index || *expected_index < *observed_index) {
    return false;
  }
  const std::size_t offset = *expected_index - *observed_index;
  return offset <= observed.inserted_bases.size() &&
         expected.inserted_bases.size() <=
             observed.inserted_bases.size() - offset &&
         std::equal(expected.inserted_bases.begin(),
                    expected.inserted_bases.end(),
                    observed.inserted_bases.begin() +
                        static_cast<std::ptrdiff_t>(offset));
}

[[nodiscard]] bool
observed_contains_mutation(const std::vector<ObservedPhyloMutation> &observed,
                           const PhyloMutation &expected) {
  auto candidate = std::lower_bound(observed.begin(), observed.end(),
                                    expected.locus.position,
                                    [](const auto &item, std::size_t position) {
                                      return item.position < position;
                                    });
  while (candidate != observed.end() &&
         candidate->position == expected.locus.position) {
    if (phylo_mutation_matches(expected, *candidate)) {
      return true;
    }
    ++candidate;
  }
  return false;
}

[[nodiscard]] bool
definition_contains_observed(const std::vector<PhyloMutation> &expected,
                             const ObservedPhyloMutation &observed) {
  auto candidate =
      std::lower_bound(expected.begin(), expected.end(), observed.position,
                       [](const auto &item, std::size_t position) {
                         return item.locus.position < position;
                       });
  while (candidate != expected.end() &&
         candidate->locus.position == observed.position) {
    if (phylo_mutation_matches(*candidate, observed)) {
      return true;
    }
    ++candidate;
  }
  return false;
}

[[nodiscard]] bool position_in_ranges(
    const std::vector<std::pair<std::size_t, std::size_t>> &ranges,
    std::size_t position) {
  const auto candidate =
      std::lower_bound(ranges.begin(), ranges.end(), position,
                       [](const auto &range, std::size_t value) {
                         return range.second < value;
                       });
  return candidate != ranges.end() && candidate->first <= position;
}

template <typename Visitor>
void for_each_mutation_in_ranges(
    const std::vector<PhyloMutation> &mutations,
    const std::vector<std::pair<std::size_t, std::size_t>> &ranges,
    Visitor &&visitor) {
  for (const auto &[start, end] : ranges) {
    auto mutation =
        std::lower_bound(mutations.begin(), mutations.end(), start,
                         [](const auto &item, std::size_t position) {
                           return item.locus.position < position;
                         });
    while (mutation != mutations.end() && mutation->locus.position <= end) {
      visitor(*mutation);
      ++mutation;
    }
  }
}

[[nodiscard]] CallablePhyloRanges
majority_callable_ranges(const std::vector<const ReadFeature *> &members) {
  CallablePhyloRanges result;
  if (members.size() == 1U) {
    result.known = members.front()->haplogroup_range_known;
    if (result.known) {
      result.ranges = members.front()->haplogroup_ranges;
    } else {
      result.ranges.emplace_back(
          1U, static_cast<std::size_t>(kDefaultReferenceLength));
    }
    return result;
  }
  const std::size_t known_range_count = static_cast<std::size_t>(
      std::count_if(members.begin(), members.end(), [](const auto *item) {
        return item->haplogroup_range_known;
      }));
  result.known = known_range_count * 2U > members.size();
  if (!result.known) {
    result.ranges.emplace_back(
        1U, static_cast<std::size_t>(kDefaultReferenceLength));
    return result;
  }

  struct RangeEvent {
    std::size_t position = 0U;
    std::int64_t delta = 0;
  };
  std::vector<RangeEvent> events;
  std::size_t range_count = 0U;
  for (const auto *member : members) {
    const auto member_ranges =
        member->haplogroup_range_known ? member->haplogroup_ranges.size() : 0U;
    if (member_ranges > std::numeric_limits<std::size_t>::max() - range_count) {
      throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                          "haplogroup callable range count overflow");
    }
    range_count += member_ranges;
  }
  if (range_count > std::numeric_limits<std::size_t>::max() / 2U) {
    throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                        "haplogroup callable event count overflow");
  }
  events.reserve(range_count * 2U);
  for (const auto *member : members) {
    if (!member->haplogroup_range_known) {
      continue;
    }
    for (const auto &[start, end] : member->haplogroup_ranges) {
      if (start == 0U || end < start ||
          end > static_cast<std::size_t>(kDefaultReferenceLength)) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "haplogroup callable range is outside the reference");
      }
      events.push_back({start, 1});
      events.push_back({end + 1U, -1});
    }
  }
  std::sort(events.begin(), events.end(), [](const auto &lhs, const auto &rhs) {
    return std::tie(lhs.position, lhs.delta) <
           std::tie(rhs.position, rhs.delta);
  });

  const auto threshold = static_cast<std::int64_t>(members.size() / 2U);
  std::int64_t active = 0;
  std::size_t cursor = 1U;
  std::size_t event_index = 0U;
  while (cursor <= static_cast<std::size_t>(kDefaultReferenceLength)) {
    while (event_index < events.size() &&
           events[event_index].position == cursor) {
      active += events[event_index].delta;
      ++event_index;
    }
    if (active < 0) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "haplogroup callable range coverage underflow");
    }
    const std::size_t next =
        event_index < events.size()
            ? std::min(events[event_index].position,
                       static_cast<std::size_t>(kDefaultReferenceLength) + 1U)
            : static_cast<std::size_t>(kDefaultReferenceLength) + 1U;
    if (active > threshold && next > cursor) {
      result.ranges.emplace_back(cursor, next - 1U);
    }
    if (next <= cursor) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "haplogroup callable range events are not ordered");
    }
    cursor = next;
  }
  result.ranges = normalize_ranges(std::move(result.ranges));
  return result;
}

[[nodiscard]] std::map<int, ClusterHaplogroupAssignment>
assign_haplogroups(const std::vector<ReadFeature> &features,
                   const std::vector<HaplogroupDefinition> &definitions,
                   const std::unordered_map<std::string, double> &weights,
                   const std::vector<PhyloAlignmentRule> &alignment_rules,
                   const AnalysisConfig &config) {
  std::map<int, std::vector<const ReadFeature *>> clusters;
  for (const auto &feature : features) {
    if (!feature.filtered_numt) {
      clusters[feature.cluster_id].push_back(&feature);
    }
  }

  std::map<int, ClusterHaplogroupAssignment> assignments;
  std::unordered_map<std::size_t, std::vector<std::size_t>>
      definitions_by_position;
  definitions_by_position.reserve(
      static_cast<std::size_t>(kDefaultReferenceLength));
  for (std::size_t definition_index = 0; definition_index < definitions.size();
       ++definition_index) {
    for (const auto &mutation : definitions[definition_index].mutations) {
      definitions_by_position[mutation.locus.position].push_back(
          definition_index);
    }
  }
  std::vector<std::size_t> candidate_marks(definitions.size(), 0U);
  std::size_t candidate_epoch = 0;
  std::unordered_map<std::string, ClusterHaplogroupAssignment> assignment_cache;
  assignment_cache.reserve(clusters.size());
  for (const auto &[cluster_id, members] : clusters) {
    throw_if_cancelled(config);
    struct ObservationCount {
      ObservedPhyloMutation mutation;
      std::size_t count = 0;
    };
    std::unordered_map<std::string_view, ObservationCount> observed_counts;
    for (const auto *member : members) {
      std::unordered_set<std::string_view> molecule_mutations;
      molecule_mutations.reserve(member->haplogroup_markers.size());
      for (const auto &mutation : member->haplogroup_markers) {
        if (!molecule_mutations.insert(mutation.encoded).second) {
          continue;
        }
        auto [count, inserted] = observed_counts.try_emplace(
            mutation.encoded, ObservationCount{mutation, 0U});
        static_cast<void>(inserted);
        ++count->second.count;
      }
    }
    std::vector<ObservedPhyloMutation> observed;
    for (const auto &[_, counted] : observed_counts) {
      if (counted.count * 2U <= members.size()) {
        continue;
      }
      observed.push_back(counted.mutation);
    }
    std::sort(observed.begin(), observed.end(), observed_phylo_mutation_less);

    auto callable = majority_callable_ranges(members);
    if (callable.known) {
      observed.erase(std::remove_if(observed.begin(), observed.end(),
                                    [&](const auto &mutation) {
                                      return !position_in_ranges(
                                          callable.ranges, mutation.position);
                                    }),
                     observed.end());
    }
    apply_phylo_alignment_rules(observed, alignment_rules);
    if (callable.known) {
      observed.erase(std::remove_if(observed.begin(), observed.end(),
                                    [&](const auto &mutation) {
                                      return !position_in_ranges(
                                          callable.ranges, mutation.position);
                                    }),
                     observed.end());
    }
    const auto &callable_ranges = callable.ranges;
    std::string observed_signature;
    for (const auto &mutation : observed) {
      observed_signature.append(mutation.encoded);
      observed_signature.push_back('\0');
    }
    for (const auto &[start, end] : callable_ranges) {
      observed_signature.append("range:");
      observed_signature.append(std::to_string(start));
      observed_signature.push_back('-');
      observed_signature.append(std::to_string(end));
      observed_signature.push_back('\0');
    }
    if (const auto cached = assignment_cache.find(observed_signature);
        cached != assignment_cache.end()) {
      assignments.emplace(cluster_id, cached->second);
      continue;
    }

    struct RankedDefinition {
      const HaplogroupDefinition *definition = nullptr;
      double score = 0.0;
      std::size_t matched_count = 0;
    };
    double observed_weight = 0.0;
    for (const auto &mutation : observed) {
      observed_weight += mutation_weight(weights, mutation.encoded);
    }
    std::vector<RankedDefinition> ranked;
    std::vector<std::size_t> candidate_indices;
    if (candidate_epoch == std::numeric_limits<std::size_t>::max()) {
      std::fill(candidate_marks.begin(), candidate_marks.end(), 0U);
      candidate_epoch = 0;
    }
    ++candidate_epoch;
    for (const auto &mutation : observed) {
      const auto posting = definitions_by_position.find(mutation.position);
      if (posting == definitions_by_position.end()) {
        continue;
      }
      for (const auto definition_index : posting->second) {
        if (candidate_marks[definition_index] != candidate_epoch) {
          candidate_marks[definition_index] = candidate_epoch;
          candidate_indices.push_back(definition_index);
        }
      }
    }
    ranked.reserve(candidate_indices.size());
    for (const auto definition_index : candidate_indices) {
      const auto &definition = definitions[definition_index];
      double expected_weight = 0.0;
      double matched_weight = 0.0;
      std::size_t matched_count = 0;
      for_each_mutation_in_ranges(
          definition.mutations, callable_ranges, [&](const auto &mutation) {
            expected_weight += mutation.weight;
            if (observed_contains_mutation(observed, mutation)) {
              matched_weight += mutation.weight;
              ++matched_count;
            }
          });
      double score = 0.0;
      if (expected_weight == 0.0 && observed_weight == 0.0) {
        score = 100.0;
      } else if (expected_weight > 0.0 && observed_weight > 0.0) {
        score = 50.0 * ((matched_weight / expected_weight) +
                        (matched_weight / observed_weight));
      }
      if (matched_count > 0 && score > 0.0) {
        ranked.push_back({&definition, score, matched_count});
      }
    }
    const auto ranked_less = [](const auto &lhs, const auto &rhs) {
      if (lhs.score != rhs.score) {
        return lhs.score > rhs.score;
      }
      if (lhs.definition->source_order != rhs.definition->source_order) {
        return lhs.definition->source_order < rhs.definition->source_order;
      }
      return lhs.definition->name < rhs.definition->name;
    };
    const auto keep = std::min<std::size_t>(3, ranked.size());
    std::partial_sort(ranked.begin(),
                      ranked.begin() + static_cast<std::ptrdiff_t>(keep),
                      ranked.end(), ranked_less);
    ClusterHaplogroupAssignment assignment;
    assignment.callable_ranges = callable_ranges;
    assignment.observed_markers.reserve(observed.size());
    for (const auto &mutation : observed) {
      assignment.observed_markers.push_back(mutation.encoded);
    }
    for (std::size_t i = 0; i < keep; ++i) {
      HaplogroupCandidate candidate;
      candidate.name = ranked[i].definition->name;
      candidate.score = ranked[i].score;
      for_each_mutation_in_ranges(
          ranked[i].definition->mutations, callable_ranges,
          [&](const auto &mutation) {
            if (observed_contains_mutation(observed, mutation)) {
              candidate.matched.push_back(mutation.encoded);
            } else {
              candidate.missing.push_back(mutation.encoded);
            }
          });
      for (const auto &mutation : observed) {
        if (!definition_contains_observed(ranked[i].definition->mutations,
                                          mutation)) {
          candidate.extra.push_back(mutation.encoded);
        }
      }
      assignment.candidates.push_back(std::move(candidate));
    }
    if (!assignment.candidates.empty()) {
      assignment.best = assignment.candidates.front().name;
      assignment.quality = assignment.candidates.front().score;
    }
    if (assignment.candidates.size() > 1 &&
        assignment.candidates[1].score > 0.0 &&
        assignment.candidates.front().score - assignment.candidates[1].score <
            3.0 &&
        macro_haplogroup(assignment.candidates.front().name) !=
            macro_haplogroup(assignment.candidates[1].name)) {
      assignment.contamination_warning = true;
    }
    assignment_cache.emplace(observed_signature, assignment);
    assignments.emplace(cluster_id, std::move(assignment));
  }
  return assignments;
}

[[nodiscard]] std::string
alignment_fragment_id_string(const AlignmentFragmentId id) {
  return "alignment:" + std::to_string(id.value);
}

[[nodiscard]] std::string_view
observation_state_name(const ObservationState state) {
  switch (state) {
  case ObservationState::reference:
    return "REFERENCE";
  case ObservationState::alternate:
    return "ALTERNATE";
  case ObservationState::event_absent:
    return "EVENT_ABSENT";
  case ObservationState::not_callable:
    return "NOT_CALLABLE";
  case ObservationState::low_quality:
    return "LOW_QUALITY";
  case ObservationState::conflict:
    return "CONFLICT";
  }
  throw AnalysisError(AnalysisErrorCode::internal_error,
                      "unknown schema 0.6 observation state");
}

[[nodiscard]] bool is_unified_variant_event(const EvidenceEvent &event) {
  return event.type == "SNV" || event.type == "SMALL_INSERTION" ||
         event.type == "SMALL_DELETION";
}

void add_variant_allele_observation(
    VariantAlleleEvidenceSummary &summary,
    const EvidenceObservation &observation) {
  ++summary.count;
  ++summary.strand_depths[observation.strand == '+' ? 0U : 1U];
  summary.mapping_quality_sum +=
      static_cast<double>(observation.mapping_quality);
  summary.mapping_quality_min =
      std::min(summary.mapping_quality_min, observation.mapping_quality);
  summary.mapping_quality_max =
      std::max(summary.mapping_quality_max, observation.mapping_quality);
  if (observation.base_quality) {
    ++summary.base_quality_count;
    summary.base_quality_sum +=
        static_cast<double>(*observation.base_quality);
    summary.base_quality_min =
        std::min(summary.base_quality_min, *observation.base_quality);
    summary.base_quality_max =
        std::max(summary.base_quality_max, *observation.base_quality);
  }
  if (observation.center_proximity) {
    ++summary.read_position_count;
    summary.read_position_sum += *observation.center_proximity;
  }
}

[[nodiscard]] std::vector<UnifiedVariantEvidenceSummary>
summarize_unified_variant_evidence(const SparseEvidenceStore &store) {
  std::vector<UnifiedVariantEvidenceSummary> summaries(store.events.size());
  for (const auto &observation : store.observations) {
    if (observation.event_index >= store.events.size()) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "variant projection references an unknown event");
    }
    const auto &event = store.events[observation.event_index];
    if (!is_unified_variant_event(event)) {
      continue;
    }
    auto &summary = summaries[observation.event_index];
    switch (observation.state) {
    case ObservationState::alternate:
      add_variant_allele_observation(summary.alternate, observation);
      break;
    case ObservationState::reference:
      add_variant_allele_observation(summary.reference, observation);
      break;
    case ObservationState::event_absent:
      ++summary.event_absent;
      if (event.type == "SNV" ||
          observation.evidence_source == "cigar_alternative_small_indel") {
        add_variant_allele_observation(summary.other, observation);
      } else {
        add_variant_allele_observation(summary.reference, observation);
      }
      break;
    case ObservationState::low_quality:
      ++summary.low_quality;
      break;
    case ObservationState::conflict:
      ++summary.conflict;
      break;
    case ObservationState::not_callable:
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "unified variant projection received explicit NOT_CALLABLE");
    }
  }

  std::map<std::pair<std::string, std::size_t>, std::size_t> locus_counts;
  for (const auto &event : store.events) {
    if (is_unified_variant_event(event)) {
      ++locus_counts[{event.type, event.start}];
    }
  }
  for (std::size_t index = 0U; index < store.events.size(); ++index) {
    const auto &event = store.events[index];
    if (is_unified_variant_event(event)) {
      summaries[index].multi_allelic =
          locus_counts[{event.type, event.start}] > 1U;
    }
  }
  return summaries;
}

[[nodiscard]] std::pair<double, double>
wilson_interval(const std::size_t support, const std::size_t depth) {
  if (depth == 0U) {
    return {0.0, 0.0};
  }
  constexpr double z = 1.959963984540054;
  const double n = static_cast<double>(depth);
  const double proportion = static_cast<double>(support) / n;
  const double z2 = z * z;
  const double denominator = 1.0 + z2 / n;
  const double center = (proportion + z2 / (2.0 * n)) / denominator;
  const double margin =
      z * std::sqrt((proportion * (1.0 - proportion) / n) +
                    (z2 / (4.0 * n * n))) /
      denominator;
  return {std::max(0.0, center - margin),
          std::min(1.0, center + margin)};
}

[[nodiscard]] std::size_t
circular_reference_homopolymer_run(const std::string &reference,
                                   const std::size_t position) {
  if (reference.empty() || position == 0U || position > reference.size()) {
    return 0U;
  }
  const char base = reference[position - 1U];
  std::size_t run = 1U;
  for (std::size_t offset = 1U; offset < reference.size(); ++offset) {
    const auto index =
        (position - 1U + reference.size() - offset) % reference.size();
    if (reference[index] != base) {
      break;
    }
    ++run;
  }
  for (std::size_t offset = 1U; run < reference.size(); ++offset) {
    const auto index = (position - 1U + offset) % reference.size();
    if (reference[index] != base) {
      break;
    }
    ++run;
  }
  return run;
}

void write_mapping_quality_summary(JsonBuffer &out,
                                   const VariantAlleleEvidenceSummary &summary) {
  out << "{\"count\":" << summary.count << ",\"mean\":";
  if (summary.count == 0U) {
    out << "null,\"min\":null,\"max\":null}";
    return;
  }
  out << (summary.mapping_quality_sum /
          static_cast<double>(summary.count))
      << ",\"min\":"
      << static_cast<unsigned int>(summary.mapping_quality_min)
      << ",\"max\":"
      << static_cast<unsigned int>(summary.mapping_quality_max) << "}";
}

void write_base_quality_summary(JsonBuffer &out,
                                const VariantAlleleEvidenceSummary &summary) {
  write_quality_summary(out, summary.base_quality_count,
                        summary.base_quality_sum, summary.base_quality_min,
                        summary.base_quality_max);
}

[[nodiscard]] std::vector<std::string> variant_qc_flags(
    const EvidenceEvent &event, const UnifiedVariantEvidenceSummary &summary,
    const std::size_t homopolymer_run, const bool numt_assessable) {
  std::vector<std::string> flags;
  if (!numt_assessable) {
    flags.emplace_back("NUMT_NOT_ASSESSABLE");
  }
  if (!event.absence_assessable) {
    flags.emplace_back("SUPPORT_ONLY_EVENT");
  }
  if (summary.conflict != 0U) {
    flags.emplace_back("CONFLICTING_MOLECULE");
  }
  if (summary.low_quality != 0U) {
    flags.emplace_back("LOW_QUALITY_OBSERVATIONS");
  }
  if (summary.alternate.count != 0U &&
      (summary.alternate.strand_depths[0] == 0U ||
       summary.alternate.strand_depths[1] == 0U)) {
    flags.emplace_back("SINGLE_STRAND_ALT_SUPPORT");
  }
  if (summary.multi_allelic) {
    flags.emplace_back("MULTI_ALLELIC_LOCUS");
  }
  if (homopolymer_run > 1U) {
    flags.emplace_back("HOMOPOLYMER_CONTEXT");
  }
  return flags;
}

[[nodiscard]] std::string
render_json(const std::string &input_path, const std::string &reference_path,
            const AnalysisConfig &config, std::string_view clustering_backend,
            bool alignment_input, bool has_nuclear_contigs,
            std::size_t input_record_count, const std::string &reference,
            const MoleculeAssemblyResult &assembly,
            const std::vector<ReadRecord> &reads,
            const std::vector<ReadFeature> &features,
            const std::map<std::string, SvCall> &svs,
            const std::map<std::string, ComplexSvCall> &complex_events,
            const CoverageResult &coverage_result,
            const std::map<std::string, SnpAggregate> &variants,
            const SparseEvidenceStore &evidence_store,
            const CallabilityResult &callability,
            const std::vector<PhaseLink> &phase_links,
            const std::map<int, ClusterHaplogroupAssignment> &haplogroups,
            const std::vector<ResourceRecord> &resources) {
  throw_if_cancelled(config);
  const auto &coverage = coverage_result.bins;
  std::map<int, std::vector<const ReadFeature *>> clusters;
  for (const auto &feature : features) {
    if (!feature.filtered_numt) {
      clusters[feature.cluster_id].push_back(&feature);
    }
  }

  std::vector<std::array<std::size_t, 5>> molecule_evidence_counts;
  std::vector<std::vector<std::string>> molecule_alternate_event_ids;
  if (config.result_schema == ResultSchema::v0_6) {
    molecule_evidence_counts.resize(assembly.molecules.size());
    molecule_alternate_event_ids.resize(assembly.molecules.size());
    for (const auto &observation : evidence_store.observations) {
      if (observation.molecule_index.value >= assembly.molecules.size() ||
          observation.event_index >= evidence_store.events.size()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "molecule projection received an unresolved observation");
      }
      auto &counts =
          molecule_evidence_counts[observation.molecule_index.value];
      switch (observation.state) {
      case ObservationState::alternate:
        ++counts[0];
        molecule_alternate_event_ids[observation.molecule_index.value]
            .push_back(evidence_store.events[observation.event_index].id);
        break;
      case ObservationState::reference:
        ++counts[1];
        break;
      case ObservationState::event_absent:
        ++counts[2];
        break;
      case ObservationState::low_quality:
        ++counts[3];
        break;
      case ObservationState::conflict:
        ++counts[4];
        break;
      case ObservationState::not_callable:
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "molecule projection received explicit NOT_CALLABLE");
      }
    }
  }

  const std::size_t read_snp_count =
      std::accumulate(features.begin(), features.end(), std::size_t{0},
                      [](const std::size_t count, const ReadFeature &feature) {
                        return count + feature.snps.size();
                      });
  constexpr std::size_t maximum_initial_reserve = 128U * 1024U * 1024U;
  std::size_t estimated_json_bytes = 64U * 1024U;
  const auto add_reserve_estimate = [&](const std::size_t count,
                                        const std::size_t bytes_each) {
    const std::size_t available =
        maximum_initial_reserve - estimated_json_bytes;
    estimated_json_bytes +=
        std::min(count, available / bytes_each) * bytes_each;
  };
  add_reserve_estimate(reads.size(), 512U);
  add_reserve_estimate(read_snp_count, 48U);
  add_reserve_estimate(variants.size(), 640U);
  JsonBuffer out(estimated_json_bytes);
  out << "{";
  out << "\"metadata\":{";
  out << "\"schema_version\":"
      << quoted(config.result_schema == ResultSchema::v0_6
                    ? kEvidenceGraphSchemaVersion
                    : kResultSchemaVersion)
      << ",";
  out << "\"sv_event_schema_version\":" << quoted(kSvEventSchemaVersion) << ",";
  out << "\"complex_sv_event_schema_version\":"
      << quoted(kComplexSvEventSchemaVersion) << ",";
  out << "\"clinical_annotation_schema_version\":"
      << quoted(kClinicalAnnotationSchemaVersion) << ",";
  out << "\"engine_version\":" << quoted(kEngineVersion) << ",";
  out << "\"sample\":"
      << quoted(config.sample_name.empty()
                    ? std::filesystem::path(input_path).stem().string()
                    : config.sample_name)
      << ",";
  out << "\"input_path\":" << quoted(input_path) << ",";
  out << "\"reference_path\":"
      << quoted(reference_path.empty() ? bundled_rcrs_path() : reference_path)
      << ",";
  out << "\"reference_accession\":"
      << quoted(reference_path.empty() ? "NC_012920.1" : "custom") << ",";
  out << "\"reference_length\":" << reference.size() << ",";
  out << "\"threads\":" << effective_thread_count(config.threads, reads.size())
      << ",";
  out << "\"requested_threads\":" << std::max<std::size_t>(1, config.threads)
      << ",";
  out << "\"calling_parameters\":{";
  out << "\"min_mapping_quality\":"
      << static_cast<unsigned int>(config.min_mapping_quality) << ",";
  out << "\"min_base_quality\":"
      << static_cast<unsigned int>(config.min_base_quality) << ",";
  out << "\"excluded_snp_flags\":" << config.excluded_snp_flags << ",";
  out << "\"numt_threshold\":" << config.numt_threshold << ",";
  out << "\"development_tags_enabled\":"
      << (config.allow_development_tags ? "true" : "false");
  if (config.result_schema == ResultSchema::v0_6) {
    out << ",\"max_evidence_observations\":" << config.max_evidence_observations
        << ","
        << "\"max_phase_links\":" << config.max_phase_links << ","
        << "\"evidence_page_size\":" << config.evidence_page_size << ","
        << "\"molecule_id_tag\":" << quoted(config.molecule_id_tag) << ","
        << "\"umi_tag\":" << quoted(config.umi_tag) << ","
        << "\"duplex_tag\":" << quoted(config.duplex_tag);
  }
  out << "},";
  out << "\"resources\":[";
  for (std::size_t i = 0; i < resources.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    const auto &resource = resources[i];
    out << "{\"name\":" << quoted(resource.name) << ","
        << "\"version\":" << quoted(resource.version) << ","
        << "\"path\":" << quoted(resource.path) << ","
        << "\"sha256\":" << quoted(resource.sha256) << ","
        << "\"source\":" << quoted(resource.source) << ","
        << "\"license\":" << quoted(resource.license) << ","
        << "\"retrieved\":" << quoted(resource.retrieved) << "}";
  }
  out << "],";
  out << "\"algorithm_notes\":["
      << quoted("FASTQ/SAM parser with optional htslib BAM/CRAM reader") << ","
      << quoted("Aligned reads use CIGAR-aware reference comparison for SNP "
                "extraction")
      << ","
      << quoted("NUMT specificity requires competitive "
                "nuclear-plus-mitochondrial alignment")
      << ","
      << quoted("Haplogroups use weighted PhyloTree 17.3 lineage matching")
      << ","
      << quoted("Allele QC metrics are observational and are not unvalidated "
                "hard filters")
      << ","
      << quoted("Clinical schema 1.0 preserves source assertions and "
                "deterministic conflicts")
      << ","
      << quoted("Clinical assertion payloads are emitted once per aggregate "
                "variant")
      << "," << quoted("SV event IDs use canonical mtDNA adjacency schema 1.0")
      << ","
      << quoted("Split and supplementary alignments are reconstructed as "
                "molecule event edges")
      << ","
      << quoted("Multi-junction split paths are coalesced with "
                "strand-invariant exact IDs")
      << "," << quoted(clustering_backend) << "]";
  out << "},";

  out << "\"filter_stats\":{";
  const std::size_t filtered = static_cast<std::size_t>(
      std::count_if(features.begin(), features.end(),
                    [](const auto &feature) { return feature.filtered_numt; }));
  out << "\"input_reads\":" << reads.size() << ",";
  out << "\"input_alignment_records\":" << input_record_count << ",";
  out << "\"input_molecules\":" << reads.size() << ",";
  out << "\"passed_reads\":" << (features.size() - filtered) << ",";
  out << "\"numt_filtered_reads\":" << filtered << ",";
  if (config.result_schema == ResultSchema::v0_6) {
    std::size_t evidence_eligible = 0U;
    for (std::size_t index = 0U; index < assembly.molecules.size(); ++index) {
      if (assembly.molecules[index].analysis_eligible &&
          index < features.size() && !features[index].filtered_numt) {
        ++evidence_eligible;
      }
    }
    const auto ambiguous = static_cast<std::size_t>(
        std::count_if(assembly.molecules.begin(), assembly.molecules.end(),
                      [](const auto &molecule) { return molecule.ambiguous; }));
    out << "\"evidence_eligible_molecules\":" << evidence_eligible << ","
        << "\"ambiguous_molecules\":" << ambiguous << ",";
  }
  out << "\"numt_threshold\":" << config.numt_threshold << ",";
  out << "\"numt_assessment\":{"
      << "\"mode\":"
      << quoted(has_nuclear_contigs ? "competitive_alignment"
                                    : (alignment_input ? "mt_only_or_unknown"
                                                       : "unaligned_fastq"))
      << ",\"nuclear_contigs_present\":"
      << (has_nuclear_contigs ? "true" : "false")
      << ",\"specificity_assessable\":"
      << (has_nuclear_contigs ? "true" : "false") << "}";
  out << "},";

  if (config.result_schema == ResultSchema::v0_6) {
    const auto observation_page_count =
        evidence_store.observations.size() / config.evidence_page_size +
        (evidence_store.observations.size() % config.evidence_page_size == 0U
             ? 0U
             : 1U);
    out << "\"evidence_encoding\":{"
        << "\"layout\":" << quoted("paged_columnar_molecule_event") << ","
        << "\"scope\":" << quoted("snv_indel_sv_complex_evidence_rc2") << ","
        << "\"observation_storage\":" << quoted("embedded_columnar_pages")
        << ","
        << "\"missing_pair_state\":" << quoted("NOT_CALLABLE") << ","
        << "\"phase_molecule_policy\":" << quoted("evidence_eligible_only")
        << ",\"phase_molecule_reference\":"
        << quoted("molecules[].index") << ","
        << "\"phase_null_model\":"
        << quoted("independent_marginals_within_jointly_callable") << ","
        << "\"observation_limit\":" << config.max_evidence_observations << ","
        << "\"observation_count\":" << evidence_store.observations.size()
        << ",\"observation_page_size\":" << config.evidence_page_size << ","
        << "\"observation_page_count\":" << observation_page_count << ","
        << "\"phase_link_limit\":" << config.max_phase_links << "},";

    out << "\"alignments\":[";
    for (std::size_t index = 0; index < assembly.fragments.size(); ++index) {
      if ((index & 1023U) == 0U) {
        throw_if_cancelled(config);
      }
      if (index != 0U) {
        out << ",";
      }
      const auto &fragment = assembly.fragments[index];
      const ReadRecord *source_record = nullptr;
      if (!assembly.source_alignment_records.empty() &&
          fragment.source_record_index <
              assembly.source_alignment_records.size()) {
        source_record =
            &assembly.source_alignment_records[fragment.source_record_index];
      } else if (fragment.molecule_index.value <
                 assembly.representatives.size()) {
        source_record =
            &assembly.representatives[fragment.molecule_index.value];
      }
      if (source_record == nullptr) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "schema 0.6 alignment source record is unresolved");
      }
      out << "{\"id\":" << quoted(alignment_fragment_id_string(fragment.id))
          << ","
          << "\"source_record_index\":" << fragment.source_record_index << ","
          << "\"molecule_id\":" << quoted(fragment.molecule_id) << ","
          << "\"molecule_index\":" << fragment.molecule_index.value << ","
          << "\"role\":" << quoted(fragment.role) << ","
          << "\"selected_representative\":"
          << (fragment.selected_representative ? "true" : "false") << ","
          << "\"flags\":" << fragment.flags << ","
          << "\"strand\":" << quoted((fragment.flags & 0x10U) == 0U ? "+" : "-")
          << ","
          << "\"mapping_quality\":"
          << static_cast<unsigned int>(fragment.mapping_quality)
          << ",\"reference_name\":" << quoted(fragment.reference_name) << ","
          << "\"reference_start\":" << fragment.reference_start << ","
          << "\"cigar\":" << quoted(fragment.cigar) << ","
          << "\"query_length\":" << source_record->sequence.size() << ","
          << "\"base_qualities_available\":"
          << (!source_record->qualities.empty() ? "true" : "false") << ","
          << "\"aux_tags\":";
      write_string_map(out, source_record->aux_tags);
      out << "}";
    }
    out << "],";

    out << "\"molecules\":[";
    for (std::size_t index = 0; index < assembly.molecules.size(); ++index) {
      if ((index & 1023U) == 0U) {
        throw_if_cancelled(config);
      }
      if (index != 0U) {
        out << ",";
      }
      const auto &molecule = assembly.molecules[index];
      const auto &molecule_callability = callability.molecules[index];
      const auto &feature = features[index];
      out << "{\"id\":" << quoted(molecule.id) << ","
          << "\"index\":" << molecule.index.value << ","
          << "\"identity_policy\":" << quoted(molecule.identity_policy) << ","
          << "\"assembly_status\":" << quoted(molecule.assembly_status) << ","
          << "\"primary_candidate_count\":" << molecule.primary_candidate_count
          << ","
          << "\"ambiguous\":" << (molecule.ambiguous ? "true" : "false") << ","
          << "\"analysis_eligible\":"
          << (molecule.analysis_eligible ? "true" : "false") << ","
          << "\"evidence_eligible\":"
          << (molecule.analysis_eligible && index < features.size() &&
                      !features[index].filtered_numt
                  ? "true"
                  : "false")
          << ","
          << "\"callability_status\":" << quoted(molecule_callability.status)
          << ","
          << "\"callable_bases\":" << molecule_callability.callable_bases << ","
          << "\"callable_fraction\":"
          << (reference.empty()
                  ? 0.0
                  : static_cast<double>(molecule_callability.callable_bases) /
                        static_cast<double>(reference.size()))
          << ","
          << "\"query_length\":" << feature.length << ","
          << "\"mean_base_quality\":" << feature.mean_quality << ","
          << "\"mapping_quality\":"
          << static_cast<unsigned int>(feature.mapping_quality) << ","
          << "\"numt_score\":" << feature.numt_score << ","
          << "\"numt_evidence\":";
      write_string_array(out, feature.numt_evidence);
      out << ",\"cluster_id\":" << feature.cluster_id << ","
          << "\"alternate_event_ids\":";
      write_string_array(out, molecule_alternate_event_ids[index]);
      const auto &evidence_counts = molecule_evidence_counts[index];
      out << ",\"evidence_state_counts\":{\"alternate\":"
          << evidence_counts[0] << ",\"reference\":" << evidence_counts[1]
          << ",\"event_absent\":" << evidence_counts[2]
          << ",\"low_quality\":" << evidence_counts[3]
          << ",\"conflict\":" << evidence_counts[4] << "},"
          << "\"representative_alignment_id\":"
          << quoted(alignment_fragment_id_string(
                 molecule.representative_fragment_id))
          << ","
          << "\"source_qnames\":";
      write_string_array(out, molecule.source_qnames);
      out << ",\"protocol_metadata\":";
      write_string_map(out, molecule.protocol_metadata);
      out << ",\"protocol_flags\":";
      write_string_array(out, molecule.protocol_flags);
      out << ",\"exclusion_reasons\":";
      write_string_array(out, molecule.exclusion_reasons);
      out << ","
          << "\"alignment_ids\":[";
      for (std::size_t fragment_index = 0;
           fragment_index < molecule.fragment_ids.size(); ++fragment_index) {
        if (fragment_index != 0U) {
          out << ",";
        }
        out << quoted(alignment_fragment_id_string(
            molecule.fragment_ids[fragment_index]));
      }
      out << "],\"warnings\":";
      write_string_array(out, molecule.warnings);
      out << "}";
    }
    out << "],";

    out << "\"callability\":[";
    for (std::size_t index = 0U; index < callability.molecules.size();
         ++index) {
      if ((index & 1023U) == 0U) {
        throw_if_cancelled(config);
      }
      if (index != 0U) {
        out << ",";
      }
      const auto &summary = callability.molecules[index];
      out << "{\"molecule_id\":" << quoted(assembly.molecules[index].id) << ","
          << "\"status\":" << quoted(summary.status) << ","
          << "\"known\":" << (summary.known ? "true" : "false") << ","
          << "\"basis\":" << quoted("passing_aligned_reference_bases") << ","
          << "\"callable_bases\":" << summary.callable_bases << ","
          << "\"callable_fraction\":"
          << (reference.empty() ? 0.0
                                : static_cast<double>(summary.callable_bases) /
                                      static_cast<double>(reference.size()))
          << ",\"ranges\":[";
      for (std::size_t range_index = 0U; range_index < summary.ranges.size();
           ++range_index) {
        if (range_index != 0U) {
          out << ",";
        }
        out << "{\"start\":" << summary.ranges[range_index].first
            << ",\"end\":" << summary.ranges[range_index].second << "}";
      }
      out << "],\"alignments\":[";
      for (std::size_t alignment_index = 0U;
           alignment_index < summary.alignments.size(); ++alignment_index) {
        if (alignment_index != 0U) {
          out << ",";
        }
        const auto &alignment = summary.alignments[alignment_index];
        out << "{\"alignment_id\":"
            << quoted(alignment_fragment_id_string(
                   alignment.alignment_fragment_id))
            << ","
            << "\"eligible\":" << (alignment.eligible ? "true" : "false") << ","
            << "\"status\":" << quoted(alignment.status) << ","
            << "\"callable_bases\":" << alignment.callable_bases << ","
            << "\"inserted_query_bases\":" << alignment.inserted_query_bases
            << ","
            << "\"soft_clipped_query_bases\":"
            << alignment.soft_clipped_query_bases << ","
            << "\"reference_exclusion_counts\":{";
        std::size_t exclusion_index = 0U;
        for (const auto &[reason, count] :
             alignment.reference_exclusion_counts) {
          if (exclusion_index++ != 0U) {
            out << ",";
          }
          out << quoted(reason) << ":" << count;
        }
        out << "},\"disrupted_adjacency_anchors\":[";
        for (std::size_t anchor_index = 0U;
             anchor_index < alignment.disrupted_adjacency_anchors.size();
             ++anchor_index) {
          if (anchor_index != 0U) {
            out << ",";
          }
          out << alignment.disrupted_adjacency_anchors[anchor_index];
        }
        out << "],\"ranges\":[";
        for (std::size_t range_index = 0U;
             range_index < alignment.ranges.size(); ++range_index) {
          if (range_index != 0U) {
            out << ",";
          }
          out << "{\"start\":" << alignment.ranges[range_index].first
              << ",\"end\":" << alignment.ranges[range_index].second << "}";
        }
        out << "]}";
      }
      out << "]}";
    }
    out << "],";

    struct SerializedEventCounts {
      std::size_t alternate = 0U;
      std::size_t reference = 0U;
      std::size_t event_absent = 0U;
      std::size_t low_quality = 0U;
      std::size_t conflict = 0U;
    };
    std::vector<SerializedEventCounts> event_counts(
        evidence_store.events.size());
    for (const auto &observation : evidence_store.observations) {
      auto &counts = event_counts[observation.event_index];
      switch (observation.state) {
      case ObservationState::alternate:
        ++counts.alternate;
        break;
      case ObservationState::reference:
        ++counts.reference;
        break;
      case ObservationState::event_absent:
        ++counts.event_absent;
        break;
      case ObservationState::low_quality:
        ++counts.low_quality;
        break;
      case ObservationState::conflict:
        ++counts.conflict;
        break;
      case ObservationState::not_callable:
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 serializer received explicit NOT_CALLABLE");
      }
    }

    out << "\"events\":[";
    for (std::size_t index = 0; index < evidence_store.events.size(); ++index) {
      if ((index & 1023U) == 0U) {
        throw_if_cancelled(config);
      }
      if (index != 0U) {
        out << ",";
      }
      const auto &event = evidence_store.events[index];
      const auto &counts = event_counts[index];
      out << "{\"id\":" << quoted(event.id) << ","
          << "\"index\":" << index << ","
          << "\"type\":" << quoted(event.type) << ",\"start\":";
      if (event.start == 0U) {
        out << "null";
      } else {
        out << event.start;
      }
      out << ",\"end\":";
      if (event.end == 0U) {
        out << "null";
      } else {
        out << event.end;
      }
      out << ",\"length\":" << event.length << ",\"ref\":";
      if (event.reference.empty()) {
        out << "null";
      } else {
        out << quoted(event.reference);
      }
      out << ",\"alt\":";
      if (event.alternate.empty()) {
        out << "null";
      } else {
        out << quoted(event.alternate);
      }
      out << ",\"normalization\":" << quoted(event.normalization) << ","
          << "\"source_projection\":" << quoted(event.source_projection) << ","
          << "\"negative_evidence_rule\":"
          << quoted(event.negative_evidence_rule) << ","
          << "\"assessability\":"
          << quoted(event.absence_assessable ? "REFERENCE_AND_ALTERNATE"
                                             : "ALTERNATE_SUPPORT_ONLY")
          << ",\"component_event_ids\":";
      write_string_array(out, event.component_event_ids);
      out << ",\"supporting_molecule_ids\":";
      write_string_array(out, event.supporting_molecules);
      out << ",\"evidence_counts\":{"
          << "\"alternate\":" << counts.alternate << ","
          << "\"reference\":" << counts.reference << ","
          << "\"event_absent\":" << counts.event_absent << ","
          << "\"callable\":"
          << (counts.alternate + counts.reference + counts.event_absent)
          << ",\"low_quality\":" << counts.low_quality << ","
          << "\"conflict\":" << counts.conflict << "}}";
    }
    out << "],";

    const auto write_observation_column = [&](const std::size_t begin,
                                              const std::size_t end,
                                              const auto &write_value) {
      out << "[";
      for (std::size_t index = begin; index < end; ++index) {
        if (index != begin) {
          out << ",";
        }
        write_value(evidence_store.observations[index]);
      }
      out << "]";
    };
    out << "\"observation_pages\":[";
    for (std::size_t page_index = 0U; page_index < observation_page_count;
         ++page_index) {
      throw_if_cancelled(config);
      if (page_index != 0U) {
        out << ",";
      }
      const auto begin = page_index * config.evidence_page_size;
      const auto end = std::min(evidence_store.observations.size(),
                                begin + config.evidence_page_size);
      out << "{\"index\":" << page_index << ",\"offset\":" << begin
          << ",\"count\":" << (end - begin) << ",\"columns\":{";
      out << "\"molecule_id\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(assembly.molecules[observation.molecule_index.value].id);
      });
      out << ",\"event_id\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(evidence_store.events[observation.event_index].id);
      });
      out << ",\"alignment_id\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(
            alignment_fragment_id_string(observation.alignment_fragment_id));
      });
      out << ",\"state\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(observation_state_name(observation.state));
      });
      out << ",\"observed_allele\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        if (observation.observed_allele.empty()) {
          out << "null";
        } else {
          out << quoted(observation.observed_allele);
        }
      });
      out << ",\"base_quality\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        if (observation.base_quality) {
          out << static_cast<unsigned int>(*observation.base_quality);
        } else {
          out << "null";
        }
      });
      out << ",\"mapping_quality\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << static_cast<unsigned int>(observation.mapping_quality);
      });
      out << ",\"strand\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(std::string(1, observation.strand));
      });
      out << ",\"evidence_source\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(observation.evidence_source);
      });
      out << ",\"read_position\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        if (observation.center_proximity) {
          out << *observation.center_proximity;
        } else {
          out << "null";
        }
      });
      out << "}}";
    }
    out << "],\"phase_links\":[";
    for (std::size_t index = 0U; index < phase_links.size(); ++index) {
      if ((index & 1023U) == 0U) {
        throw_if_cancelled(config);
      }
      if (index != 0U) {
        out << ",";
      }
      const auto &link = phase_links[index];
      const auto &event_a = evidence_store.events[link.event_a_index];
      const auto &event_b = evidence_store.events[link.event_b_index];
      if (link.supporting_molecule_indices.size() != link.both_alternate ||
          link.uncertain_molecule_indices.size() != link.jointly_uncertain ||
          std::ranges::any_of(link.supporting_molecule_indices,
                              [&](const auto molecule_index) {
                                return molecule_index >= assembly.molecules.size();
                              }) ||
          std::ranges::any_of(link.uncertain_molecule_indices,
                              [&](const auto molecule_index) {
                                return molecule_index >= assembly.molecules.size();
                              })) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "phase molecule-index traceability is inconsistent");
      }
      out << "{\"id\":" << quoted("phase:" + event_a.id + "|" + event_b.id)
          << ","
          << "\"event_a_id\":" << quoted(event_a.id) << ","
          << "\"event_b_id\":" << quoted(event_b.id) << ","
          << "\"assessability\":"
          << quoted(link.complete_callability ? "COMPLETE_FOR_BOTH_EVENTS"
                                              : "SUPPORT_CONDITIONED")
          << ",\"jointly_callable\":" << link.jointly_callable << ","
          << "\"jointly_uncertain\":" << link.jointly_uncertain << ","
          << "\"both_alternate\":" << link.both_alternate << ","
          << "\"a_alternate_b_absent\":" << link.a_alternate_b_absent << ","
          << "\"a_absent_b_alternate\":" << link.a_absent_b_alternate << ","
          << "\"neither_alternate\":" << link.neither_alternate << ","
          << "\"co_alternate_fraction\":" << link.co_alternate_fraction << ","
          << "\"co_alternate_ci95_low\":" << link.ci95_low << ","
          << "\"co_alternate_ci95_high\":" << link.ci95_high << ","
          << "\"expected_co_alternate_fraction\":"
          << link.expected_co_alternate_fraction << ","
          << "\"linkage_delta\":" << link.linkage_delta << ","
          << "\"supporting_molecule_indices\":";
      write_integer_array(out, link.supporting_molecule_indices);
      out << ",\"uncertain_molecule_indices\":";
      write_integer_array(out, link.uncertain_molecule_indices);
      out << ",\"qc_flags\":[";
      bool wrote_phase_flag = false;
      if (!link.complete_callability) {
        out << quoted("SUPPORT_CONDITIONED");
        wrote_phase_flag = true;
      }
      if (link.jointly_uncertain != 0U) {
        if (wrote_phase_flag) {
          out << ",";
        }
        out << quoted("UNCERTAIN_COOCCURRENCE_EXCLUDED");
      }
      out << "]}";
    }
    out << "],\"architectures\":[],";
  }

  out << "\"genes\":[";
  const auto genes = default_genes();
  for (std::size_t i = 0; i < genes.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << "{"
        << "\"name\":" << quoted(genes[i].name) << ","
        << "\"start\":" << genes[i].start << ","
        << "\"end\":" << genes[i].end << ","
        << "\"strand\":" << quoted(genes[i].strand) << ","
        << "\"biotype\":" << quoted(genes[i].biotype) << "}";
  }
  out << "],";

  out << "\"coverage\":[";
  for (std::size_t i = 0; i < coverage.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << "{"
        << "\"start\":" << coverage[i].start << ","
        << "\"end\":" << coverage[i].end << ","
        << "\"depth\":" << coverage[i].depth << "}";
  }
  out << "],";

  const auto bins_gt20 = static_cast<std::size_t>(
      std::count_if(coverage.begin(), coverage.end(),
                    [](const auto &bin) { return bin.depth > 20; }));
  const double pct_bins_gt20 = coverage.empty()
                                   ? 0.0
                                   : (static_cast<double>(bins_gt20) * 100.0) /
                                         static_cast<double>(coverage.size());
  std::map<std::uint8_t, std::size_t> mapq_histogram;
  for (std::size_t i = 0; i < reads.size(); ++i) {
    if (i < features.size() && !features[i].filtered_numt &&
        (reads[i].flags & 0x4U) == 0U && !reads[i].reference_name.empty() &&
        reads[i].reference_name != "*") {
      if (!looks_like_nuclear_contig(reads[i].reference_name)) {
        ++mapq_histogram[reads[i].mapping_quality];
      }
    }
  }
  out << "\"coverage_metrics\":{";
  out << "\"mean_depth\":" << coverage_result.mean_depth << ",";
  out << "\"pct_sites_gt20x\":" << coverage_result.pct_sites_gt20x << ",";
  out << "\"pct_bins_gt20x\":" << pct_bins_gt20 << ",";
  out << "\"max_depth\":" << coverage_result.max_depth << ",";
  out << "\"mapping_quality_histogram\":[";
  std::size_t mapq_index = 0;
  for (const auto &[mapq, count] : mapq_histogram) {
    if (mapq_index++ != 0) {
      out << ",";
    }
    out << "{"
        << "\"mapq\":" << static_cast<int>(mapq) << ","
        << "\"count\":" << count << "}";
  }
  out << "]";
  out << "},";

  out << "\"svs\":[";
  std::size_t sv_index = 0;
  for (const auto &[_, sv] : svs) {
    if (sv_index++ != 0) {
      out << ",";
    }
    out << "{"
        << "\"id\":" << quoted(sv.id) << ",";
    if (config.result_schema == ResultSchema::v0_6) {
      out << "\"event_id\":" << quoted("sv:" + sv.id) << ",";
    }
    out << "\"type\":" << quoted(sv.type) << ","
        << "\"start\":" << sv.start << ","
        << "\"end\":" << sv.end << ","
        << "\"length\":" << sv.length << ","
        << "\"known_event\":" << (sv.known_event ? "true" : "false") << ","
        << "\"evidence_source\":"
        << quoted(sv.evidence_sources.size() == 1U ? sv.evidence_sources.front()
                                                   : "combined")
        << ",\"evidence_sources\":";
    write_string_array(out, sv.evidence_sources);
    out << ",\"segment_count\":" << sv.segment_count << ",";
    if (!sv.orientations.empty()) {
      out << "\"orientation\":"
          << quoted(sv.orientations.size() == 1U ? sv.orientations.front()
                                                 : "mixed")
          << ",\"orientations\":";
      write_string_array(out, sv.orientations);
      out << ",";
    }
    out << "\"supporting_reads\":[";
    for (std::size_t i = 0; i < sv.supporting_reads.size(); ++i) {
      if (i != 0) {
        out << ",";
      }
      out << quoted(sv.supporting_reads[i]);
    }
    out << "]}";
  }
  out << "],";

  out << "\"complex_events\":[";
  std::size_t complex_event_index = 0;
  for (const auto &[_, event] : complex_events) {
    if (complex_event_index++ != 0U) {
      out << ",";
    }
    out << "{"
        << "\"id\":" << quoted(event.id) << ",";
    if (config.result_schema == ResultSchema::v0_6) {
      out << "\"event_id\":" << quoted(event.id) << ",";
    }
    out << "\"junction_count\":" << event.junction_ids.size() << ","
        << "\"segment_count\":" << event.segment_count << ","
        << "\"canonicalization\":" << quoted("strand_invariant_path") << ","
        << "\"junction_ids\":";
    write_string_array(out, event.junction_ids);
    out << ",\"junction_orientations\":";
    write_string_array(out, event.junction_orientations);
    out << ",\"supporting_reads\":";
    write_string_array(out, event.supporting_reads);
    out << "}";
  }
  out << "],";

  out << "\"variants\":[";
  std::size_t variant_index = 0;
  if (config.result_schema == ResultSchema::v0_6) {
    const auto unified_summaries =
        summarize_unified_variant_evidence(evidence_store);
    for (std::size_t event_index = 0U;
         event_index < evidence_store.events.size(); ++event_index) {
      const auto &event = evidence_store.events[event_index];
      if (!is_unified_variant_event(event)) {
        continue;
      }
      if (variant_index++ != 0U) {
        out << ",";
      }
      const auto &summary = unified_summaries[event_index];
      const auto callable_depth = summary.alternate.count +
                                  summary.reference.count +
                                  summary.other.count;
      const double heteroplasmy =
          callable_depth == 0U
              ? 0.0
              : static_cast<double>(summary.alternate.count) /
                    static_cast<double>(callable_depth);
      const auto [ci95_low, ci95_high] =
          wilson_interval(summary.alternate.count, callable_depth);
      const auto homopolymer_run =
          circular_reference_homopolymer_run(reference, event.start);
      const auto qc_flags = variant_qc_flags(
          event, summary, homopolymer_run, has_nuclear_contigs);
      const bool circular_origin_event =
          event.type == "SMALL_DELETION" &&
          (event.start == 1U || event.end < event.start);
      const bool vcf_representable = !circular_origin_event;
      const auto vcf_position =
          event.type == "SMALL_DELETION" && event.start > 1U
              ? event.start - 1U
              : event.start;
      const SnpCall *snp = nullptr;
      if (event.type == "SNV" && event.reference.size() == 1U &&
          event.alternate.size() == 1U) {
        const auto found = variants.find(snp_key(
            event.start, event.reference.front(), event.alternate.front()));
        if (found == variants.end()) {
          throw AnalysisError(
              AnalysisErrorCode::internal_error,
              "schema 0.6 SNV event has no aggregate annotation projection");
        }
        snp = &found->second.call;
      }

      out << "{\"event_id\":" << quoted(event.id) << ","
          << "\"type\":" << quoted(event.type) << ","
          << "\"position\":" << event.start << ","
          << "\"start\":" << event.start << ","
          << "\"end\":" << event.end << ","
          << "\"length\":" << event.length << ","
          << "\"ref\":" << quoted(event.reference) << ","
          << "\"alt\":" << quoted(event.alternate) << ","
          << "\"normalization\":" << quoted(event.normalization) << ","
          << "\"negative_evidence_rule\":"
          << quoted(event.negative_evidence_rule) << ","
          << "\"assessability\":"
          << quoted(event.absence_assessable ? "REFERENCE_AND_ALTERNATE"
                                             : "ALTERNATE_SUPPORT_ONLY")
          << ",\"vcf_position\":" << vcf_position << ","
          << "\"vcf_representable\":"
          << (vcf_representable ? "true" : "false") << ","
          << "\"alt_depth\":" << summary.alternate.count << ","
          << "\"ref_depth\":" << summary.reference.count << ","
          << "\"other_depth\":" << summary.other.count << ","
          << "\"event_absent_depth\":" << summary.event_absent << ","
          << "\"low_quality_depth\":" << summary.low_quality << ","
          << "\"conflict_depth\":" << summary.conflict << ","
          << "\"callable_depth\":" << callable_depth << ","
          << "\"heteroplasmy\":" << heteroplasmy << ","
          << "\"ci95_low\":" << ci95_low << ","
          << "\"ci95_high\":" << ci95_high << ","
          << "\"filter_status\":" << quoted("NOT_CALIBRATED") << ","
          << "\"qc_flags\":";
      write_string_array(out, qc_flags);
      out << ",\"numt_assessability\":"
          << quoted(has_nuclear_contigs ? "ASSESSABLE"
                                       : "NOT_ASSESSABLE")
          << ",\"multi_allelic\":"
          << (summary.multi_allelic ? "true" : "false") << ","
          << "\"homopolymer_context\":{\"reference_base\":";
      if (event.start == 0U || event.start > reference.size()) {
        out << "null";
      } else {
        out << quoted(std::string(1U, reference[event.start - 1U]));
      }
      out << ",\"run_length\":" << homopolymer_run << "},"
          << "\"molecule_support\":{\"alternate\":"
          << summary.alternate.count << ",\"reference\":"
          << summary.reference.count << ",\"other\":"
          << summary.other.count << ",\"callable\":" << callable_depth
          << "},\"strand_support\":{\"alt_forward\":"
          << summary.alternate.strand_depths[0] << ",\"alt_reverse\":"
          << summary.alternate.strand_depths[1] << ",\"ref_forward\":"
          << summary.reference.strand_depths[0] << ",\"ref_reverse\":"
          << summary.reference.strand_depths[1] << ",\"other_forward\":"
          << summary.other.strand_depths[0] << ",\"other_reverse\":"
          << summary.other.strand_depths[1]
          << "},\"strand_bias_delta\":";
      if (summary.alternate.count == 0U || summary.reference.count == 0U) {
        out << "null";
      } else {
        const double alternate_forward_fraction =
            static_cast<double>(summary.alternate.strand_depths[0]) /
            static_cast<double>(summary.alternate.count);
        const double reference_forward_fraction =
            static_cast<double>(summary.reference.strand_depths[0]) /
            static_cast<double>(summary.reference.count);
        out << std::abs(alternate_forward_fraction -
                        reference_forward_fraction);
      }
      out << ",\"allele_quality\":{\"alternate\":";
      write_base_quality_summary(out, summary.alternate);
      out << ",\"reference\":";
      write_base_quality_summary(out, summary.reference);
      out << ",\"other\":";
      write_base_quality_summary(out, summary.other);
      out << "},\"mapping_quality\":{\"alternate\":";
      write_mapping_quality_summary(out, summary.alternate);
      out << ",\"reference\":";
      write_mapping_quality_summary(out, summary.reference);
      out << ",\"other\":";
      write_mapping_quality_summary(out, summary.other);
      out << "},\"read_position\":{\"definition\":"
          << quoted("normalized_center_proximity")
          << ",\"alternate_mean\":";
      write_nullable_mean(out, summary.alternate.read_position_sum,
                          summary.alternate.read_position_count);
      out << ",\"reference_mean\":";
      write_nullable_mean(out, summary.reference.read_position_sum,
                          summary.reference.read_position_count);
      out << ",\"other_mean\":";
      write_nullable_mean(out, summary.other.read_position_sum,
                          summary.other.read_position_count);
      out << ",\"bias_delta\":";
      if (summary.alternate.read_position_count == 0U ||
          summary.reference.read_position_count == 0U) {
        out << "null";
      } else {
        const double alternate_mean =
            summary.alternate.read_position_sum /
            static_cast<double>(summary.alternate.read_position_count);
        const double reference_mean =
            summary.reference.read_position_sum /
            static_cast<double>(summary.reference.read_position_count);
        out << std::abs(alternate_mean - reference_mean);
      }
      out << "},\"supporting_molecule_ids\":";
      write_string_array(out, event.supporting_molecules);
      out << ",\"supporting_reads\":";
      write_string_array(out, event.supporting_molecules);
      if (snp != nullptr && !snp->gene.empty()) {
        out << ",\"gene\":" << quoted(snp->gene);
      }
      if (snp != nullptr && !snp->consequence.empty()) {
        out << ",\"consequence\":" << quoted(snp->consequence);
      }
      if (snp != nullptr && !snp->protein.empty()) {
        out << ",\"protein\":" << quoted(snp->protein);
      }
      if (snp != nullptr && !snp->residue.empty()) {
        out << ",\"residue\":" << quoted(snp->residue);
      }
      if (snp != nullptr && has_clinical_annotation(*snp)) {
        out << ",\"annotation\":";
        write_clinical_annotation(out, *snp);
      }
      if (snp != nullptr && !snp->structure_id.empty()) {
        out << ",\"structure\":{\"structure_id\":"
            << quoted(snp->structure_id);
        if (!snp->structure_chain.empty()) {
          out << ",\"chain\":" << quoted(snp->structure_chain);
        }
        if (snp->structure_residue != 0) {
          out << ",\"residue_index\":" << snp->structure_residue;
        }
        if (!snp->structure_complex.empty()) {
          out << ",\"complex\":" << quoted(snp->structure_complex);
        }
        out << "}";
      }
      out << "}";
    }
  } else {
    for (const auto &[_, aggregate] : variants) {
      if (variant_index++ != 0) {
        out << ",";
      }
      const auto &snp = aggregate.call;
      out << "{\"position\":" << snp.position << ","
        << "\"ref\":" << quoted(std::string(1, snp.reference)) << ","
        << "\"alt\":" << quoted(std::string(1, snp.alternate)) << ","
        << "\"alt_depth\":" << aggregate.alternate_depth << ","
        << "\"ref_depth\":" << aggregate.reference_depth << ","
        << "\"other_depth\":" << aggregate.other_depth << ","
        << "\"callable_depth\":" << aggregate.callable_depth << ","
        << "\"heteroplasmy\":" << aggregate.heteroplasmy << ","
        << "\"ci95_low\":" << aggregate.ci95_low << ","
        << "\"ci95_high\":" << aggregate.ci95_high << ","
        << "\"molecule_support\":{"
        << "\"alternate\":" << aggregate.alternate_depth << ","
        << "\"reference\":" << aggregate.reference_depth << ","
        << "\"other\":" << aggregate.other_depth << ","
        << "\"callable\":" << aggregate.callable_depth << "},"
        << "\"strand_support\":{"
        << "\"alt_forward\":" << aggregate.alternate_strand_depths[0] << ","
        << "\"alt_reverse\":" << aggregate.alternate_strand_depths[1] << ","
        << "\"ref_forward\":" << aggregate.reference_strand_depths[0] << ","
        << "\"ref_reverse\":" << aggregate.reference_strand_depths[1] << ","
        << "\"other_forward\":" << aggregate.other_strand_depths[0] << ","
        << "\"other_reverse\":" << aggregate.other_strand_depths[1] << "},"
        << "\"strand_bias_delta\":";
    if (aggregate.alternate_depth == 0U || aggregate.reference_depth == 0U) {
      out << "null";
    } else {
      const double alternate_forward_fraction =
          static_cast<double>(aggregate.alternate_strand_depths[0]) /
          static_cast<double>(aggregate.alternate_depth);
      const double reference_forward_fraction =
          static_cast<double>(aggregate.reference_strand_depths[0]) /
          static_cast<double>(aggregate.reference_depth);
      out << std::abs(alternate_forward_fraction - reference_forward_fraction);
    }
    out << ",\"allele_quality\":{\"alternate\":";
    write_quality_summary(
        out, aggregate.alternate_depth, aggregate.alternate_quality_sum,
        aggregate.alternate_quality_min, aggregate.alternate_quality_max);
    out << ",\"reference\":";
    write_quality_summary(
        out, aggregate.reference_depth, aggregate.reference_quality_sum,
        aggregate.reference_quality_min, aggregate.reference_quality_max);
    out << ",\"other\":";
    write_quality_summary(
        out, aggregate.other_depth, aggregate.other_quality_sum,
        aggregate.other_quality_min, aggregate.other_quality_max);
    out << "},\"read_position\":{"
        << "\"definition\":" << quoted("normalized_center_proximity") << ","
        << "\"alternate_mean\":";
    write_nullable_mean(out, aggregate.alternate_read_position_sum,
                        aggregate.alternate_depth);
    out << ",\"reference_mean\":";
    write_nullable_mean(out, aggregate.reference_read_position_sum,
                        aggregate.reference_depth);
    out << ",\"other_mean\":";
    write_nullable_mean(out, aggregate.other_read_position_sum,
                        aggregate.other_depth);
    out << ",\"bias_delta\":";
    if (aggregate.alternate_depth == 0U || aggregate.reference_depth == 0U) {
      out << "null";
    } else {
      const double alternate_mean =
          aggregate.alternate_read_position_sum /
          static_cast<double>(aggregate.alternate_depth);
      const double reference_mean =
          aggregate.reference_read_position_sum /
          static_cast<double>(aggregate.reference_depth);
      out << std::abs(alternate_mean - reference_mean);
    }
    out << "},\"supporting_reads\":";
    write_string_array(out, aggregate.supporting_reads);
    if (!snp.gene.empty()) {
      out << ",\"gene\":" << quoted(snp.gene);
    }
    if (!snp.consequence.empty()) {
      out << ",\"consequence\":" << quoted(snp.consequence);
    }
    if (!snp.protein.empty()) {
      out << ",\"protein\":" << quoted(snp.protein);
    }
    if (!snp.residue.empty()) {
      out << ",\"residue\":" << quoted(snp.residue);
    }
    if (has_clinical_annotation(snp)) {
      out << ",\"annotation\":";
      write_clinical_annotation(out, snp);
    }
    if (!snp.structure_id.empty()) {
      out << ",\"structure\":{"
          << "\"structure_id\":" << quoted(snp.structure_id);
      if (!snp.structure_chain.empty()) {
        out << ",\"chain\":" << quoted(snp.structure_chain);
      }
      if (snp.structure_residue != 0) {
        out << ",\"residue_index\":" << snp.structure_residue;
      }
      if (!snp.structure_complex.empty()) {
        out << ",\"complex\":" << quoted(snp.structure_complex);
      }
      out << "}";
    }
      out << "}";
    }
  }
  out << "],";

  out << "\"clusters\":[";
  std::size_t cluster_index = 0;
  for (const auto &[cluster_id, members] : clusters) {
    if (cluster_index++ != 0) {
      out << ",";
    }
    std::map<std::string, std::size_t> sv_counts;
    std::map<std::string, std::size_t> complex_event_counts;
    std::vector<std::string> read_ids;
    for (const auto *member : members) {
      read_ids.push_back(member->id);
      for (const auto &sv_id : member->sv_ids) {
        ++sv_counts[sv_id];
      }
      for (const auto &complex_event_id : member->complex_event_ids) {
        ++complex_event_counts[complex_event_id];
      }
    }
    const auto assignment_it = haplogroups.find(cluster_id);
    const ClusterHaplogroupAssignment *assignment =
        assignment_it == haplogroups.end() ? nullptr : &assignment_it->second;
    out << "{"
        << "\"id\":" << cluster_id << ","
        << "\"label\":"
        << quoted(cluster_id < 0 ? "Outliers"
                                 : "Cluster " + std::to_string(cluster_id + 1))
        << ","
        << "\"haplogroup\":"
        << quoted(assignment == nullptr ? "unassigned" : assignment->best)
        << ","
        << "\"haplogroup_assignment\":{";
    out << "\"resource\":" << quoted("phylotree-rcrs@17.3") << ","
        << "\"quality\":" << (assignment == nullptr ? 0.0 : assignment->quality)
        << ","
        << "\"contamination_warning\":"
        << (assignment != nullptr && assignment->contamination_warning
                ? "true"
                : "false")
        << ","
        << "\"observed_markers\":";
    if (assignment == nullptr) {
      out << "[]";
    } else {
      write_string_array(out, assignment->observed_markers);
    }
    out << ",\"callable_ranges\":[";
    if (assignment != nullptr) {
      for (std::size_t range_index = 0;
           range_index < assignment->callable_ranges.size(); ++range_index) {
        if (range_index != 0U) {
          out << ",";
        }
        const auto &[start, end] = assignment->callable_ranges[range_index];
        out << "{\"start\":" << start << ",\"end\":" << end << "}";
      }
    }
    out << "],\"candidates\":[";
    if (assignment != nullptr) {
      for (std::size_t candidate_index = 0;
           candidate_index < assignment->candidates.size(); ++candidate_index) {
        if (candidate_index != 0) {
          out << ",";
        }
        const auto &candidate = assignment->candidates[candidate_index];
        out << "{\"name\":" << quoted(candidate.name)
            << ",\"score\":" << candidate.score << ",\"matched\":";
        write_string_array(out, candidate.matched);
        out << ",\"missing\":";
        write_string_array(out, candidate.missing);
        out << ",\"extra\":";
        write_string_array(out, candidate.extra);
        out << "}";
      }
    }
    out << "]},"
        << "\"size\":" << members.size() << ","
        << "\"outlier\":" << (cluster_id < 0 ? "true" : "false") << ","
        << "\"consensus_haplotype\":"
        << quoted(cluster_id < 0
                      ? "noise"
                      : "feature-consensus-" + std::to_string(cluster_id))
        << ","
        << "\"reads\":[";
    for (std::size_t i = 0; i < read_ids.size(); ++i) {
      if (i != 0) {
        out << ",";
      }
      out << quoted(read_ids[i]);
    }
    out << "],\"sv_signature\":[";
    std::size_t signature_index = 0;
    for (const auto &[sv_id, count] : sv_counts) {
      if (signature_index++ != 0) {
        out << ",";
      }
      out << "{"
          << "\"sv_id\":" << quoted(sv_id) << ","
          << "\"support\":" << count << "}";
    }
    out << "],\"complex_event_signature\":[";
    std::size_t complex_signature_index = 0;
    for (const auto &[event_id, count] : complex_event_counts) {
      if (complex_signature_index++ != 0U) {
        out << ",";
      }
      out << "{"
          << "\"event_id\":" << quoted(event_id) << ","
          << "\"support\":" << count << "}";
    }
    out << "]}";
  }
  out << "],";

  out << "\"reads\":[";
  for (std::size_t i = 0; i < features.size(); ++i) {
    if ((i & 1023U) == 0U) {
      throw_if_cancelled(config);
    }
    if (i != 0) {
      out << ",";
    }
    const auto &feature = features[i];
    out << "{"
        << "\"id\":" << quoted(feature.id) << ","
        << "\"length\":" << feature.length << ","
        << "\"mean_quality\":" << feature.mean_quality << ","
        << "\"numt_score\":" << feature.numt_score << ","
        << "\"filtered_numt\":" << (feature.filtered_numt ? "true" : "false")
        << ","
        << "\"numt_evidence\":";
    write_string_array(out, feature.numt_evidence);
    out << ","
        << "\"cluster_id\":" << feature.cluster_id << ","
        << "\"outlier\":" << (feature.outlier ? "true" : "false") << ","
        << "\"mapping_quality\":" << static_cast<int>(feature.mapping_quality)
        << ","
        << "\"flags\":" << feature.flags << ",";
    if (!feature.reference_name.empty()) {
      out << "\"reference_name\":" << quoted(feature.reference_name) << ",";
    }
    if (!feature.aux_tags.empty()) {
      out << "\"aux_tags\":";
      write_string_map(out, feature.aux_tags);
      out << ",";
    }
    out << "\"snps\":[";
    for (std::size_t snp_index = 0; snp_index < feature.snps.size();
         ++snp_index) {
      if (snp_index != 0) {
        out << ",";
      }
      const auto &snp = feature.snps[snp_index];
      out << "{"
          << "\"position\":" << snp.position << ","
          << "\"ref\":" << quoted(std::string(1, snp.reference)) << ","
          << "\"alt\":" << quoted(std::string(1, snp.alternate));
      out << "}";
    }
    out << "],\"haplogroup_markers\":[";
    for (std::size_t marker_index = 0;
         marker_index < feature.haplogroup_markers.size(); ++marker_index) {
      if (marker_index != 0U) {
        out << ",";
      }
      out << quoted(feature.haplogroup_markers[marker_index].encoded);
    }
    out << "],\"haplogroup_range_known\":"
        << (feature.haplogroup_range_known ? "true" : "false")
        << ",\"haplogroup_callable_ranges\":[";
    for (std::size_t range_index = 0;
         range_index < feature.haplogroup_ranges.size(); ++range_index) {
      if (range_index != 0U) {
        out << ",";
      }
      const auto &[start, end] = feature.haplogroup_ranges[range_index];
      out << "{\"start\":" << start << ",\"end\":" << end << "}";
    }
    out << "],\"sv_ids\":[";
    for (std::size_t sv_id_index = 0; sv_id_index < feature.sv_ids.size();
         ++sv_id_index) {
      if (sv_id_index != 0) {
        out << ",";
      }
      out << quoted(feature.sv_ids[sv_id_index]);
    }
    out << "],\"complex_event_ids\":";
    write_string_array(out, feature.complex_event_ids);
    out << "}";
  }
  out << "]";
  out << "}";
  return std::move(out).take();
}

using AnalysisClock = std::chrono::steady_clock;

[[nodiscard]] std::uint64_t
elapsed_microseconds(const AnalysisClock::time_point start) {
  return static_cast<std::uint64_t>(
      std::chrono::duration_cast<std::chrono::microseconds>(
          AnalysisClock::now() - start)
          .count());
}

[[nodiscard]] std::string analyze_impl(const std::string &input_path,
                                       const std::string &reference_path,
                                       const AnalysisConfig &config,
                                       AnalysisPhaseTimings *timings) {
  if (timings != nullptr) {
    *timings = {};
  }
  const auto total_start = AnalysisClock::now();
  validate_config(config);
  throw_if_cancelled(config);

  const auto reference_start = AnalysisClock::now();
  const auto reference = load_reference(reference_path);
  if (timings != nullptr) {
    timings->reference_load_us = elapsed_microseconds(reference_start);
  }
  throw_if_cancelled(config);

  const auto input_start = AnalysisClock::now();
  auto input = read_input(input_path, config);
  if (input.reads.empty()) {
    throw AnalysisError(AnalysisErrorCode::input_empty,
                        "input contains no read records: " + input_path);
  }
  const bool alignment_input = input.alignment_input;
  const bool has_nuclear_contigs = input.has_nuclear_contigs;
  const std::size_t input_record_count = input.reads.size();
  auto assembly = assemble_molecules(input, config, reference.size());
  const auto &reads = assembly.representatives;
  if (timings != nullptr) {
    timings->input_ingest_us = elapsed_microseconds(input_start);
  }
  throw_if_cancelled(config);

  std::vector<ReadFeature> features;
  features.reserve(reads.size());
  std::map<std::string, SvCall> svs;
  std::map<std::string, ComplexSvCall> complex_events;
  const auto resource_start = AnalysisClock::now();
  const auto clinical_annotations = load_clinical_annotations();
  static const auto base_resources = load_resource_manifest();
  auto resources = base_resources;
  apply_runtime_resource_overrides(resources);
  static const auto haplogroup_weights = load_phylo_weights();
  static const auto haplogroup_alignment_rules = load_phylo_alignment_rules();
  static const auto haplogroup_definitions = [] {
    auto definitions = load_haplogroups();
    for (auto &definition : definitions) {
      for (auto &mutation : definition.mutations) {
        mutation.weight = mutation_weight(haplogroup_weights, mutation.encoded);
      }
    }
    return definitions;
  }();
  if (timings != nullptr) {
    timings->resource_load_us = elapsed_microseconds(resource_start);
  }

  const auto extraction_start = AnalysisClock::now();
  auto extraction_results =
      extract_features_parallel(reads, reference, config, clinical_annotations);
  if (timings != nullptr) {
    timings->feature_extraction_us = elapsed_microseconds(extraction_start);
  }

  const auto event_merge_start = AnalysisClock::now();
  for (auto &result : extraction_results) {
    throw_if_cancelled(config);
    for (auto &sv : result.svs) {
      merge_sv(svs, std::move(sv));
    }
    for (auto &complex_event : result.complex_events) {
      merge_complex_sv(complex_events, std::move(complex_event));
    }
    features.push_back(std::move(result.feature));
  }
  if (timings != nullptr) {
    timings->event_merge_us = elapsed_microseconds(event_merge_start);
  }

  std::string_view clustering_backend =
      "Read clusters use DBSCAN over SNP/SV feature tokens";
  const auto clustering_start = AnalysisClock::now();
#ifdef MITO_HAS_HDBSCAN_CPP
  if (assign_hdbscan_clusters(features, config.min_cluster_size, config)) {
    clustering_backend =
        "Read clusters use HDBSCAN-C++ over dense SNP/SV feature vectors";
  } else {
    assign_dbscan_clusters(features, config.cluster_epsilon,
                           config.min_cluster_size, config);
    clustering_backend =
        "Read clusters use DBSCAN fallback after HDBSCAN-C++ adapter failure";
  }
#else
  assign_dbscan_clusters(features, config.cluster_epsilon,
                         config.min_cluster_size, config);
#endif
  if (timings != nullptr) {
    timings->clustering_us = elapsed_microseconds(clustering_start);
  }

  throw_if_cancelled(config);
  const auto aggregation_start = AnalysisClock::now();
  const auto coverage =
      compute_coverage(reads, features, config, reference.size());
  auto evidence = aggregate_snps(assembly, features, config, reference.size(),
                                 reference, clinical_annotations);
  CallabilityResult callability;
  std::vector<PhaseLink> phase_links;
  if (config.result_schema == ResultSchema::v0_6) {
    callability = build_callability(assembly, features, reference, config);
    append_small_indel_evidence(evidence.store, assembly, features, callability,
                                reference, config);
    append_structural_evidence(evidence.store, assembly, features, callability,
                               svs, complex_events, reference.size(), config);
    validate_unified_evidence_graph(assembly, callability, evidence.store,
                                    reference.size());
    phase_links = build_phase_links(evidence.store, assembly, config);
  }
  if (timings != nullptr) {
    timings->evidence_aggregation_us = elapsed_microseconds(aggregation_start);
  }

  const auto haplogroup_start = AnalysisClock::now();
  const auto haplogroups =
      assign_haplogroups(features, haplogroup_definitions, haplogroup_weights,
                         haplogroup_alignment_rules, config);
  if (timings != nullptr) {
    timings->haplogroup_assignment_us = elapsed_microseconds(haplogroup_start);
  }

  const auto serialization_start = AnalysisClock::now();
  auto json = render_json(
      input_path, reference_path, config, clustering_backend, alignment_input,
      has_nuclear_contigs, input_record_count, reference, assembly, reads,
      features, svs, complex_events, coverage, evidence.variants,
      evidence.store, callability, phase_links, haplogroups, resources);
  if (timings != nullptr) {
    timings->serialization_us = elapsed_microseconds(serialization_start);
    timings->total_us = elapsed_microseconds(total_start);
  }
  return json;
}

} // namespace

std::string AnalysisEngine::analyze(const std::string &input_path,
                                    const std::string &reference_path,
                                    const AnalysisConfig &config) const {
  return analyze_impl(input_path, reference_path, config, nullptr);
}

ProfiledAnalysis
AnalysisEngine::analyze_profiled(const std::string &input_path,
                                 const std::string &reference_path,
                                 const AnalysisConfig &config) const {
  ProfiledAnalysis result;
  result.json =
      analyze_impl(input_path, reference_path, config, &result.timings);
  return result;
}

} // namespace mito
