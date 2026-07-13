#include "mito/analysis_engine.hpp"

#include "mito/version.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <charconv>
#include <cmath>
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
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <tuple>
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

std::string_view analysis_error_code_name(const AnalysisErrorCode code) noexcept {
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
};

struct InputRecords {
  std::vector<ReadRecord> reads;
  bool alignment_input = false;
  bool has_nuclear_contigs = false;
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
  std::vector<std::string> sv_ids;
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
  std::vector<std::pair<std::size_t, std::string>> mutations;
  double expected_weight = 0.0;
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
  while (first < value.size() && std::isspace(static_cast<unsigned char>(value[first])) != 0) {
    ++first;
  }
  std::size_t last = value.size();
  while (last > first && std::isspace(static_cast<unsigned char>(value[last - 1])) != 0) {
    --last;
  }
  return std::string(value.substr(first, last - first));
}

void strip_trailing_carriage_return(std::string& value) {
  if (!value.empty() && value.back() == '\r') {
    value.pop_back();
  }
}

[[nodiscard]] std::vector<std::string> split_semicolon(std::string_view line) {
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (start <= line.size()) {
    const std::size_t next = line.find(';', start);
    const auto token =
        trim_copy(line.substr(start, next == std::string_view::npos ? line.size() - start : next - start));
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

[[nodiscard]] std::string parse_reference_fasta(std::istream& in, const std::string& source) {
  std::string reference;
  std::string line;
  std::size_t header_count = 0;
  while (std::getline(in, line)) {
    if (!line.empty() && line.front() == '>') {
      ++header_count;
      if (header_count > 1) {
        throw AnalysisError(AnalysisErrorCode::reference_invalid,
                            "reference FASTA must contain exactly one sequence: " + source);
      }
      continue;
    }
    for (const char base : line) {
      if (std::isspace(static_cast<unsigned char>(base)) != 0) {
        continue;
      }
      const char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
      if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T' || upper == 'N') {
        reference.push_back(upper);
      } else {
        throw AnalysisError(AnalysisErrorCode::reference_invalid,
                            "reference FASTA contains an unsupported base: " + source);
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

[[nodiscard]] std::string clinical_annotations_path() {
  if (const char* override_path = std::getenv("MITO_CLINICAL_ANNOTATIONS")) {
    if (override_path[0] != '\0') {
      return override_path;
    }
  }
#ifdef MITO_CLINICAL_ANNOTATIONS_PATH
  return MITO_CLINICAL_ANNOTATIONS_PATH;
#else
  return "core/data/clinical_annotations.tsv";
#endif
}

[[nodiscard]] std::string phylotree_path() {
  if (const char* override_path = std::getenv("MITO_PHYLOTREE")) {
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
  if (const char* override_path = std::getenv("MITO_PHYLOTREE_WEIGHTS")) {
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
                        "could not open resource manifest: " + resource_manifest_path());
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
    records.push_back({fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6]});
  }
  if (records.empty()) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "resource manifest contains no resources");
  }
  return records;
}

void throw_if_cancelled(const AnalysisConfig& config) {
  if (config.should_cancel && config.should_cancel()) {
    throw AnalysisError(AnalysisErrorCode::analysis_cancelled, "analysis cancelled");
  }
}

void validate_config(const AnalysisConfig& config) {
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
}

[[nodiscard]] std::string load_reference(const std::string& reference_path) {
  const std::string effective_path = reference_path.empty() ? bundled_rcrs_path() : reference_path;
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
      {"MT-RNR1", 648, 1601, "+", "rRNA"},    {"MT-RNR2", 1671, 3229, "+", "rRNA"},
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
  const auto* begin = value.data();
  const auto* end = value.data() + value.size();
  const auto result = std::from_chars(begin, end, parsed);
  if (result.ec != std::errc{} || result.ptr != end) {
    return std::nullopt;
  }
  return parsed;
}

[[nodiscard]] std::vector<CigarOperation> parse_cigar(std::string_view cigar,
                                                      std::string_view read_id) {
  std::vector<CigarOperation> operations;
  if (cigar.empty() || cigar == "*") {
    return operations;
  }

  std::size_t number_start = 0;
  for (std::size_t i = 0; i < cigar.size(); ++i) {
    if (std::isdigit(static_cast<unsigned char>(cigar[i])) != 0) {
      continue;
    }

    const auto length = parse_size(cigar.substr(number_start, i - number_start));
    constexpr std::string_view valid_operations = "MIDNSHP=X";
    if (!length || *length == 0 || valid_operations.find(cigar[i]) == std::string_view::npos) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid CIGAR for read '" + std::string(read_id) + "': " +
                              std::string(cigar));
    }
    operations.push_back({*length, cigar[i]});
    number_start = i + 1;
  }

  if (number_start != cigar.size() || operations.empty()) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "invalid CIGAR for read '" + std::string(read_id) + "': " +
                            std::string(cigar));
  }
  return operations;
}

[[nodiscard]] std::map<std::string, std::string> parse_sam_aux_tags(const std::vector<std::string>& fields) {
  std::map<std::string, std::string> tags;
  for (std::size_t i = 11; i < fields.size(); ++i) {
    if (fields[i].size() < 5 || fields[i][2] != ':' || fields[i][4] != ':') {
      continue;
    }
    tags.emplace(fields[i].substr(0, 2), fields[i].substr(5));
  }
  return tags;
}

[[nodiscard]] std::string snp_key(std::size_t position, char reference, char alternate) {
  std::string key = std::to_string(position);
  key.push_back(':');
  key.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(reference))));
  key.push_back(':');
  key.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(alternate))));
  return key;
}

[[nodiscard]] std::map<std::string, SnpCall> load_clinical_annotations() {
  std::map<std::string, SnpCall> annotations;
  std::ifstream input(clinical_annotations_path());
  if (!input) {
    return annotations;
  }

  std::string line;
  bool header = true;
  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }
    if (header) {
      header = false;
      continue;
    }

    const auto fields = split_tab(line);
    if (fields.size() < 11) {
      continue;
    }
    const auto position = parse_size(fields[0]);
    if (!position || fields[1].empty() || fields[2].empty()) {
      continue;
    }

    SnpCall call;
    call.position = *position;
    call.reference = static_cast<char>(std::toupper(static_cast<unsigned char>(fields[1].front())));
    call.alternate = static_cast<char>(std::toupper(static_cast<unsigned char>(fields[2].front())));
    call.gene = fields[3];
    call.consequence = fields[4];
    call.protein = fields[5];
    call.residue = fields[6];
    call.phenotype = fields[7];
    call.pathogenicity = fields[8];
    call.references = split_semicolon(fields[9]);
    call.sources = split_semicolon(fields[10]);
    if (fields.size() > 11) {
      call.structure_id = fields[11];
    }
    if (fields.size() > 12) {
      call.structure_chain = fields[12];
    }
    if (fields.size() > 13) {
      call.structure_residue = parse_size(fields[13]).value_or(0);
    }
    if (fields.size() > 14) {
      call.structure_complex = fields[14];
    }
    if (fields.size() > 15) {
      call.clinvar_allele_id = fields[15];
    }
    if (fields.size() > 16) {
      call.mitomap_url = fields[16];
    }
    const auto key = snp_key(call.position, call.reference, call.alternate);
    auto [it, inserted] = annotations.try_emplace(key, call);
    if (!inserted) {
      auto& existing = it->second;
      const auto replace_if_present = [](std::string& target, const std::string& incoming) {
        if (!incoming.empty()) {
          target = incoming;
        }
      };
      replace_if_present(existing.gene, call.gene);
      replace_if_present(existing.consequence, call.consequence);
      replace_if_present(existing.protein, call.protein);
      replace_if_present(existing.residue, call.residue);
      replace_if_present(existing.phenotype, call.phenotype);
      replace_if_present(existing.pathogenicity, call.pathogenicity);
      replace_if_present(existing.structure_id, call.structure_id);
      replace_if_present(existing.structure_chain, call.structure_chain);
      replace_if_present(existing.structure_complex, call.structure_complex);
      replace_if_present(existing.clinvar_allele_id, call.clinvar_allele_id);
      replace_if_present(existing.mitomap_url, call.mitomap_url);
      if (call.structure_residue != 0) {
        existing.structure_residue = call.structure_residue;
      }
      existing.references.insert(existing.references.end(), call.references.begin(), call.references.end());
      existing.sources.insert(existing.sources.end(), call.sources.begin(), call.sources.end());
      for (auto* values : {&existing.references, &existing.sources}) {
        std::sort(values->begin(), values->end());
        values->erase(std::unique(values->begin(), values->end()), values->end());
      }
    }
  }

  return annotations;
}

[[nodiscard]] InputRecords parse_fastq(std::istream& input,
                                       const AnalysisConfig& config) {
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
      throw AnalysisError(AnalysisErrorCode::input_parse_failed, "FASTQ record is truncated");
    }
    strip_trailing_carriage_return(sequence);
    strip_trailing_carriage_return(plus);
    strip_trailing_carriage_return(quality);
    if (plus.empty() || plus.front() != '+') {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "FASTQ parser expected a separator starting with '+' for read '" +
                              header.substr(1) + "'");
    }
    if (sequence.empty()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "FASTQ sequence is empty for read '" + header.substr(1) + "'");
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

[[nodiscard]] InputRecords parse_sam(std::istream& input,
                                     const AnalysisConfig& config) {
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
      for (const auto& field : header_fields) {
        if (field.rfind("SN:", 0) == 0 && looks_like_nuclear_contig(field.substr(3))) {
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
                          "invalid SAM FLAG at line " + std::to_string(line_number));
    }
    if (!reference_start) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid SAM POS at line " + std::to_string(line_number));
    }
    if (!mapping_quality || *mapping_quality > std::numeric_limits<std::uint8_t>::max()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid SAM MAPQ at line " + std::to_string(line_number));
    }
    if (fields[9] != "*" && fields[10] != "*" && fields[9].size() != fields[10].size()) {
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
      for (const auto& operation : record.cigar_operations) {
        if (operation.code == 'M' || operation.code == 'I' || operation.code == 'S' ||
            operation.code == '=' || operation.code == 'X') {
          if (operation.length > std::numeric_limits<std::size_t>::max() - query_length) {
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
    record.aux_tags = parse_sam_aux_tags(fields);
    records.reads.push_back(std::move(record));
  }
  return records;
}

#ifdef MITO_HAS_HTSLIB
class MitoFileReader {
public:
  explicit MitoFileReader(std::string path) : path_(std::move(path)) {}

  [[nodiscard]] InputRecords read_all(const AnalysisConfig& config) const {
    htsFile* raw_file = sam_open(path_.c_str(), "r");
    if (raw_file == nullptr) {
      throw AnalysisError(AnalysisErrorCode::input_open_failed,
                          "htslib could not open alignment file: " + path_);
    }
    HtsFileHandle file(raw_file);

    bam_hdr_t* raw_header = sam_hdr_read(file.get());
    if (raw_header == nullptr) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "htslib could not read alignment header: " + path_);
    }
    BamHeaderHandle header(raw_header);

    bam1_t* raw_record = bam_init1();
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
    while ((read_status = sam_read1(file.get(), header.get(), record.get())) >= 0) {
      throw_if_cancelled(config);
      records.reads.push_back(convert_record(*header.get(), *record.get()));
    }
    if (read_status < -1) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "htslib failed while reading alignment records: " + path_);
    }
    return records;
  }

private:
  struct HtsFileHandle {
    explicit HtsFileHandle(htsFile* value) : value_(value) {}
    ~HtsFileHandle() { hts_close(value_); }
    HtsFileHandle(const HtsFileHandle&) = delete;
    HtsFileHandle& operator=(const HtsFileHandle&) = delete;
    [[nodiscard]] htsFile* get() const { return value_; }
    htsFile* value_;
  };

  struct BamHeaderHandle {
    explicit BamHeaderHandle(bam_hdr_t* value) : value_(value) {}
    ~BamHeaderHandle() { bam_hdr_destroy(value_); }
    BamHeaderHandle(const BamHeaderHandle&) = delete;
    BamHeaderHandle& operator=(const BamHeaderHandle&) = delete;
    [[nodiscard]] bam_hdr_t* get() const { return value_; }
    bam_hdr_t* value_;
  };

  struct BamRecordHandle {
    explicit BamRecordHandle(bam1_t* value) : value_(value) {}
    ~BamRecordHandle() { bam_destroy1(value_); }
    BamRecordHandle(const BamRecordHandle&) = delete;
    BamRecordHandle& operator=(const BamRecordHandle&) = delete;
    [[nodiscard]] bam1_t* get() const { return value_; }
    bam1_t* value_;
  };

  [[nodiscard]] static ReadRecord convert_record(const bam_hdr_t& header, const bam1_t& record) {
    ReadRecord out;
    out.id = bam_get_qname(&record);
    out.flags = record.core.flag;
    out.mapping_quality = record.core.qual;
    out.reference_start = record.core.pos < 0 ? 1 : static_cast<std::size_t>(record.core.pos) + 1;
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

  [[nodiscard]] static std::string cigar_string(const bam1_t& record) {
    std::string cigar;
    const auto* raw_cigar = bam_get_cigar(&record);
    for (std::uint32_t i = 0; i < record.core.n_cigar; ++i) {
      cigar += std::to_string(bam_cigar_oplen(raw_cigar[i]));
      cigar.push_back(bam_cigar_opchr(raw_cigar[i]));
    }
    return cigar;
  }

  [[nodiscard]] static std::string sequence_string(const bam1_t& record) {
    std::string sequence;
    sequence.reserve(static_cast<std::size_t>(record.core.l_qseq));
    const auto* seq = bam_get_seq(&record);
    for (int i = 0; i < record.core.l_qseq; ++i) {
      sequence.push_back(seq_nt16_str[bam_seqi(seq, i)]);
    }
    return sequence;
  }

  [[nodiscard]] static std::string quality_string(const bam1_t& record) {
    std::string qualities;
    qualities.reserve(static_cast<std::size_t>(record.core.l_qseq));
    const auto* qual = bam_get_qual(&record);
    for (int i = 0; i < record.core.l_qseq; ++i) {
      qualities.push_back(qual[i] == 0xFFU ? '!' : static_cast<char>(qual[i] + 33U));
    }
    return qualities;
  }

  [[nodiscard]] static std::map<std::string, std::string> aux_tags(const bam1_t& record) {
    std::map<std::string, std::string> tags;
    for (const char* tag : {"NM", "MD", "AS", "SA"}) {
      const auto* value = bam_aux_get(&record, tag);
      if (value == nullptr) {
        continue;
      }
      const char type = static_cast<char>(*value);
      if (type == 'Z' || type == 'H') {
        if (const char* text = bam_aux2Z(value)) {
          tags.emplace(tag, text);
        }
      } else if (type == 'A') {
        tags.emplace(tag, std::string(1, bam_aux2A(value)));
      } else if (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' ||
                 type == 'I') {
        tags.emplace(tag, std::to_string(bam_aux2i(value)));
      } else if (type == 'f' || type == 'd') {
        tags.emplace(tag, std::to_string(bam_aux2f(value)));
      }
    }
    return tags;
  }

  std::string path_;
};
#endif

[[nodiscard]] bool has_extension(const std::filesystem::path& path,
                                 std::initializer_list<std::string_view> extensions) {
  auto ext = path.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char value) {
    return static_cast<char>(std::tolower(value));
  });
  return std::any_of(extensions.begin(), extensions.end(), [&](std::string_view value) {
    return ext == value;
  });
}

[[nodiscard]] InputRecords read_input(const std::string& input_path,
                                      const AnalysisConfig& config) {
  const std::filesystem::path path(input_path);
  if (has_extension(path, {".fa", ".fas", ".fasta", ".fna"})) {
    throw AnalysisError(AnalysisErrorCode::input_format_unsupported,
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
    throw AnalysisError(
        AnalysisErrorCode::dependency_unavailable,
        "BAM/CRAM support requires htslib; install htslib development headers and rebuild");
  }
  return parse_fastq(input, config);
}

[[nodiscard]] std::string alignment_as_sa(const ReadRecord& read) {
  const char strand = (read.flags & 0x10U) != 0U ? '-' : '+';
  const auto nm = read.aux_tags.find("NM");
  return read.reference_name + "," + std::to_string(read.reference_start) + "," + strand + "," +
         (read.cigar.empty() ? "*" : read.cigar) + "," +
         std::to_string(static_cast<unsigned int>(read.mapping_quality)) + "," +
         (nm == read.aux_tags.end() ? "0" : nm->second) + ";";
}

[[nodiscard]] std::vector<ReadRecord> collapse_molecule_alignments(InputRecords& input) {
  if (!input.alignment_input) {
    return std::move(input.reads);
  }

  std::map<std::string, std::vector<std::size_t>> groups;
  for (std::size_t i = 0; i < input.reads.size(); ++i) {
    groups[input.reads[i].id].push_back(i);
  }

  std::vector<ReadRecord> molecules;
  molecules.reserve(groups.size());
  for (const auto& [_, indices] : groups) {
    std::size_t primary_index = indices.front();
    for (const auto index : indices) {
      const auto flags = input.reads[index].flags;
      if ((flags & (0x100U | 0x800U)) == 0U) {
        primary_index = index;
        break;
      }
    }

    ReadRecord primary = input.reads[primary_index];
    std::string sa = primary.aux_tags.contains("SA") ? primary.aux_tags.at("SA") : std::string{};
    for (const auto index : indices) {
      if (index == primary_index || (input.reads[index].flags & 0x800U) == 0U ||
          input.reads[index].reference_name.empty() || input.reads[index].reference_name == "*") {
        continue;
      }
      const auto encoded = alignment_as_sa(input.reads[index]);
      if (sa.find(encoded) == std::string::npos) {
        sa += encoded;
      }
    }
    if (!sa.empty()) {
      primary.aux_tags["SA"] = std::move(sa);
    }
    primary.aux_tags["MA"] = std::to_string(indices.size());
    molecules.push_back(std::move(primary));
  }
  return molecules;
}

[[nodiscard]] double mean_quality(const std::string& qualities) {
  if (qualities.empty()) {
    return 0.0;
  }
  const auto total = std::accumulate(
      qualities.begin(), qualities.end(), 0.0,
      [](double sum, char c) { return sum + static_cast<double>(std::max(0, c - 33)); });
  return total / static_cast<double>(qualities.size());
}

[[nodiscard]] double gc_fraction(std::string_view sequence) {
  if (sequence.empty()) {
    return 0.0;
  }
  std::size_t gc = 0;
  for (const char c : sequence) {
    const char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
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

[[nodiscard]] bool has_nuclear_supplementary_alignment(const ReadRecord& read) {
  const auto sa = read.aux_tags.find("SA");
  if (sa == read.aux_tags.end()) {
    return false;
  }
  std::size_t start = 0;
  while (start < sa->second.size()) {
    const auto comma = sa->second.find(',', start);
    if (comma == std::string::npos) {
      break;
    }
    if (looks_like_nuclear_contig(std::string_view(sa->second).substr(start, comma - start))) {
      return true;
    }
    const auto semicolon = sa->second.find(';', comma);
    if (semicolon == std::string::npos) {
      break;
    }
    start = semicolon + 1;
  }
  return false;
}

[[nodiscard]] std::vector<std::string> numt_evidence(const ReadRecord& read,
                                                     std::size_t reference_length) {
  std::vector<std::string> evidence;
  const std::string id_lower = [&] {
    std::string out = read.id;
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
      return static_cast<char>(std::tolower(c));
    });
    return out;
  }();

  if (id_lower.find("numt") != std::string::npos || id_lower.find("nuclear") != std::string::npos) {
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

[[nodiscard]] double numt_score(const std::vector<std::string>& evidence) {
  double score = 0.05;
  for (const auto& item : evidence) {
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

[[nodiscard]] bool passes_snp_alignment_filters(const ReadRecord& read,
                                                const AnalysisConfig& config) {
  return (read.flags & 0x4U) == 0U &&
         (read.flags & config.excluded_snp_flags) == 0U &&
         read.mapping_quality >= config.min_mapping_quality && read.reference_start != 0 &&
         !read.reference_name.empty() && read.reference_name != "*" &&
         !looks_like_nuclear_contig(read.reference_name);
}

[[nodiscard]] std::optional<std::uint8_t> phred_quality_at(const ReadRecord& read,
                                                           std::size_t query_index) {
  if (query_index >= read.sequence.size() || query_index >= read.qualities.size()) {
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

[[nodiscard]] std::vector<SnpCall> call_snps_from_alignment(const ReadRecord& read,
                                                            const std::string& reference,
                                                            const AnalysisConfig& config) {
  std::vector<SnpCall> snps;
  if (reference.empty() || read.sequence.empty() || read.cigar_operations.empty() ||
      !passes_snp_alignment_filters(read, config)) {
    return snps;
  }

  std::size_t read_cursor = 0;
  std::size_t reference_cursor = read.reference_start == 0 ? 1 : read.reference_start;

  for (const auto& operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;

    if (op == 'M' || op == '=' || op == 'X') {
      for (std::size_t offset = 0; offset < len && read_cursor + offset < read.sequence.size(); ++offset) {
        const std::size_t position = ((reference_cursor - 1 + offset) % reference.size()) + 1;
        const char reference_base = reference[position - 1];
        const char alternate_base = static_cast<char>(
            std::toupper(static_cast<unsigned char>(read.sequence[read_cursor + offset])));
        const auto quality = phred_quality_at(read, read_cursor + offset);
        if (quality && *quality >= config.min_base_quality && is_base(alternate_base) &&
            is_base(reference_base) && alternate_base != reference_base) {
          SnpCall call;
          call.position = position;
          call.reference = reference_base;
          call.alternate = alternate_base;
          snps.push_back(std::move(call));
        }
      }
      read_cursor += len;
      reference_cursor =
          ((reference_cursor - 1U + (len % reference.size())) % reference.size()) + 1U;
    } else if (op == 'I' || op == 'S') {
      read_cursor += len;
    } else if (op == 'D' || op == 'N') {
      reference_cursor =
          ((reference_cursor - 1U + (len % reference.size())) % reference.size()) + 1U;
    } else if (op == 'H' || op == 'P') {
      // Does not consume read or reference sequence.
    }
  }

  std::sort(snps.begin(), snps.end(), [](const auto& lhs, const auto& rhs) {
    return snp_key(lhs.position, lhs.reference, lhs.alternate) <
           snp_key(rhs.position, rhs.reference, rhs.alternate);
  });
  snps.erase(std::unique(snps.begin(), snps.end(),
                         [](const auto& lhs, const auto& rhs) {
                           return lhs.position == rhs.position && lhs.reference == rhs.reference &&
                                  lhs.alternate == rhs.alternate;
                         }),
             snps.end());

  return snps;
}

[[nodiscard]] std::vector<SnpCall> call_snps_from_reference_span(const ReadRecord& read,
                                                                 const std::string& reference,
                                                                 const AnalysisConfig& config) {
  std::vector<SnpCall> snps;
  if (reference.empty() || read.sequence.empty() || read.reference_start == 0 ||
      read.reference_name.empty() || !passes_snp_alignment_filters(read, config)) {
    return snps;
  }

  for (std::size_t i = 0; i < read.sequence.size(); ++i) {
    const std::size_t position = ((read.reference_start - 1 + i) % reference.size()) + 1;
    const char reference_base = reference[position - 1];
    const char alternate_base =
        static_cast<char>(std::toupper(static_cast<unsigned char>(read.sequence[i])));
    const auto quality = phred_quality_at(read, i);
    if (quality && *quality >= config.min_base_quality && is_base(alternate_base) &&
        is_base(reference_base) && alternate_base != reference_base) {
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
  const auto arrow = token.find('>', colon == std::string_view::npos ? 0 : colon + 1);
  if (colon == std::string_view::npos || arrow == std::string_view::npos ||
      arrow + 1 >= token.size()) {
    return std::nullopt;
  }

  const auto position = parse_size(token.substr(0, colon));
  if (!position || *position == 0) {
    return std::nullopt;
  }

  const char reference =
      static_cast<char>(std::toupper(static_cast<unsigned char>(token[colon + 1])));
  const char alternate =
      static_cast<char>(std::toupper(static_cast<unsigned char>(token[arrow + 1])));
  const auto valid_base = [](char base) {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
  };
  if (!valid_base(reference) || !valid_base(alternate) || reference == alternate) {
    return std::nullopt;
  }

  SnpCall call;
  call.position = *position;
  call.reference = reference;
  call.alternate = alternate;
  return call;
}

[[nodiscard]] std::vector<SnpCall> snp_tags_from_header(const ReadRecord& read) {
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

void merge_snps(std::vector<SnpCall>& snps, std::vector<SnpCall> tagged_snps) {
  snps.insert(snps.end(), std::make_move_iterator(tagged_snps.begin()),
              std::make_move_iterator(tagged_snps.end()));
  std::sort(snps.begin(), snps.end(), [](const auto& lhs, const auto& rhs) {
    return snp_key(lhs.position, lhs.reference, lhs.alternate) <
           snp_key(rhs.position, rhs.reference, rhs.alternate);
  });
  snps.erase(std::unique(snps.begin(), snps.end(),
                         [](const auto& lhs, const auto& rhs) {
                           return lhs.position == rhs.position && lhs.reference == rhs.reference &&
                                  lhs.alternate == rhs.alternate;
                         }),
             snps.end());
}

void apply_clinical_annotations(std::vector<SnpCall>& snps,
                                const std::map<std::string, SnpCall>& annotations) {
  for (auto& snp : snps) {
    const auto it = annotations.find(snp_key(snp.position, snp.reference, snp.alternate));
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

void write_string_array(std::ostringstream& out, const std::vector<std::string>& values) {
  out << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << quoted(values[i]);
  }
  out << "]";
}

void write_string_map(std::ostringstream& out, const std::map<std::string, std::string>& values) {
  out << "{";
  std::size_t index = 0;
  for (const auto& [key, value] : values) {
    if (index++ != 0) {
      out << ",";
    }
    out << quoted(key) << ":" << quoted(value);
  }
  out << "}";
}

[[nodiscard]] std::optional<SvCall> parse_sv_tag(std::string_view read_id, std::string_view token) {
  const auto colon = token.find(':');
  const auto dash = token.find('-', colon == std::string_view::npos ? 0 : colon + 1);
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
  std::transform(sv.type.begin(), sv.type.end(), sv.type.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
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
  sv.known_event = sv.type == "deletion" && sv.start >= 8400 && sv.start <= 8500 &&
                   sv.end >= 13400 && sv.end <= 13550;
  sv.id = sv.type + ":" + std::to_string(sv.start) + "-" + std::to_string(sv.end);
  return sv;
}

[[nodiscard]] std::vector<SvCall> sv_tags_from_header(const ReadRecord& read) {
  std::vector<SvCall> calls;
  std::istringstream iss(read.id);
  std::string token;
  while (iss >> token) {
    constexpr std::string_view prefix = "sv=";
    if (token.rfind(prefix, 0) == 0) {
      auto call = parse_sv_tag(read.id, std::string_view(token).substr(prefix.size()));
      if (call) {
        calls.push_back(std::move(*call));
      }
    }
  }
  return calls;
}

[[nodiscard]] std::vector<SvCall> parse_cigar_svs(const ReadRecord& read,
                                                  std::size_t sv_min_length,
                                                  std::size_t reference_length) {
  std::vector<SvCall> calls;
  if (read.cigar_operations.empty() || (read.flags & 0x4U) != 0U ||
      read.reference_start == 0 || read.reference_name.empty() || read.reference_name == "*" ||
      looks_like_nuclear_contig(read.reference_name)) {
    return calls;
  }

  std::size_t reference_cursor = read.reference_start == 0 ? 1 : read.reference_start;
  for (const auto& operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;
    if ((op == 'D' || op == 'N') && len >= sv_min_length) {
      if (len - 1U > std::numeric_limits<std::size_t>::max() - reference_cursor) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "CIGAR structural-variant coordinate overflow for read '" + read.id +
                                "'");
      }
      SvCall sv;
      sv.type = "deletion";
      sv.start = reference_cursor;
      sv.end = reference_cursor + len - 1;
      sv.length = len;
      sv.supporting_reads.emplace_back(read.id);
      sv.evidence_sources.emplace_back("cigar");
      sv.known_event = sv.start >= 8400 && sv.start <= 8500 && sv.end >= 13400 &&
                       sv.end <= 13550;
      sv.id = "deletion:" + std::to_string(sv.start) + "-" + std::to_string(sv.end);
      calls.push_back(std::move(sv));
    } else if (op == 'I' && len >= sv_min_length) {
      SvCall sv;
      sv.type = "insertion";
      sv.start = reference_cursor > 1U ? reference_cursor - 1U : reference_length;
      sv.end = sv.start;
      sv.length = len;
      sv.supporting_reads.emplace_back(read.id);
      sv.evidence_sources.emplace_back("cigar");
      sv.id = "insertion:" + std::to_string(sv.start) + "+" + std::to_string(len);
      calls.push_back(std::move(sv));
    } else if (op == 'S' && len >= sv_min_length) {
      SvCall sv;
      const bool leading =
          reference_cursor == (read.reference_start == 0 ? 1 : read.reference_start);
      sv.type = leading ? "soft_clip_left" : "soft_clip_right";
      sv.start = leading ? reference_cursor
                         : (reference_cursor > 1U ? reference_cursor - 1U : reference_length);
      sv.end = reference_cursor;
      sv.length = len;
      sv.supporting_reads.emplace_back(read.id);
      sv.evidence_sources.emplace_back("cigar");
      sv.id = sv.type + ":" + std::to_string(sv.start) + "+" + std::to_string(len);
      calls.push_back(std::move(sv));
    }

    if (op == 'M' || op == '=' || op == 'X' || op == 'D' || op == 'N') {
      if (len > std::numeric_limits<std::size_t>::max() - reference_cursor) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "CIGAR reference coordinate overflow for read '" + read.id + "'");
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
  std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return normalized == "mt" || normalized == "chrm" || normalized == "m" ||
         normalized == "nc_012920.1" || normalized == "nc_012920";
}

[[nodiscard]] std::optional<AlignmentSegment> make_alignment_segment(
    std::string reference_name,
    std::size_t reference_start,
    char strand,
    std::uint8_t mapping_quality,
    const std::vector<CigarOperation>& operations) {
  if (!is_mitochondrial_contig(reference_name) || reference_start == 0 || operations.empty()) {
    return std::nullopt;
  }
  std::size_t leading_clip = 0;
  std::size_t query_extent = 0;
  std::size_t query_aligned = 0;
  std::size_t reference_span = 0;
  bool seen_alignment = false;
  for (const auto& operation : operations) {
    if (operation.code == 'M' || operation.code == 'I' || operation.code == 'S' ||
        operation.code == 'H' || operation.code == '=' || operation.code == 'X') {
      if (operation.length > std::numeric_limits<std::size_t>::max() - query_extent) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment query extent overflow");
      }
      query_extent += operation.length;
    }
    if (!seen_alignment && (operation.code == 'S' || operation.code == 'H')) {
      if (operation.length > std::numeric_limits<std::size_t>::max() - leading_clip) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment query coordinate overflow");
      }
      leading_clip += operation.length;
      continue;
    }
    if (operation.code == 'M' || operation.code == '=' || operation.code == 'X') {
      seen_alignment = true;
      if (operation.length > std::numeric_limits<std::size_t>::max() - query_aligned ||
          operation.length > std::numeric_limits<std::size_t>::max() - reference_span) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment segment length overflow");
      }
      query_aligned += operation.length;
      reference_span += operation.length;
    } else if (operation.code == 'I') {
      seen_alignment = true;
      if (operation.length > std::numeric_limits<std::size_t>::max() - query_aligned) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment query length overflow");
      }
      query_aligned += operation.length;
    } else if (operation.code == 'D' || operation.code == 'N') {
      seen_alignment = true;
      if (operation.length > std::numeric_limits<std::size_t>::max() - reference_span) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment reference length overflow");
      }
      reference_span += operation.length;
    }
  }
  if (query_aligned == 0 || reference_span == 0) {
    return std::nullopt;
  }
  if (leading_clip > query_extent || query_aligned > query_extent - leading_clip) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "alignment clipping exceeds query extent");
  }
  const auto query_start =
      strand == '+' ? leading_clip : query_extent - leading_clip - query_aligned;
  if (reference_span - 1U > std::numeric_limits<std::size_t>::max() - reference_start ||
      query_aligned - 1U > std::numeric_limits<std::size_t>::max() - query_start) {
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

[[nodiscard]] std::vector<AlignmentSegment> alignment_segments(const ReadRecord& read) {
  std::vector<AlignmentSegment> segments;
  if (const auto primary = make_alignment_segment(
          read.reference_name, read.reference_start, (read.flags & 0x10U) != 0U ? '-' : '+',
          read.mapping_quality, read.cigar_operations)) {
    segments.push_back(*primary);
  }
  const auto sa = read.aux_tags.find("SA");
  if (sa == read.aux_tags.end()) {
    return segments;
  }
  for (const auto& encoded : split_semicolon(sa->second)) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= encoded.size()) {
      const auto comma = encoded.find(',', start);
      fields.emplace_back(encoded.substr(start, comma == std::string::npos ? std::string::npos
                                                                          : comma - start));
      if (comma == std::string::npos) {
        break;
      }
      start = comma + 1;
    }
    if (fields.size() < 5) {
      continue;
    }
    const auto position = parse_size(fields[1]);
    const auto mapq = parse_size(fields[4]);
    if (!position || !mapq || *mapq > std::numeric_limits<std::uint8_t>::max() ||
        (fields[2] != "+" && fields[2] != "-")) {
      continue;
    }
    try {
      const auto operations = parse_cigar(fields[3], read.id + " SA");
      if (const auto segment = make_alignment_segment(fields[0], *position, fields[2].front(),
                                                      static_cast<std::uint8_t>(*mapq), operations)) {
        segments.push_back(*segment);
      }
    } catch (const std::runtime_error&) {
      // Malformed optional SA evidence is ignored; the primary alignment remains usable.
    }
  }
  std::sort(segments.begin(), segments.end(), [](const auto& lhs, const auto& rhs) {
    return std::tie(lhs.query_start, lhs.query_end, lhs.reference_start, lhs.strand) <
           std::tie(rhs.query_start, rhs.query_end, rhs.reference_start, rhs.strand);
  });
  segments.erase(std::unique(segments.begin(), segments.end(), [](const auto& lhs, const auto& rhs) {
                   return lhs.query_start == rhs.query_start && lhs.query_end == rhs.query_end &&
                          lhs.reference_start == rhs.reference_start && lhs.reference_end == rhs.reference_end &&
                          lhs.strand == rhs.strand;
                 }),
                 segments.end());
  return segments;
}

[[nodiscard]] std::vector<SvCall> split_alignment_svs(const ReadRecord& read,
                                                      std::size_t reference_length,
                                                      std::size_t sv_min_length) {
  std::vector<SvCall> calls;
  const auto segments = alignment_segments(read);
  if (segments.size() < 2 || reference_length == 0) {
    return calls;
  }
  for (std::size_t i = 1; i < segments.size(); ++i) {
    const auto& previous = segments[i - 1];
    const auto& current = segments[i];
    const auto query_gap = static_cast<std::int64_t>(current.query_start) -
                           static_cast<std::int64_t>(previous.query_end) - 1;
    SvCall event;
    event.supporting_reads.push_back(read.id);
    event.evidence_sources.emplace_back("split_alignment");
    event.segment_count = segments.size();
    event.orientations.push_back(std::string(1, previous.strand) + "/" + current.strand);

    const auto leaving_position = previous.strand == '+' ? previous.reference_end
                                                         : previous.reference_start;
    const auto entering_position = current.strand == '+' ? current.reference_start
                                                         : current.reference_end;
    if (leaving_position > reference_length || entering_position > reference_length) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
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
          forward ? entering_position > leaving_position : entering_position < leaving_position;
      const auto linear_gap = follows_linear_reference
                                  ? (forward ? entering_position - leaving_position - 1U
                                             : leaving_position - entering_position - 1U)
                                  : 0U;
      const auto overlap = follows_linear_reference
                               ? 0U
                               : (forward ? leaving_position - entering_position + 1U
                                          : entering_position - leaving_position + 1U);
      const auto circular_gap = follows_linear_reference
                                    ? reference_length
                                    : (forward ? reference_length - leaving_position +
                                                     entering_position - 1U
                                               : leaving_position - 1U + reference_length -
                                                     entering_position);
      const bool crosses_origin = !follows_linear_reference && circular_gap < overlap;
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
                        event.start <= 8500 && event.end >= 13400 && event.end <= 13550;
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

void merge_sv(std::map<std::string, SvCall>& svs, SvCall sv) {
  const auto existing = svs.find(sv.id);
  if (existing == svs.end()) {
    for (auto* values : {&sv.supporting_reads, &sv.evidence_sources, &sv.orientations}) {
      std::sort(values->begin(), values->end());
      values->erase(std::unique(values->begin(), values->end()), values->end());
    }
    svs.emplace(sv.id, std::move(sv));
    return;
  }

  auto& merged = existing->second;
  if (merged.type != sv.type || merged.start != sv.start || merged.end != sv.end ||
      merged.length != sv.length) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "canonical SV ID collision for '" + sv.id + "'");
  }
  const auto merge_unique = [](std::vector<std::string>& target,
                               const std::vector<std::string>& incoming) {
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

[[nodiscard]] std::size_t effective_thread_count(std::size_t requested_threads,
                                                 std::size_t work_items) {
  if (work_items <= 1) {
    return 1;
  }
  const auto hardware_threads = std::max(1U, std::thread::hardware_concurrency());
  const std::size_t requested = requested_threads == 0 ? 1 : requested_threads;
  return std::max<std::size_t>(
      1, std::min<std::size_t>({requested, work_items, static_cast<std::size_t>(hardware_threads)}));
}

[[nodiscard]] FeatureExtractionResult extract_read_feature(
    const ReadRecord& read,
    const std::string& reference,
    const AnalysisConfig& config,
    const std::map<std::string, SnpCall>& clinical_annotations) {
  throw_if_cancelled(config);
  FeatureExtractionResult result;
  auto& feature = result.feature;
  feature.id = read.id;
  feature.length = read.sequence.size();
  feature.mean_quality = mean_quality(read.qualities);
  feature.numt_evidence = numt_evidence(read, reference.size());
  feature.numt_score = numt_score(feature.numt_evidence);
  feature.filtered_numt = config.filter_numt && feature.numt_score > config.numt_threshold;
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
    if (std::any_of(feature.snps.begin(), feature.snps.end(), [&](const auto& snp) {
          return snp.position == 0 || snp.position > reference.size();
        })) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "SNP coordinate is outside the reference for read '" + read.id + "'");
    }
    apply_clinical_annotations(feature.snps, clinical_annotations);
    if (config.allow_development_tags) {
      result.svs = sv_tags_from_header(read);
    }
    auto cigar_svs = parse_cigar_svs(read, config.sv_min_length, reference.size());
    auto split_svs = split_alignment_svs(read, reference.size(), config.sv_min_length);
    if (!split_svs.empty()) {
      cigar_svs.erase(std::remove_if(cigar_svs.begin(), cigar_svs.end(), [](const auto& sv) {
                         return sv.type == "soft_clip_left" || sv.type == "soft_clip_right";
                       }),
                       cigar_svs.end());
    }
    result.svs.insert(result.svs.end(), std::make_move_iterator(cigar_svs.begin()),
                      std::make_move_iterator(cigar_svs.end()));
    result.svs.insert(result.svs.end(), std::make_move_iterator(split_svs.begin()),
                      std::make_move_iterator(split_svs.end()));
    if (std::any_of(result.svs.begin(), result.svs.end(), [&](const auto& sv) {
          return sv.start == 0 || sv.end > reference.size();
        })) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "SV coordinate is outside the reference for read '" + read.id + "'");
    }

    feature.sv_ids.reserve(result.svs.size());
    for (const auto& sv : result.svs) {
      feature.sv_ids.push_back(sv.id);
    }
    std::sort(feature.sv_ids.begin(), feature.sv_ids.end());
    feature.sv_ids.erase(std::unique(feature.sv_ids.begin(), feature.sv_ids.end()),
                         feature.sv_ids.end());
  }

  return result;
}

[[nodiscard]] std::vector<FeatureExtractionResult> extract_features_parallel(
    const std::vector<ReadRecord>& reads,
    const std::string& reference,
    const AnalysisConfig& config,
    const std::map<std::string, SnpCall>& clinical_annotations) {
  std::vector<FeatureExtractionResult> results(reads.size());
  if (reads.empty()) {
    return results;
  }

  const std::size_t worker_count = effective_thread_count(config.threads, reads.size());
  if (worker_count == 1) {
    for (std::size_t i = 0; i < reads.size(); ++i) {
      throw_if_cancelled(config);
      results[i] = extract_read_feature(reads[i], reference, config, clinical_annotations);
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
          const std::size_t begin = next_index.fetch_add(kChunkSize, std::memory_order_relaxed);
          if (begin >= reads.size()) {
            break;
          }
          const std::size_t end = std::min(reads.size(), begin + kChunkSize);
          for (std::size_t i = begin; i < end; ++i) {
            throw_if_cancelled(config);
            results[i] = extract_read_feature(reads[i], reference, config, clinical_annotations);
          }
        }
      } catch (...) {
        errors[worker_id] = std::current_exception();
        failed.store(true, std::memory_order_relaxed);
      }
    });
  }

  for (auto& worker : workers) {
    worker.join();
  }
  for (const auto& error : errors) {
    if (error) {
      std::rethrow_exception(error);
    }
  }

  return results;
}

using FeatureTokens = std::vector<std::string>;

[[nodiscard]] FeatureTokens feature_tokens(const ReadFeature& feature) {
  FeatureTokens tokens;
  tokens.reserve(feature.snps.size() + feature.sv_ids.size() + 1U);
  for (const auto& snp : feature.snps) {
    tokens.push_back("snp:" + snp_key(snp.position, snp.reference, snp.alternate));
  }
  for (const auto& sv_id : feature.sv_ids) {
    tokens.push_back("sv:" + sv_id);
  }
  if (tokens.empty()) {
    tokens.push_back("length:" + std::to_string(feature.length / 100U));
  }
  std::sort(tokens.begin(), tokens.end());
  tokens.erase(std::unique(tokens.begin(), tokens.end()), tokens.end());
  return tokens;
}

[[nodiscard]] double token_distance(const FeatureTokens& lhs, const FeatureTokens& rhs) {
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
  return 1.0 - (static_cast<double>(intersection) / static_cast<double>(union_size));
}

struct TokenProfile {
  FeatureTokens tokens;
  std::vector<std::size_t> feature_indices;
};

struct Neighborhood {
  std::vector<std::size_t> profiles;
  std::size_t read_count = 0;
};

using InvertedTokenIndex = std::unordered_map<std::string_view, std::vector<std::size_t>>;

[[nodiscard]] InvertedTokenIndex build_inverted_token_index(
    const std::vector<TokenProfile>& profiles) {
  InvertedTokenIndex index;
  for (std::size_t profile_index = 0; profile_index < profiles.size(); ++profile_index) {
    for (const auto& token : profiles[profile_index].tokens) {
      index[token].push_back(profile_index);
    }
  }
  return index;
}

[[nodiscard]] Neighborhood region_query(const std::vector<TokenProfile>& profiles,
                                        const InvertedTokenIndex& inverted_index,
                                        std::vector<std::size_t>& candidate_marks,
                                        std::size_t& mark_epoch,
                                        std::size_t index,
                                        double epsilon) {
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
    for (const auto& token : profiles[index].tokens) {
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
    if (token_distance(profiles[index].tokens, profiles[candidate].tokens) <= epsilon) {
      neighborhood.profiles.push_back(candidate);
      neighborhood.read_count += profiles[candidate].feature_indices.size();
    }
  }
  return neighborhood;
}

void append_new_neighbors(std::vector<std::size_t>& queue,
                          std::vector<unsigned char>& queued,
                          const std::vector<std::size_t>& candidates) {
  for (const std::size_t candidate : candidates) {
    if (candidate >= queued.size() || queued[candidate] != 0U) {
      continue;
    }
    queued[candidate] = 1U;
    queue.push_back(candidate);
  }
}

void assign_dbscan_clusters(std::vector<ReadFeature>& features,
                            double epsilon,
                            std::size_t min_cluster_size,
                            const AnalysisConfig& config) {
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

    auto neighborhood = region_query(profiles, inverted_index, candidate_marks, mark_epoch, point,
                                     epsilon);
    if (neighborhood.read_count < min_points) {
      labels[point] = -1;
      continue;
    }

    const int cluster_id = next_cluster_id++;
    labels[point] = cluster_id;
    auto& neighbors = neighborhood.profiles;
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
      auto expanded = region_query(profiles, inverted_index, candidate_marks, mark_epoch, neighbor,
                                   epsilon);
      if (expanded.read_count >= min_points) {
        append_new_neighbors(neighbors, queued, expanded.profiles);
      }
    }
  }

  for (std::size_t i = 0; i < profiles.size(); ++i) {
    for (const auto feature_index : profiles[i].feature_indices) {
      auto& feature = features[feature_index];
      feature.cluster_id = labels[i];
      feature.outlier = labels[i] < 0;
    }
  }
}

#ifdef MITO_HAS_HDBSCAN_CPP
[[nodiscard]] bool assign_hdbscan_clusters(std::vector<ReadFeature>& features,
                                           std::size_t min_cluster_size,
                                           const AnalysisConfig& config) {
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
    for (const auto& token : token_sets.back()) {
      if (!token_columns.contains(token)) {
        token_columns.emplace(token, token_columns.size());
      }
    }
  }

  std::vector<std::vector<double>> dataset(
      token_sets.size(), std::vector<double>(std::max<std::size_t>(1, token_columns.size()), 0.0));
  for (std::size_t row = 0; row < token_sets.size(); ++row) {
    for (const auto& token : token_sets[row]) {
      dataset[row][token_columns.at(token)] = 1.0;
    }
  }

  try {
    throw_if_cancelled(config);
    Hdbscan hdbscan("");
    hdbscan.dataset = std::move(dataset);
    const auto min_points = static_cast<int>(std::max<std::size_t>(1, min_cluster_size));
    hdbscan.execute(min_points, min_points, "Euclidean");
    if (hdbscan.normalizedLabels_.size() != active_indices.size()) {
      return false;
    }
    for (std::size_t i = 0; i < active_indices.size(); ++i) {
      const int label = hdbscan.normalizedLabels_[i];
      auto& feature = features[active_indices[i]];
      feature.cluster_id = label <= 0 ? -1 : label - 1;
      feature.outlier = label <= 0;
    }
    return true;
  } catch (const std::bad_alloc&) {
    throw;
  } catch (...) {
    throw_if_cancelled(config);
    return false;
  }
}
#endif

void add_circular_coverage_span(std::vector<std::int64_t>& difference,
                                std::uint64_t& uniform_depth,
                                std::size_t start,
                                std::size_t length) {
  const std::size_t reference_length = difference.size() - 1U;
  if (reference_length == 0 || length == 0) {
    return;
  }

  const auto complete_cycles = static_cast<std::uint64_t>(length / reference_length);
  if (complete_cycles > std::numeric_limits<std::uint64_t>::max() - uniform_depth) {
    throw AnalysisError(AnalysisErrorCode::resource_exhausted, "coverage depth overflow");
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

[[nodiscard]] CoverageResult compute_coverage(const std::vector<ReadRecord>& reads,
                                              const std::vector<ReadFeature>& features,
                                              const AnalysisConfig& config,
                                              std::size_t reference_length,
                                              std::size_t bin_count = 180) {
  CoverageResult result;
  if (reference_length == 0 || bin_count == 0) {
    return result;
  }

  std::vector<std::int64_t> difference(reference_length + 1U, 0);
  std::uint64_t uniform_depth = 0;
  for (std::size_t read_index = 0; read_index < reads.size(); ++read_index) {
    throw_if_cancelled(config);
    const auto& read = reads[read_index];
    if (read_index >= features.size() || features[read_index].filtered_numt ||
        (read.flags & 0x4U) != 0U || read.reference_start == 0 ||
        read.reference_name.empty() || read.reference_name == "*" ||
        read.cigar_operations.empty() || looks_like_nuclear_contig(read.reference_name)) {
      continue;
    }

    std::size_t reference_cursor = ((read.reference_start - 1U) % reference_length) + 1U;
    for (const auto& operation : read.cigar_operations) {
      if (operation.code == 'M' || operation.code == '=' || operation.code == 'X') {
        add_circular_coverage_span(difference, uniform_depth, reference_cursor, operation.length);
      }
      if (operation.code == 'M' || operation.code == '=' || operation.code == 'X' ||
          operation.code == 'D' || operation.code == 'N') {
        reference_cursor =
            ((reference_cursor - 1U + (operation.length % reference_length)) % reference_length) +
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
    if (variable_depth > std::numeric_limits<std::uint64_t>::max() - uniform_depth) {
      throw AnalysisError(AnalysisErrorCode::resource_exhausted, "coverage depth overflow");
    }
    const auto depth = uniform_depth + variable_depth;
    if (depth > std::numeric_limits<std::size_t>::max()) {
      throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                          "coverage depth exceeds platform size limit");
    }
    site_depth[i] = depth;
    total_depth += static_cast<long double>(depth);
    result.max_depth = std::max(result.max_depth, static_cast<std::size_t>(depth));
    if (depth > 20U) {
      ++sites_gt20;
    }
  }

  result.mean_depth = static_cast<double>(total_depth / static_cast<long double>(reference_length));
  result.pct_sites_gt20x =
      (static_cast<double>(sites_gt20) * 100.0) / static_cast<double>(reference_length);

  auto& bins = result.bins;
  bins.reserve(bin_count);
  const std::size_t bin_width = std::max<std::size_t>(1, reference_length / bin_count);
  for (std::size_t i = 0; i < bin_count; ++i) {
    const std::size_t start = i * bin_width + 1;
    if (start > reference_length) {
      break;
    }
    const std::size_t end =
        i == bin_count - 1 ? reference_length
                           : std::min(reference_length, start + bin_width - 1U);
    const auto first = site_depth.begin() + static_cast<std::ptrdiff_t>(start - 1U);
    const auto last = site_depth.begin() + static_cast<std::ptrdiff_t>(end);
    const long double bin_total = std::accumulate(
        first, last, 0.0L,
        [](long double sum, std::uint64_t depth) { return sum + static_cast<long double>(depth); });
    const auto site_count = end - start + 1U;
    const auto mean_depth = static_cast<std::size_t>(
        (bin_total / static_cast<long double>(site_count)) + 0.5L);
    bins.push_back({start, end, mean_depth});
  }

  return result;
}

[[nodiscard]] std::map<std::string, SnpAggregate> aggregate_snps(
    const std::vector<ReadRecord>& reads,
    const std::vector<ReadFeature>& features,
    const AnalysisConfig& config,
    std::size_t reference_length) {
  std::map<std::string, SnpAggregate> aggregates;
  std::unordered_set<std::size_t> target_positions;
  for (const auto& feature : features) {
    if (feature.filtered_numt) {
      continue;
    }
    for (const auto& snp : feature.snps) {
      const auto key = snp_key(snp.position, snp.reference, snp.alternate);
      auto [it, inserted] = aggregates.try_emplace(key);
      if (inserted) {
        it->second.call = snp;
      }
      target_positions.insert(snp.position);
    }
  }

  std::unordered_map<std::size_t, std::vector<SnpAggregate*>> aggregates_by_position;
  aggregates_by_position.reserve(target_positions.size());
  for (auto& [_, aggregate] : aggregates) {
    aggregates_by_position[aggregate.call.position].push_back(&aggregate);
  }

  std::unordered_map<std::size_t, std::array<std::size_t, 4>> allele_depths;
  allele_depths.reserve(target_positions.size());
  for (std::size_t read_index = 0; read_index < reads.size(); ++read_index) {
    throw_if_cancelled(config);
    if (read_index >= features.size() || features[read_index].filtered_numt) {
      continue;
    }
    const auto& read = reads[read_index];
    if (!passes_snp_alignment_filters(read, config) || read.cigar_operations.empty() ||
        reference_length == 0) {
      continue;
    }

    std::unordered_map<std::size_t, char> molecule_alleles;
    molecule_alleles.reserve(target_positions.size());
    std::size_t query_cursor = 0;
    std::size_t reference_cursor = ((read.reference_start - 1U) % reference_length) + 1U;
    for (const auto& operation : read.cigar_operations) {
      const auto len = operation.length;
      const char op = operation.code;
      if (op == 'M' || op == '=' || op == 'X') {
        for (std::size_t offset = 0; offset < len && query_cursor + offset < read.sequence.size();
             ++offset) {
          const std::size_t position = ((reference_cursor - 1U + offset) % reference_length) + 1U;
          if (!target_positions.contains(position)) {
            continue;
          }
          const auto quality = phred_quality_at(read, query_cursor + offset);
          const char base = static_cast<char>(
              std::toupper(static_cast<unsigned char>(read.sequence[query_cursor + offset])));
          if (!quality || *quality < config.min_base_quality || !base_index(base)) {
            continue;
          }
          auto [it, inserted] = molecule_alleles.try_emplace(position, base);
          if (!inserted && it->second != base) {
            it->second = 'N';
          }
        }
        query_cursor += len;
        reference_cursor =
            ((reference_cursor - 1U + (len % reference_length)) % reference_length) + 1U;
      } else if (op == 'I' || op == 'S') {
        query_cursor += len;
      } else if (op == 'D' || op == 'N') {
        reference_cursor =
            ((reference_cursor - 1U + (len % reference_length)) % reference_length) + 1U;
      }
    }

    for (const auto& [position, base] : molecule_alleles) {
      if (const auto index = base_index(base)) {
        ++allele_depths[position][*index];
        const auto position_aggregates = aggregates_by_position.find(position);
        if (position_aggregates != aggregates_by_position.end()) {
          for (auto* aggregate : position_aggregates->second) {
            if (base == aggregate->call.alternate) {
              aggregate->supporting_reads.push_back(features[read_index].id);
            }
          }
        }
      }
    }
  }

  constexpr double z = 1.959963984540054;
  for (auto& [_, aggregate] : aggregates) {
    auto& support = aggregate.supporting_reads;
    std::sort(support.begin(), support.end());
    support.erase(std::unique(support.begin(), support.end()), support.end());
    const auto depths_it = allele_depths.find(aggregate.call.position);
    if (depths_it == allele_depths.end()) {
      continue;
    }
    const auto& depths = depths_it->second;
    aggregate.callable_depth = std::accumulate(depths.begin(), depths.end(), std::size_t{0});
    if (const auto alternate = base_index(aggregate.call.alternate)) {
      aggregate.alternate_depth = depths[*alternate];
    }
    if (const auto reference = base_index(aggregate.call.reference)) {
      aggregate.reference_depth = depths[*reference];
    }
    aggregate.other_depth = aggregate.callable_depth - aggregate.alternate_depth - aggregate.reference_depth;
    if (aggregate.callable_depth == 0) {
      continue;
    }
    const double n = static_cast<double>(aggregate.callable_depth);
    const double proportion = static_cast<double>(aggregate.alternate_depth) / n;
    aggregate.heteroplasmy = proportion;
    const double z2 = z * z;
    const double denominator = 1.0 + z2 / n;
    const double center = (proportion + z2 / (2.0 * n)) / denominator;
    const double margin =
        z * std::sqrt((proportion * (1.0 - proportion) / n) + (z2 / (4.0 * n * n))) /
        denominator;
    aggregate.ci95_low = std::max(0.0, center - margin);
    aggregate.ci95_high = std::min(1.0, center + margin);
  }
  return aggregates;
}

[[nodiscard]] std::optional<std::pair<std::size_t, char>> phylo_snv(std::string_view token) {
  std::size_t digit_count = 0;
  while (digit_count < token.size() &&
         std::isdigit(static_cast<unsigned char>(token[digit_count])) != 0) {
    ++digit_count;
  }
  if (digit_count == 0 || digit_count >= token.size()) {
    return std::nullopt;
  }
  const auto position = parse_size(token.substr(0, digit_count));
  const char base = static_cast<char>(std::toupper(static_cast<unsigned char>(token[digit_count])));
  if (!position || *position == 0 || !is_base(base)) {
    return std::nullopt;
  }
  const auto suffix = token.substr(digit_count + 1U);
  if (!suffix.empty() && suffix != "!") {
    return std::nullopt;
  }
  return std::pair{*position, base};
}

[[nodiscard]] std::vector<HaplogroupDefinition> load_haplogroups() {
  std::ifstream input(phylotree_path());
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open PhyloTree resource: " + phylotree_path());
  }
  std::vector<HaplogroupDefinition> definitions;
  std::vector<HaplogroupDefinition> stack;
  std::string line;
  while (std::getline(input, line)) {
    const auto opening = line.find("<haplogroup name=\"");
    if (opening != std::string::npos) {
      const auto name_start = opening + std::string_view("<haplogroup name=\"").size();
      const auto name_end = line.find('"', name_start);
      if (name_end == std::string::npos) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "malformed haplogroup name in PhyloTree resource");
      }
      HaplogroupDefinition node;
      node.name = line.substr(name_start, name_end - name_start);
      if (!stack.empty()) {
        node.mutations = stack.back().mutations;
      }
      stack.push_back(std::move(node));
    }

    const auto poly_start = line.find("<poly>");
    if (poly_start != std::string::npos && !stack.empty()) {
      const auto value_start = poly_start + std::string_view("<poly>").size();
      const auto value_end = line.find("</poly>", value_start);
      if (value_end != std::string::npos) {
        const auto mutation = trim_copy(std::string_view(line).substr(value_start, value_end - value_start));
        if (const auto snv = phylo_snv(mutation)) {
          auto mutation_it = std::lower_bound(
              stack.back().mutations.begin(), stack.back().mutations.end(), snv->first,
              [](const auto& item, std::size_t position) { return item.first < position; });
          if (mutation.ends_with('!')) {
            if (mutation_it != stack.back().mutations.end() && mutation_it->first == snv->first) {
              stack.back().mutations.erase(mutation_it);
            }
          } else {
            const auto encoded = std::to_string(snv->first) + std::string(1, snv->second);
            if (mutation_it != stack.back().mutations.end() && mutation_it->first == snv->first) {
              mutation_it->second = encoded;
            } else {
              stack.back().mutations.insert(mutation_it, {snv->first, encoded});
            }
          }
        }
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
  if (!stack.empty() || definitions.size() < 5000U) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "PhyloTree resource is incomplete or unbalanced");
  }
  return definitions;
}

[[nodiscard]] std::unordered_map<std::string, double> load_phylo_weights() {
  std::ifstream input(phylotree_weights_path());
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open PhyloTree weights: " + phylotree_weights_path());
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
    } catch (const std::exception&) {
      // Ignore a malformed optional weight and use the deterministic default.
    }
  }
  return weights;
}

[[nodiscard]] double mutation_weight(const std::unordered_map<std::string, double>& weights,
                                     const std::string& mutation) {
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

[[nodiscard]] std::map<int, ClusterHaplogroupAssignment> assign_haplogroups(
    const std::vector<ReadFeature>& features,
    const std::vector<HaplogroupDefinition>& definitions,
    const std::unordered_map<std::string, double>& weights,
    const AnalysisConfig& config) {
  std::map<int, std::vector<const ReadFeature*>> clusters;
  for (const auto& feature : features) {
    if (!feature.filtered_numt) {
      clusters[feature.cluster_id].push_back(&feature);
    }
  }
  std::map<int, ClusterHaplogroupAssignment> assignments;
  std::unordered_map<std::string, std::vector<std::size_t>> definitions_by_mutation;
  definitions_by_mutation.reserve(definitions.size());
  for (std::size_t definition_index = 0; definition_index < definitions.size();
       ++definition_index) {
    for (const auto& [_, mutation] : definitions[definition_index].mutations) {
      definitions_by_mutation[mutation].push_back(definition_index);
    }
  }
  std::vector<std::size_t> candidate_marks(definitions.size(), 0U);
  std::size_t candidate_epoch = 0;
  std::map<std::vector<std::pair<std::size_t, std::string>>, ClusterHaplogroupAssignment>
      assignment_cache;
  for (const auto& [cluster_id, members] : clusters) {
    throw_if_cancelled(config);
    std::map<std::string, std::size_t> observed_counts;
    for (const auto* member : members) {
      std::unordered_set<std::string> molecule_mutations;
      for (const auto& snp : member->snps) {
        molecule_mutations.insert(std::to_string(snp.position) +
                                  std::string(1, snp.alternate));
      }
      for (const auto& mutation : molecule_mutations) {
        ++observed_counts[mutation];
      }
    }
    std::map<std::size_t, std::string> observed;
    for (const auto& [mutation, count] : observed_counts) {
      if (count * 2U < members.size()) {
        continue;
      }
      if (const auto parsed = phylo_snv(mutation)) {
        observed[parsed->first] = mutation;
      }
    }
    const std::vector<std::pair<std::size_t, std::string>> observed_signature(observed.begin(),
                                                                              observed.end());
    if (const auto cached = assignment_cache.find(observed_signature);
        cached != assignment_cache.end()) {
      assignments.emplace(cluster_id, cached->second);
      continue;
    }

    struct RankedDefinition {
      const HaplogroupDefinition* definition = nullptr;
      double score = 0.0;
      std::size_t matched_count = 0;
    };
    double observed_weight = 0.0;
    for (const auto& [_, mutation] : observed) {
      observed_weight += mutation_weight(weights, mutation);
    }
    std::vector<RankedDefinition> ranked;
    std::vector<std::size_t> candidate_indices;
    if (candidate_epoch == std::numeric_limits<std::size_t>::max()) {
      std::fill(candidate_marks.begin(), candidate_marks.end(), 0U);
      candidate_epoch = 0;
    }
    ++candidate_epoch;
    for (const auto& [_, mutation] : observed) {
      const auto posting = definitions_by_mutation.find(mutation);
      if (posting == definitions_by_mutation.end()) {
        continue;
      }
      for (const auto definition_index : posting->second) {
        if (candidate_marks[definition_index] != candidate_epoch) {
          candidate_marks[definition_index] = candidate_epoch;
          candidate_indices.push_back(definition_index);
        }
      }
    }
    std::sort(candidate_indices.begin(), candidate_indices.end());
    ranked.reserve(candidate_indices.size());
    for (const auto definition_index : candidate_indices) {
      const auto& definition = definitions[definition_index];
      const double expected_weight = definition.expected_weight;
      double matched_weight = 0.0;
      std::size_t matched_count = 0;
      for (const auto& [position, mutation] : definition.mutations) {
        const auto observed_it = observed.find(position);
        if (observed_it != observed.end() && observed_it->second == mutation) {
          matched_weight += mutation_weight(weights, mutation);
          ++matched_count;
        }
      }
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
    std::sort(ranked.begin(), ranked.end(), [](const auto& lhs, const auto& rhs) {
      if (std::abs(lhs.score - rhs.score) > 1e-9) {
        return lhs.score > rhs.score;
      }
      if (lhs.matched_count != rhs.matched_count) {
        return lhs.matched_count > rhs.matched_count;
      }
      if (lhs.definition->mutations.size() != rhs.definition->mutations.size()) {
        return lhs.definition->mutations.size() < rhs.definition->mutations.size();
      }
      return lhs.definition->name < rhs.definition->name;
    });
    ClusterHaplogroupAssignment assignment;
    const auto keep = std::min<std::size_t>(3, ranked.size());
    for (std::size_t i = 0; i < keep; ++i) {
      HaplogroupCandidate candidate;
      candidate.name = ranked[i].definition->name;
      candidate.score = ranked[i].score;
      for (const auto& [position, mutation] : ranked[i].definition->mutations) {
        const auto observed_it = observed.find(position);
        if (observed_it != observed.end() && observed_it->second == mutation) {
          candidate.matched.push_back(mutation);
        } else {
          candidate.missing.push_back(mutation);
        }
      }
      for (const auto& [position, mutation] : observed) {
        const auto& expected_mutations = ranked[i].definition->mutations;
        const auto expected = std::lower_bound(
            expected_mutations.begin(), expected_mutations.end(), position,
            [](const auto& item, std::size_t value) { return item.first < value; });
        if (expected == expected_mutations.end() || expected->first != position ||
            expected->second != mutation) {
          candidate.extra.push_back(mutation);
        }
      }
      assignment.candidates.push_back(std::move(candidate));
    }
    if (!assignment.candidates.empty()) {
      assignment.best = assignment.candidates.front().name;
      assignment.quality = assignment.candidates.front().score;
    }
    if (assignment.candidates.size() > 1 && assignment.candidates[1].score > 0.0 &&
        assignment.candidates.front().score - assignment.candidates[1].score < 3.0 &&
        macro_haplogroup(assignment.candidates.front().name) !=
            macro_haplogroup(assignment.candidates[1].name)) {
      assignment.contamination_warning = true;
    }
    assignment_cache.emplace(observed_signature, assignment);
    assignments.emplace(cluster_id, std::move(assignment));
  }
  return assignments;
}

[[nodiscard]] std::string render_json(const std::string& input_path,
                                      const std::string& reference_path,
                                      const AnalysisConfig& config,
                                      std::string_view clustering_backend,
                                      bool alignment_input,
                                      bool has_nuclear_contigs,
                                      std::size_t input_record_count,
                                      const std::string& reference,
                                      const std::vector<ReadRecord>& reads,
                                      const std::vector<ReadFeature>& features,
                                      const std::map<std::string, SvCall>& svs,
                                      const CoverageResult& coverage_result,
                                      const std::map<std::string, SnpAggregate>& variants,
                                      const std::map<int, ClusterHaplogroupAssignment>& haplogroups,
                                      const std::vector<ResourceRecord>& resources) {
  const auto& coverage = coverage_result.bins;
  std::map<int, std::vector<const ReadFeature*>> clusters;
  for (const auto& feature : features) {
    if (!feature.filtered_numt) {
      clusters[feature.cluster_id].push_back(&feature);
    }
  }

  std::ostringstream out;
  out << std::fixed << std::setprecision(6);
  out << "{";
  out << "\"metadata\":{";
  out << "\"schema_version\":" << quoted(kResultSchemaVersion) << ",";
  out << "\"sv_event_schema_version\":" << quoted(kSvEventSchemaVersion) << ",";
  out << "\"engine_version\":" << quoted(kEngineVersion) << ",";
  out << "\"sample\":" << quoted(config.sample_name.empty()
                                      ? std::filesystem::path(input_path).stem().string()
                                      : config.sample_name)
      << ",";
  out << "\"input_path\":" << quoted(input_path) << ",";
  out << "\"reference_path\":" << quoted(reference_path.empty() ? bundled_rcrs_path() : reference_path)
      << ",";
  out << "\"reference_accession\":" << quoted(reference_path.empty() ? "NC_012920.1" : "custom")
      << ",";
  out << "\"reference_length\":" << reference.size() << ",";
  out << "\"threads\":" << effective_thread_count(config.threads, reads.size()) << ",";
  out << "\"requested_threads\":" << std::max<std::size_t>(1, config.threads) << ",";
  out << "\"calling_parameters\":{";
  out << "\"min_mapping_quality\":" << static_cast<unsigned int>(config.min_mapping_quality) << ",";
  out << "\"min_base_quality\":" << static_cast<unsigned int>(config.min_base_quality) << ",";
  out << "\"excluded_snp_flags\":" << config.excluded_snp_flags << ",";
  out << "\"numt_threshold\":" << config.numt_threshold << ",";
  out << "\"development_tags_enabled\":"
      << (config.allow_development_tags ? "true" : "false") << "},";
  out << "\"resources\":[";
  for (std::size_t i = 0; i < resources.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    const auto& resource = resources[i];
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
      << quoted("Aligned reads use CIGAR-aware reference comparison for SNP extraction") << ","
      << quoted("NUMT specificity requires competitive nuclear-plus-mitochondrial alignment") << ","
      << quoted("Haplogroups use weighted PhyloTree 17.3 lineage matching") << ","
      << quoted("SV event IDs use canonical mtDNA adjacency schema 1.0") << ","
      << quoted("Split and supplementary alignments are reconstructed as molecule event edges") << ","
      << quoted(clustering_backend)
      << "]";
  out << "},";

  out << "\"filter_stats\":{";
  const std::size_t filtered =
      static_cast<std::size_t>(std::count_if(features.begin(), features.end(), [](const auto& feature) {
        return feature.filtered_numt;
      }));
  out << "\"input_reads\":" << reads.size() << ",";
  out << "\"input_alignment_records\":" << input_record_count << ",";
  out << "\"input_molecules\":" << reads.size() << ",";
  out << "\"passed_reads\":" << (features.size() - filtered) << ",";
  out << "\"numt_filtered_reads\":" << filtered << ",";
  out << "\"numt_threshold\":" << config.numt_threshold << ",";
  out << "\"numt_assessment\":{"
      << "\"mode\":" << quoted(has_nuclear_contigs ? "competitive_alignment" :
                                                      (alignment_input ? "mt_only_or_unknown" : "unaligned_fastq"))
      << ",\"nuclear_contigs_present\":" << (has_nuclear_contigs ? "true" : "false")
      << ",\"specificity_assessable\":" << (has_nuclear_contigs ? "true" : "false") << "}";
  out << "},";

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

  const auto bins_gt20 =
      static_cast<std::size_t>(std::count_if(coverage.begin(), coverage.end(), [](const auto& bin) {
        return bin.depth > 20;
  }));
  const double pct_bins_gt20 =
      coverage.empty()
          ? 0.0
          : (static_cast<double>(bins_gt20) * 100.0) / static_cast<double>(coverage.size());
  std::map<std::uint8_t, std::size_t> mapq_histogram;
  for (std::size_t i = 0; i < reads.size(); ++i) {
    if (i < features.size() && !features[i].filtered_numt && (reads[i].flags & 0x4U) == 0U &&
        !reads[i].reference_name.empty() && reads[i].reference_name != "*") {
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
  for (const auto& [mapq, count] : mapq_histogram) {
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
  for (const auto& [_, sv] : svs) {
    if (sv_index++ != 0) {
      out << ",";
    }
    out << "{"
        << "\"id\":" << quoted(sv.id) << ","
        << "\"type\":" << quoted(sv.type) << ","
        << "\"start\":" << sv.start << ","
        << "\"end\":" << sv.end << ","
        << "\"length\":" << sv.length << ","
        << "\"known_event\":" << (sv.known_event ? "true" : "false") << ","
        << "\"evidence_source\":"
        << quoted(sv.evidence_sources.size() == 1U ? sv.evidence_sources.front() : "combined")
        << ",\"evidence_sources\":";
    write_string_array(out, sv.evidence_sources);
    out << ",\"segment_count\":" << sv.segment_count << ",";
    if (!sv.orientations.empty()) {
      out << "\"orientation\":"
          << quoted(sv.orientations.size() == 1U ? sv.orientations.front() : "mixed")
          << ",\"orientations\":";
      write_string_array(out, sv.orientations);
      out << ",";
    }
    out
        << "\"supporting_reads\":[";
    for (std::size_t i = 0; i < sv.supporting_reads.size(); ++i) {
      if (i != 0) {
        out << ",";
      }
      out << quoted(sv.supporting_reads[i]);
    }
    out << "]}";
  }
  out << "],";

  out << "\"variants\":[";
  std::size_t variant_index = 0;
  for (const auto& [_, aggregate] : variants) {
    if (variant_index++ != 0) {
      out << ",";
    }
    const auto& snp = aggregate.call;
    out << "{"
        << "\"position\":" << snp.position << ","
        << "\"ref\":" << quoted(std::string(1, snp.reference)) << ","
        << "\"alt\":" << quoted(std::string(1, snp.alternate)) << ","
        << "\"alt_depth\":" << aggregate.alternate_depth << ","
        << "\"ref_depth\":" << aggregate.reference_depth << ","
        << "\"other_depth\":" << aggregate.other_depth << ","
        << "\"callable_depth\":" << aggregate.callable_depth << ","
        << "\"heteroplasmy\":" << aggregate.heteroplasmy << ","
        << "\"ci95_low\":" << aggregate.ci95_low << ","
        << "\"ci95_high\":" << aggregate.ci95_high << ","
        << "\"supporting_reads\":";
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
    if (!snp.phenotype.empty() || !snp.pathogenicity.empty() || !snp.references.empty() ||
        !snp.sources.empty() || !snp.clinvar_allele_id.empty() || !snp.mitomap_url.empty()) {
      out << ",\"annotation\":{";
      bool field = false;
      const auto comma = [&]() {
        if (field) {
          out << ",";
        }
        field = true;
      };
      if (!snp.phenotype.empty()) {
        comma();
        out << "\"phenotype\":" << quoted(snp.phenotype);
      }
      if (!snp.pathogenicity.empty()) {
        comma();
        out << "\"pathogenicity\":" << quoted(snp.pathogenicity);
      }
      if (!snp.references.empty()) {
        comma();
        out << "\"references\":";
        write_string_array(out, snp.references);
      }
      if (!snp.sources.empty()) {
        comma();
        out << "\"sources\":";
        write_string_array(out, snp.sources);
        out << ",\"source\":" << quoted("local-cache");
      }
      if (!snp.clinvar_allele_id.empty()) {
        comma();
        out << "\"clinvar_allele_id\":" << quoted(snp.clinvar_allele_id);
      }
      if (!snp.mitomap_url.empty()) {
        comma();
        out << "\"mitomap_url\":" << quoted(snp.mitomap_url);
      }
      out << "}";
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
  out << "],";

  out << "\"clusters\":[";
  std::size_t cluster_index = 0;
  for (const auto& [cluster_id, members] : clusters) {
    if (cluster_index++ != 0) {
      out << ",";
    }
    std::map<std::string, std::size_t> sv_counts;
    std::vector<std::string> read_ids;
    for (const auto* member : members) {
      read_ids.push_back(member->id);
      for (const auto& sv_id : member->sv_ids) {
        ++sv_counts[sv_id];
      }
    }
    const auto assignment_it = haplogroups.find(cluster_id);
    const ClusterHaplogroupAssignment* assignment =
        assignment_it == haplogroups.end() ? nullptr : &assignment_it->second;
    out << "{"
        << "\"id\":" << cluster_id << ","
        << "\"label\":"
        << quoted(cluster_id < 0 ? "Outliers" : "Cluster " + std::to_string(cluster_id + 1)) << ","
        << "\"haplogroup\":" << quoted(assignment == nullptr ? "unassigned" : assignment->best) << ","
        << "\"haplogroup_assignment\":{";
    out << "\"resource\":" << quoted("phylotree-rcrs@17.3") << ","
        << "\"quality\":" << (assignment == nullptr ? 0.0 : assignment->quality) << ","
        << "\"contamination_warning\":"
        << (assignment != nullptr && assignment->contamination_warning ? "true" : "false") << ","
        << "\"candidates\":[";
    if (assignment != nullptr) {
      for (std::size_t candidate_index = 0; candidate_index < assignment->candidates.size();
           ++candidate_index) {
        if (candidate_index != 0) {
          out << ",";
        }
        const auto& candidate = assignment->candidates[candidate_index];
        out << "{\"name\":" << quoted(candidate.name) << ",\"score\":" << candidate.score
            << ",\"matched\":";
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
        << quoted(cluster_id < 0 ? "noise" : "feature-consensus-" + std::to_string(cluster_id))
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
    for (const auto& [sv_id, count] : sv_counts) {
      if (signature_index++ != 0) {
        out << ",";
      }
      out << "{"
          << "\"sv_id\":" << quoted(sv_id) << ","
          << "\"support\":" << count << "}";
    }
    out << "]}";
  }
  out << "],";

  out << "\"reads\":[";
  for (std::size_t i = 0; i < features.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    const auto& feature = features[i];
    out << "{"
        << "\"id\":" << quoted(feature.id) << ","
        << "\"length\":" << feature.length << ","
        << "\"mean_quality\":" << feature.mean_quality << ","
        << "\"numt_score\":" << feature.numt_score << ","
        << "\"filtered_numt\":" << (feature.filtered_numt ? "true" : "false") << ","
        << "\"numt_evidence\":";
    write_string_array(out, feature.numt_evidence);
    out << ","
        << "\"cluster_id\":" << feature.cluster_id << ","
        << "\"outlier\":" << (feature.outlier ? "true" : "false") << ","
        << "\"mapping_quality\":" << static_cast<int>(feature.mapping_quality) << ","
        << "\"flags\":" << feature.flags << ",";
    if (!feature.reference_name.empty()) {
      out << "\"reference_name\":" << quoted(feature.reference_name) << ",";
    }
    if (!feature.aux_tags.empty()) {
      out << "\"aux_tags\":";
      write_string_map(out, feature.aux_tags);
      out << ",";
    }
    out
        << "\"snps\":[";
    for (std::size_t snp_index = 0; snp_index < feature.snps.size(); ++snp_index) {
      if (snp_index != 0) {
        out << ",";
      }
      const auto& snp = feature.snps[snp_index];
      out << "{"
          << "\"position\":" << snp.position << ","
          << "\"ref\":" << quoted(std::string(1, snp.reference)) << ","
          << "\"alt\":" << quoted(std::string(1, snp.alternate));
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
      if (!snp.phenotype.empty() || !snp.pathogenicity.empty() || !snp.references.empty() ||
          !snp.sources.empty() || !snp.clinvar_allele_id.empty() || !snp.mitomap_url.empty()) {
        out << ",\"annotation\":{";
        bool annotation_field = false;
        const auto comma = [&]() {
          if (annotation_field) {
            out << ",";
          }
          annotation_field = true;
        };
        if (!snp.phenotype.empty()) {
          comma();
          out << "\"phenotype\":" << quoted(snp.phenotype);
        }
        if (!snp.pathogenicity.empty()) {
          comma();
          out << "\"pathogenicity\":" << quoted(snp.pathogenicity);
        }
        if (!snp.references.empty()) {
          comma();
          out << "\"references\":";
          write_string_array(out, snp.references);
        }
        if (!snp.sources.empty()) {
          comma();
          out << "\"sources\":";
          write_string_array(out, snp.sources);
          out << ",\"source\":" << quoted("local-cache");
        }
        if (!snp.clinvar_allele_id.empty()) {
          comma();
          out << "\"clinvar_allele_id\":" << quoted(snp.clinvar_allele_id);
        }
        if (!snp.mitomap_url.empty()) {
          comma();
          out << "\"mitomap_url\":" << quoted(snp.mitomap_url);
        }
        out << "}";
      }
      if (!snp.structure_id.empty()) {
        out << ",\"structure\":{";
        out << "\"structure_id\":" << quoted(snp.structure_id);
        if (!snp.structure_chain.empty()) {
          out << ",\"chain\":" << quoted(snp.structure_chain);
        }
        if (snp.structure_residue != 0) {
          out << ",\"residue_index\":" << snp.structure_residue;
        }
        if (!snp.residue.empty()) {
          out << ",\"residue_label\":" << quoted(snp.residue);
        }
        if (!snp.structure_complex.empty()) {
          out << ",\"complex\":" << quoted(snp.structure_complex);
        }
        out << "}";
      }
      out << "}";
    }
    out << "],\"sv_ids\":[";
    for (std::size_t sv_id_index = 0; sv_id_index < feature.sv_ids.size(); ++sv_id_index) {
      if (sv_id_index != 0) {
        out << ",";
      }
      out << quoted(feature.sv_ids[sv_id_index]);
    }
    out << "]}";
  }
  out << "]";
  out << "}";
  return out.str();
}

} // namespace

std::string AnalysisEngine::analyze(const std::string& input_path,
                                    const std::string& reference_path,
                                    const AnalysisConfig& config) const {
  validate_config(config);
  throw_if_cancelled(config);
  const auto reference = load_reference(reference_path);
  throw_if_cancelled(config);
  auto input = read_input(input_path, config);
  if (input.reads.empty()) {
    throw AnalysisError(AnalysisErrorCode::input_empty,
                        "input contains no read records: " + input_path);
  }
  const bool alignment_input = input.alignment_input;
  const bool has_nuclear_contigs = input.has_nuclear_contigs;
  const std::size_t input_record_count = input.reads.size();
  const auto reads = collapse_molecule_alignments(input);
  throw_if_cancelled(config);
  std::vector<ReadFeature> features;
  features.reserve(reads.size());
  std::map<std::string, SvCall> svs;
  static const auto clinical_annotations = load_clinical_annotations();
  static const auto resources = load_resource_manifest();
  static const auto haplogroup_weights = load_phylo_weights();
  static const auto haplogroup_definitions = [&] {
    auto definitions = load_haplogroups();
    for (auto& definition : definitions) {
      for (const auto& [_, mutation] : definition.mutations) {
        definition.expected_weight += mutation_weight(haplogroup_weights, mutation);
      }
    }
    return definitions;
  }();
  auto extraction_results = extract_features_parallel(reads, reference, config, clinical_annotations);

  for (auto& result : extraction_results) {
    throw_if_cancelled(config);
    for (auto& sv : result.svs) {
      merge_sv(svs, std::move(sv));
    }
    features.push_back(std::move(result.feature));
  }

  std::string_view clustering_backend = "Read clusters use DBSCAN over SNP/SV feature tokens";
#ifdef MITO_HAS_HDBSCAN_CPP
  if (assign_hdbscan_clusters(features, config.min_cluster_size, config)) {
    clustering_backend = "Read clusters use HDBSCAN-C++ over dense SNP/SV feature vectors";
  } else {
    assign_dbscan_clusters(features, config.cluster_epsilon, config.min_cluster_size, config);
    clustering_backend = "Read clusters use DBSCAN fallback after HDBSCAN-C++ adapter failure";
  }
#else
  assign_dbscan_clusters(features, config.cluster_epsilon, config.min_cluster_size, config);
#endif

  throw_if_cancelled(config);
  const auto coverage = compute_coverage(reads, features, config, reference.size());
  const auto variants = aggregate_snps(reads, features, config, reference.size());
  const auto haplogroups =
      assign_haplogroups(features, haplogroup_definitions, haplogroup_weights, config);
  return render_json(input_path, reference_path, config, clustering_backend, alignment_input,
                     has_nuclear_contigs, input_record_count, reference, reads, features, svs, coverage, variants,
                     haplogroups, resources);
}

} // namespace mito
