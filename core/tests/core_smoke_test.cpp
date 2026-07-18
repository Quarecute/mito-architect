#include "mito/analysis_engine.hpp"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#define REQUIRE(condition)                                                     \
  do {                                                                         \
    if (!(condition)) {                                                        \
      throw std::runtime_error("test assertion failed: " #condition);          \
    }                                                                          \
  } while (false)

int main() {
  const auto path =
      std::filesystem::temp_directory_path() / "mito_core_smoke.fastq";
  {
    std::ofstream out(path);
    out << "@read-hapA snp=8993:T>G sv=del:8483-13459\n";
    out << "ACGTACGTACGTACGTACGTACGTACGTACGT\n";
    out << "+\n";
    out << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    out << "@read-hapB\n";
    out << "TGCATGCATGCATGCATGCATGCATGCATGCA\n";
    out << "+\n";
    out << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    out << "@read-hapC\n";
    out << "GATTACAGATTACAGATTACAGATTACAGATT\n";
    out << "+\n";
    out << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    out << "@read-numt nuclear\n";
    out << "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n";
    out << "+\n";
    out << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
  }

  mito::AnalysisEngine engine;
  mito::AnalysisConfig config;
  config.sample_name = "smoke\n\"\\\x01";
  config.allow_development_tags = true;
  const std::string json = engine.analyze(path.string(), "", config);

  REQUIRE(json.find("\"schema_version\":\"0.5\"") != std::string::npos);
  REQUIRE(json.find("\"sample\":\"smoke\\n\\\"\\\\\\u0001\"") !=
          std::string::npos);
  REQUIRE(json.find("\"clinical_annotation_schema_version\":\"1.0\"") !=
          std::string::npos);
  REQUIRE(json.find("\"clusters\"") != std::string::npos);
  REQUIRE(json.find("\"numt_filtered_reads\":1") != std::string::npos);
  REQUIRE(json.find("del:8483-13459") != std::string::npos);
  REQUIRE(json.find("\"known_event\":true") != std::string::npos);
  REQUIRE(json.find("\"pathogenicity\":\"pathogenic\"") != std::string::npos);
  REQUIRE(json.find("\"conflict_status\":\"consistent\"") != std::string::npos);
  REQUIRE(json.find("\"assertions\":[{\"source\":\"ClinVar\"") !=
          std::string::npos);
  REQUIRE(json.find("\"protein\":\"MT-ATP6\"") != std::string::npos);
  REQUIRE(json.find("\"structure_id\":\"8H9S\"") != std::string::npos);
  REQUIRE(json.find("\"chain\":\"N\"") != std::string::npos);
  REQUIRE(json.find("\"name\":\"phylotree-rcrs\"") != std::string::npos);
  REQUIRE(json.find("\"version\":\"17.3\"") != std::string::npos);

  std::filesystem::remove(path);

  const auto duplicate_fastq_path = std::filesystem::temp_directory_path() /
                                    "mito_core_duplicate_names.fastq";
  {
    std::ofstream out(duplicate_fastq_path);
    out << "@duplicate\nGAT\n+\nIII\n";
    out << "@duplicate\nGAT\n+\nIII\n";
  }
  mito::AnalysisConfig duplicate_fastq_config;
  duplicate_fastq_config.result_schema = mito::ResultSchema::v0_6;
  const std::string duplicate_fastq_json =
      engine.analyze(duplicate_fastq_path.string(), "", duplicate_fastq_config);
  REQUIRE(duplicate_fastq_json.find("\"id\":\"duplicate#record:1\"") !=
          std::string::npos);
  REQUIRE(duplicate_fastq_json.find(
              "\"warnings\":[\"DUPLICATE_SOURCE_NAME_DISAMBIGUATED\"]") !=
          std::string::npos);
  std::filesystem::remove(duplicate_fastq_path);

  const auto sam_path =
      std::filesystem::temp_directory_path() / "mito_core_smoke.sam";
  {
    std::ofstream out(sam_path);
    out << "@SQ\tSN:MT\tLN:16569\n";
    out << "@SQ\tSN:chr1\tLN:248956422\n";
    out << "aligned-snp\t0\tMT\t1\t60\t10M\t*"
           "\t0\t0\tGAGCACAGGT\tIIIIIIIIII\tNM:i:1\tMD:Z:2T7\n";
    out << "aligned-cigar-sv\t0\tMT\t100\t42\t10M12D10M15I10M13S\t*\t0\t0\t"
        << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t"
        << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    out << "unmapped\t4\t*\t0\t55\t*\t*\t0\t0\tACGT\tIIII\n";
  }

  const std::string sam_json = engine.analyze(sam_path.string(), "", config);
  REQUIRE(sam_json.find("\"position\":3") != std::string::npos);
  REQUIRE(sam_json.find("\"ref\":\"T\"") != std::string::npos);
  REQUIRE(sam_json.find("\"alt\":\"G\"") != std::string::npos);
  REQUIRE(sam_json.find("\"mapping_quality\":60") != std::string::npos);
  REQUIRE(sam_json.find("\"mapq\":42") != std::string::npos);
  REQUIRE(sam_json.find("\"mapq\":55") == std::string::npos);
  REQUIRE(sam_json.find("\"NM\":\"1\"") != std::string::npos);
  REQUIRE(sam_json.find("\"MD\":\"2T7\"") != std::string::npos);
  REQUIRE(sam_json.find("deletion:110-121") != std::string::npos);
  REQUIRE(sam_json.find("insertion:131+15") != std::string::npos);
  REQUIRE(sam_json.find("soft_clip_right:141+13") != std::string::npos);
  REQUIRE(sam_json.find("\"max_depth\":1") != std::string::npos);
  REQUIRE(sam_json.find("\"variants\":[") != std::string::npos);
  REQUIRE(sam_json.find("\"position\":3") != std::string::npos);
  REQUIRE(sam_json.find("\"callable_depth\":1") != std::string::npos);

  std::filesystem::remove(sam_path);

  const auto quality_path =
      std::filesystem::temp_directory_path() / "mito_core_quality_filters.sam";
  {
    std::ofstream out(quality_path);
    out << "@SQ\tSN:MT\tLN:16569\n";
    out << "@SQ\tSN:chr1\tLN:248956422\n";
    out << "high-alt\t0\tMT\t1\t60\t10M\t*\t0\t0\tGAGCACAGGT\tIIIIIIIIII\n";
    out << "high-ref\t0\tMT\t1\t60\t10M\t*\t0\t0\tGATCACAGGT\tIIIIIIIIII\n";
    out << "low-base-alt\t0\tMT\t1\t60\t10M\t*\t0\t0\tGAGCACAGGT\tII!IIIIIII\n";
    out << "low-mapq-alt\t0\tMT\t1\t5\t10M\t*\t0\t0\tGAGCACAGGT\tIIIIIIIIII\n";
    out << "secondary-alt\t256\tMT\t1\t60\t10M\t*"
           "\t0\t0\tGAGCACAGGT\tIIIIIIIIII\n";
    out << "nuclear-primary\t0\tchr1\t1\t60\t10M\t*"
           "\t0\t0\tGAGCACAGGT\tIIIIIIIIII\n";
    out << "nuclear-SA\t0\tMT\t1\t60\t10M\t*\t0\t0\tGAGCACAGGT\tIIIIIIIIII\t"
        << "SA:Z:chr1,100,+,10M,60,0;\n";
  }
  const std::string quality_json =
      engine.analyze(quality_path.string(), "", config);
  REQUIRE(quality_json.find("\"alt_depth\":1") != std::string::npos);
  REQUIRE(quality_json.find("\"ref_depth\":1") != std::string::npos);
  REQUIRE(quality_json.find("\"callable_depth\":2") != std::string::npos);
  REQUIRE(quality_json.find("\"heteroplasmy\":0.500000") != std::string::npos);
  REQUIRE(quality_json.find("\"supporting_reads\":[\"high-alt\"]") !=
          std::string::npos);
  REQUIRE(quality_json.find(
              "\"strand_support\":{\"alt_forward\":1,\"alt_reverse\":0,"
              "\"ref_forward\":1,\"ref_reverse\":0") != std::string::npos);
  REQUIRE(quality_json.find("\"allele_quality\":{\"alternate\":{\"count\":1,"
                            "\"mean_phred\":40.000000") != std::string::npos);
  REQUIRE(quality_json.find("\"numt_filtered_reads\":2") != std::string::npos);
  REQUIRE(quality_json.find("\"mode\":\"competitive_alignment\"") !=
          std::string::npos);
  REQUIRE(quality_json.find("primary_nuclear_alignment") != std::string::npos);
  REQUIRE(quality_json.find("supplementary_nuclear_alignment") !=
          std::string::npos);
  std::filesystem::remove(quality_path);

  const auto split_path =
      std::filesystem::temp_directory_path() / "mito_core_split_alignment.sam";
  {
    std::ofstream out(split_path);
    out << "@SQ\tSN:MT\tLN:16569\n";
    out << "molecule-split\t0\tMT\t100\t60\t50M50S\t*\t0\t0\t"
        << std::string(100, 'A') << "\t" << std::string(100, 'I')
        << "\tSA:Z:MT,300,+,50S50M,60,0;\n";
    out << "molecule-split\t2048\tMT\t300\t60\t50S50M\t*\t0\t0\t"
        << std::string(100, 'A') << "\t" << std::string(100, 'I')
        << "\tNM:i:0\n";
    out << "molecule-cigar-equivalent\t0\tMT\t140\t60\t10M150D10M\t*\t0\t0\t"
        << std::string(20, 'A') << "\t" << std::string(20, 'I') << "\n";
    out << "molecule-reverse-deletion\t16\tMT\t300\t60\t50S50M\t*\t0\t0\t"
        << std::string(100, 'A') << "\t" << std::string(100, 'I')
        << "\tSA:Z:MT,100,-,50M50S,60,0;\n";
    out << "molecule-inversion\t0\tMT\t500\t60\t50M50S\t*\t0\t0\t"
        << std::string(100, 'A') << "\t" << std::string(100, 'I')
        << "\tSA:Z:MT,800,-,50M50S,60,0;\n";
    out << "molecule-duplication\t0\tMT\t1000\t60\t50M50S\t*\t0\t0\t"
        << std::string(100, 'A') << "\t" << std::string(100, 'I')
        << "\tSA:Z:MT,1020,+,50S50M,60,0;\n";
    out << "molecule-origin\t0\tMT\t16500\t60\t50M50S\t*\t0\t0\t"
        << std::string(100, 'A') << "\t" << std::string(100, 'I')
        << "\tSA:Z:MT,20,+,50S50M,60,0;\n";
    out << "molecule-reverse-origin\t16\tMT\t20\t60\t50S50M\t*\t0\t0\t"
        << std::string(100, 'A') << "\t" << std::string(100, 'I')
        << "\tSA:Z:MT,16500,-,50M50S,60,0;\n";
  }
  const std::string split_json =
      engine.analyze(split_path.string(), "", config);
  REQUIRE(split_json.find("\"input_alignment_records\":8") !=
          std::string::npos);
  REQUIRE(split_json.find("\"input_molecules\":7") != std::string::npos);
  REQUIRE(split_json.find("deletion:150-299") != std::string::npos);
  REQUIRE(split_json.find("inversion:549-849") != std::string::npos);
  REQUIRE(split_json.find("duplication:1020-1049") != std::string::npos);
  REQUIRE(split_json.find("circular_origin:16549-20") != std::string::npos);
  REQUIRE(
      split_json.find("\"evidence_sources\":[\"cigar\",\"split_alignment\"]") !=
      std::string::npos);
  REQUIRE(split_json.find("\"orientations\":[\"+/+\",\"-/-\"]") !=
          std::string::npos);
  REQUIRE(
      split_json.find("\"supporting_reads\":[\"molecule-cigar-equivalent\","
                      "\"molecule-reverse-deletion\",\"molecule-split\"]") !=
      std::string::npos);
  REQUIRE(split_json.find("\"segment_count\":2") != std::string::npos);
  std::filesystem::remove(split_path);

  const auto complex_path = std::filesystem::temp_directory_path() /
                            "mito_core_complex_alignment.sam";
  {
    std::ofstream out(complex_path);
    out << "@SQ\tSN:MT\tLN:16569\n";
    out << "complex-forward\t0\tMT\t100\t60\t10M20S\t*\t0\t0\t"
        << std::string(30, 'A') << "\t" << std::string(30, 'I')
        << "\tSA:Z:MT,300,+,10S10M10S,60,0;MT,600,+,20S10M,60,0;\n";
    out << "complex-reverse\t16\tMT\t600\t60\t20S10M\t*\t0\t0\t"
        << std::string(30, 'A') << "\t" << std::string(30, 'I')
        << "\tSA:Z:MT,300,-,10S10M10S,60,0;MT,100,-,10M20S,60,0;\n";
  }
  const std::string complex_json =
      engine.analyze(complex_path.string(), "", config);
  REQUIRE(complex_json.find("\"complex_sv_event_schema_version\":\"1.0\"") !=
          std::string::npos);
  REQUIRE(
      complex_json.find("complex:deletion:110-299@+/+|deletion:310-599@+/+") !=
      std::string::npos);
  REQUIRE(complex_json.find("\"junction_ids\":[\"deletion:110-299\","
                            "\"deletion:310-599\"]") != std::string::npos);
  REQUIRE(complex_json.find("\"junction_orientations\":[\"+/+\",\"+/+\"]") !=
          std::string::npos);
  REQUIRE(complex_json.find("\"supporting_reads\":[\"complex-forward\","
                            "\"complex-reverse\"]") != std::string::npos);
  REQUIRE(complex_json.find("\"complex_event_signature\":[{") !=
          std::string::npos);
  REQUIRE(
      complex_json.find("\"complex_event_ids\":[\"complex:deletion:110-299@+/+|"
                        "deletion:310-599@+/+\"]") != std::string::npos);
  std::filesystem::remove(complex_path);

  const auto haplogroup_path =
      std::filesystem::temp_directory_path() / "mito_core_haplogroup.fastq";
  {
    std::ofstream out(haplogroup_path);
    out << "@haplogroup-test snp=15314:G>A\nACGT\n+\nIIII\n";
  }
  const std::string haplogroup_json =
      engine.analyze(haplogroup_path.string(), "", config);
  REQUIRE(haplogroup_json.find("\"haplogroup\":\"H2a2a1a\"") !=
          std::string::npos);
  REQUIRE(haplogroup_json.find("\"resource\":\"phylotree-rcrs@17.3\"") !=
          std::string::npos);
  REQUIRE(haplogroup_json.find("\"quality\":100.000000") != std::string::npos);
  REQUIRE(haplogroup_json.find("\"matched\":[\"15314A\"]") !=
          std::string::npos);
  const auto profiled_haplogroup =
      engine.analyze_profiled(haplogroup_path.string(), "", config);
  REQUIRE(profiled_haplogroup.json == haplogroup_json);
  REQUIRE(profiled_haplogroup.timings.total_us > 0U);
  std::filesystem::remove(haplogroup_path);

  const auto normalized_indel_path =
      std::filesystem::temp_directory_path() /
      "mito_core_haplogroup_indel_normalization.sam";
  {
    std::ofstream out(normalized_indel_path);
    out << "@SQ\tSN:NC_012920.1\tLN:16569\n";
    out << "right-normalized\t0\tNC_012920.1\t519\t60\t4M2D49M4I7M\t*\t0\t0\t"
           "ACACCGCTGCTAACCCCATACCCCGAACCAACCAAACCCCAAAGACACCCCCCCCCCACAGTTT\t"
           "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    out << "left-aligned\t0\tNC_012920.1\t512\t60\t4M2D50M4I13M\t*\t0\t0\t"
           "AGCACACACACCGCTGCTAACCCCATACCCCGAACCAACCAAACCCCAAAGACACCCCCCCCCCACA"
           "GTTT\t"
           "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
           "IIII\n";
  }
  const std::string normalized_indel_json =
      engine.analyze(normalized_indel_path.string(), "", config);
  REQUIRE(normalized_indel_json.find("\"input_molecules\":2") !=
          std::string::npos);
  REQUIRE(normalized_indel_json.find(
              "\"observed_markers\":[\"523d\",\"524d\",\"573.1CCCC\"]") !=
          std::string::npos);
  REQUIRE(normalized_indel_json.find(
              "\"callable_ranges\":[{\"start\":519,\"end\":580}]") !=
          std::string::npos);
  REQUIRE(normalized_indel_json.find("516d") == std::string::npos);
  REQUIRE(normalized_indel_json.find("567.1CCCC") == std::string::npos);
  std::filesystem::remove(normalized_indel_path);

  const auto rotated_indel_path = std::filesystem::temp_directory_path() /
                                  "mito_core_haplogroup_rotated_indel.sam";
  {
    std::ofstream out(rotated_indel_path);
    out << "@SQ\tSN:NC_012920.1\tLN:16569\n";
    out << "rotated-insertion\t0\tNC_012920.1\t8266\t60\t5M9I25M\t*\t0\t0\t"
           "ATAGCACCCCCTCTACCCCCTCTACCCCCTCTAGAGCCC\t"
           "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
  }
  const std::string rotated_indel_json =
      engine.analyze(rotated_indel_path.string(), "", config);
  REQUIRE(rotated_indel_json.find("8289.1CCCCCTCTA") != std::string::npos);
  REQUIRE(rotated_indel_json.find("8270.1ACCCCCTCT") == std::string::npos);
  std::filesystem::remove(rotated_indel_path);

  const auto callable_majority_path = std::filesystem::temp_directory_path() /
                                      "mito_core_callable_majority.fastq";
  {
    std::ofstream out(callable_majority_path);
    out << "@range-a phylo=15314A phylo_range=100-110\nACGT\n+\nIIII\n";
    out << "@range-b phylo=15314A phylo_range=105-115\nACGT\n+\nIIII\n";
    out << "@range-c phylo=15314A phylo_range=108-120\nACGT\n+\nIIII\n";
  }
  const std::string callable_majority_json =
      engine.analyze(callable_majority_path.string(), "", config);
  REQUIRE(callable_majority_json.find(
              "\"callable_ranges\":[{\"start\":105,\"end\":115}]") !=
          std::string::npos);
  std::filesystem::remove(callable_majority_path);

  const auto callable_origin_path =
      std::filesystem::temp_directory_path() / "mito_core_callable_origin.sam";
  {
    std::ofstream out(callable_origin_path);
    out << "@SQ\tSN:MT\tLN:16569\n";
    out << "origin-callable\t0\tMT\t16565\t60\t10M\t*"
           "\t0\t0\tCGATGGATCA\tIIIIIIIIII\n";
  }
  const std::string callable_origin_json =
      engine.analyze(callable_origin_path.string(), "", config);
  REQUIRE(callable_origin_json.find(
              "\"callable_ranges\":[{\"start\":1,\"end\":5},{\"start\":16565,"
              "\"end\":16569}]") != std::string::npos);
  std::filesystem::remove(callable_origin_path);

  const auto duplicate_sv_path =
      std::filesystem::temp_directory_path() / "mito_core_duplicate_sv.fastq";
  {
    std::ofstream out(duplicate_sv_path);
    out << "@read-a sv=del:100-120\nACGT\n+\nIIII\n";
    out << "@read-b sv=del:100-120\nTGCA\n+\nIIII\n";
  }
  const std::string duplicate_sv_json =
      engine.analyze(duplicate_sv_path.string(), "", config);
  REQUIRE(duplicate_sv_json.find(
              "\"supporting_reads\":[\"read-a sv=del:100-120\",\"read-b "
              "sv=del:100-120\"]") != std::string::npos);
  std::filesystem::remove(duplicate_sv_path);

  const auto development_tag_path = std::filesystem::temp_directory_path() /
                                    "mito_core_development_tags.fastq";
  {
    std::ofstream out(development_tag_path);
    out << "@tagged snp=8993:T>G sv=del:100-120\nACGT\n+\nIIII\n";
  }
  mito::AnalysisConfig production_config;
  const std::string production_json =
      engine.analyze(development_tag_path.string(), "", production_config);
  REQUIRE(production_json.find("\"development_tags_enabled\":false") !=
          std::string::npos);
  REQUIRE(production_json.find("\"svs\":[]") != std::string::npos);
  REQUIRE(production_json.find("\"evidence_source\":\"development_tag\"") ==
          std::string::npos);
  REQUIRE(production_json.find("\"position\":8993") == std::string::npos);
  std::filesystem::remove(development_tag_path);

  const auto threshold_path =
      std::filesystem::temp_directory_path() / "mito_core_numt_threshold.fastq";
  {
    std::ofstream out(threshold_path);
    out << "@possible-numt nuclear\nACGTACGT\n+\nIIIIIIII\n";
  }
  mito::AnalysisConfig permissive_numt_config;
  permissive_numt_config.numt_threshold = 0.90;
  const std::string permissive_numt_json =
      engine.analyze(threshold_path.string(), "", permissive_numt_config);
  REQUIRE(permissive_numt_json.find("\"numt_filtered_reads\":0") !=
          std::string::npos);
  REQUIRE(permissive_numt_json.find("\"numt_threshold\":0.900000") !=
          std::string::npos);
  std::filesystem::remove(threshold_path);

  const auto edge_case_path =
      std::filesystem::temp_directory_path() / "mito_core_edge_cases.sam";
  {
    std::ofstream out(edge_case_path);
    out << "@SQ\tSN:NC_012920.1\tLN:16569\n";
    out << "multi-g\t0\tNC_012920.1\t1\t60\t4M\t*\t0\t0\tGAGC\tIIII\n";
    out << "multi-c\t0\tNC_012920.1\t1\t60\t4M\t*\t0\t0\tGACC\tIIII\n";
    out << "multi-ref\t0\tNC_012920.1\t1\t60\t4M\t*\t0\t0\tGATC\tIIII\n";
    out << "origin-wrap\t0\tNC_012920.1\t16568\t60\t4M\t*\t0\t0\tTGAA\tIIII\n";
  }
  const std::string edge_case_json =
      engine.analyze(edge_case_path.string(), "", production_config);
  REQUIRE(
      edge_case_json.find("\"position\":3,\"ref\":\"T\",\"alt\":\"C\","
                          "\"alt_depth\":1,\"ref_depth\":1,\"other_depth\":1,"
                          "\"callable_depth\":3") != std::string::npos);
  REQUIRE(
      edge_case_json.find("\"position\":3,\"ref\":\"T\",\"alt\":\"G\","
                          "\"alt_depth\":1,\"ref_depth\":1,\"other_depth\":1,"
                          "\"callable_depth\":3") != std::string::npos);
  REQUIRE(edge_case_json.find("\"position\":1,\"ref\":\"G\",\"alt\":\"A\"") !=
          std::string::npos);
  std::filesystem::remove(edge_case_path);

  const auto molecule_case_path = std::filesystem::temp_directory_path() /
                                  "mito_core_molecule_evidence.sam";
  {
    std::ofstream out(molecule_case_path);
    out << "@SQ\tSN:NC_012920.1\tLN:16569\n";
    out << "conflict\t0\tNC_012920.1\t1\t60\t3M16566D3M\t*"
           "\t0\t0\tGAAGAT\tIIIIII\n";
    out << "alt\t0\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAA\tIII\n";
    out << "reverse\t16\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAA\tIII\n";
    out << "reference\t0\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAT\tIII\n";
    out << "missing-quality\t0\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAA\t*\n";
    out << "secondary\t256\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAA\tIII\n";
    out << "qc-fail\t512\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAA\tIII\n";
    out << "duplicate\t1024\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAA\tIII\n";
    out << "supplementary\t2048\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAA\tIII\n";
  }
  const std::string molecule_case_json =
      engine.analyze(molecule_case_path.string(), "", production_config);
  REQUIRE(molecule_case_json.find(
              "\"position\":3,\"ref\":\"T\",\"alt\":\"A\","
              "\"alt_depth\":2,\"ref_depth\":1,\"other_depth\":0,"
              "\"callable_depth\":3") != std::string::npos);
  REQUIRE(
      molecule_case_json.find("\"supporting_reads\":[\"alt\",\"reverse\"]") !=
      std::string::npos);
  REQUIRE(molecule_case_json.find("\"strand_bias_delta\":0.500000") !=
          std::string::npos);

  mito::AnalysisConfig evidence_config = production_config;
  evidence_config.result_schema = mito::ResultSchema::v0_6;
  evidence_config.threads = 1;
  evidence_config.evidence_page_size = 2;
  const std::string evidence_json =
      engine.analyze(molecule_case_path.string(), "", evidence_config);
  REQUIRE(evidence_json.find("\"schema_version\":\"0.6\"") !=
          std::string::npos);
  REQUIRE(evidence_json.find("\"layout\":\"paged_columnar_molecule_event\"") !=
          std::string::npos);
  REQUIRE(evidence_json.find("\"observation_storage\":"
                             "\"embedded_columnar_pages\"") !=
          std::string::npos);
  REQUIRE(evidence_json.find("\"observation_pages\":[") != std::string::npos);
  REQUIRE(evidence_json.find("\"observations\":[") == std::string::npos);
  REQUIRE(
      evidence_json.find("\"scope\":\"snv_indel_sv_complex_evidence_rc2\"") !=
      std::string::npos);
  REQUIRE(evidence_json.find("\"callability\":[") != std::string::npos);
  REQUIRE(evidence_json.find("\"missing_pair_state\":\"NOT_CALLABLE\"") !=
          std::string::npos);
  REQUIRE(evidence_json.find("\"id\":\"snv:3:T:A\"") != std::string::npos);
  REQUIRE(evidence_json.find("\"ALTERNATE\"") != std::string::npos);
  REQUIRE(evidence_json.find("\"REFERENCE\"") != std::string::npos);
  REQUIRE(evidence_json.find("\"LOW_QUALITY\"") != std::string::npos);
  REQUIRE(evidence_json.find("\"CONFLICT\"") != std::string::npos);
  REQUIRE(evidence_json.find("\"phase_links\":[") != std::string::npos);
  REQUIRE(evidence_json.find("\"architectures\":[]") != std::string::npos);
  REQUIRE(evidence_json.find("\"event_id\":\"snv:3:T:A\",\"type\":\"SNV\","
                             "\"position\":3") != std::string::npos);
  REQUIRE(evidence_json.find("\"alternate_event_ids\":[") !=
          std::string::npos);
  evidence_config.threads = 4;
  const std::string parallel_evidence_json =
      engine.analyze(molecule_case_path.string(), "", evidence_config);
  const auto evidence_slice = [](const std::string &value) {
    const auto start = value.find("\"evidence_encoding\":");
    const auto end = value.find("\"genes\":", start);
    REQUIRE(start != std::string::npos);
    REQUIRE(end != std::string::npos);
    return value.substr(start, end - start);
  };
  REQUIRE(evidence_slice(parallel_evidence_json) ==
          evidence_slice(evidence_json));

  evidence_config.max_evidence_observations = 1;
  bool rejected_evidence_overflow = false;
  try {
    static_cast<void>(
        engine.analyze(molecule_case_path.string(), "", evidence_config));
  } catch (const mito::AnalysisError &error) {
    rejected_evidence_overflow =
        error.code() == mito::AnalysisErrorCode::resource_exhausted;
  }
  REQUIRE(rejected_evidence_overflow);
  std::filesystem::remove(molecule_case_path);

  const auto phase_case_path =
      std::filesystem::temp_directory_path() / "mito_core_callable_phase.sam";
  {
    std::ofstream out(phase_case_path);
    out << "@SQ\tSN:NC_012920.1\tLN:16569\n";
    out << "phase-both\t0\tNC_012920.1\t1\t60\t12M\t*"
           "\t0\t0\tGAACGCAGGTCT\tIIIIIIIIIIII\n";
    out << "phase-a-only\t0\tNC_012920.1\t1\t60\t12M\t*"
           "\t0\t0\tGAACACAGGTCT\tIIIIIIIIIIII\n";
    out << "reference\t0\tNC_012920.1\t1\t60\t12M\t*"
           "\t0\t0\tGATCACAGGTCT\tIIIIIIIIIIII\n";
    out << "insertion-phase\t0\tNC_012920.1\t1\t60\t5M1I7M\t*"
           "\t0\t0\tGAACAGCAGGTCT\tIIIIIIIIIIIII\n";
    out << "insertion-other\t0\tNC_012920.1\t1\t60\t5M1I7M\t*"
           "\t0\t0\tGAACATCAGGTCT\tIIIIIIIIIIIII\n";
    out << "deletion-phase\t0\tNC_012920.1\t1\t60\t4M1D7M\t*"
           "\t0\t0\tGAACCAGGTCT\tIIIIIIIIIII\n";
    out << "split-gap\t0\tNC_012920.1\t1\t60\t5M\t*"
           "\t0\t0\tGAACA\tIIIII\n";
    out << "split-gap\t2048\tNC_012920.1\t6\t60\t7M\t*"
           "\t0\t0\tCAGGTCT\tIIIIIII\n";
  }
  mito::AnalysisConfig phase_config;
  phase_config.result_schema = mito::ResultSchema::v0_6;
  const std::string phase_json =
      engine.analyze(phase_case_path.string(), "", phase_config);
  REQUIRE(phase_json.find("\"id\":\"indel:insertion:5:G\"") !=
          std::string::npos);
  REQUIRE(phase_json.find("\"id\":\"indel:insertion:5:T\"") !=
          std::string::npos);
  REQUIRE(phase_json.find("\"id\":\"indel:deletion:5-5:A\"") !=
          std::string::npos);
  REQUIRE(phase_json.find("\"event_id\":\"indel:insertion:5:G\","
                          "\"type\":\"SMALL_INSERTION\"") !=
          std::string::npos);
  REQUIRE(phase_json.find("\"event_a_id\":\"snv:3:T:A\"") != std::string::npos);
  REQUIRE(phase_json.find("\"assessability\":\"COMPLETE_FOR_BOTH_EVENTS\"") !=
          std::string::npos);
  REQUIRE(phase_json.find("\"both_alternate\":1") != std::string::npos);
  REQUIRE(phase_json.find("\"supporting_molecule_indices\":[") !=
          std::string::npos);
  REQUIRE(phase_json.find("\"uncertain_molecule_indices\":[") !=
          std::string::npos);
  REQUIRE(phase_json.find("\"qc_flags\":[") != std::string::npos);
  REQUIRE(phase_json.find("\"phase_null_model\":"
                          "\"independent_marginals_within_jointly_callable\"") !=
          std::string::npos);
  REQUIRE(phase_json.find("\"cigar_alternative_small_indel\"") !=
          std::string::npos);
  phase_config.max_phase_links = 1U;
  bool rejected_phase_overflow = false;
  try {
    static_cast<void>(
        engine.analyze(phase_case_path.string(), "", phase_config));
  } catch (const mito::AnalysisError &error) {
    rejected_phase_overflow =
        error.code() == mito::AnalysisErrorCode::resource_exhausted;
  }
  REQUIRE(rejected_phase_overflow);
  std::filesystem::remove(phase_case_path);

  const auto assembly_case_path = std::filesystem::temp_directory_path() /
                                  "mito_core_molecule_assembly.sam";
  {
    std::ofstream out(assembly_case_path);
    out << "@SQ\tSN:NC_012920.1\tLN:16569\n";
    out << "multiple-primary\t0\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAT\tIII\n";
    out << "multiple-primary\t16\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAT\tIII\n";
    out << "secondary-only\t256\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAT\tIII\n";
  }
  mito::AnalysisConfig assembly_config;
  assembly_config.result_schema = mito::ResultSchema::v0_6;
  const std::string assembly_json =
      engine.analyze(assembly_case_path.string(), "", assembly_config);
  REQUIRE(assembly_json.find("\"input_alignment_records\":3") !=
          std::string::npos);
  REQUIRE(assembly_json.find("\"input_molecules\":2") != std::string::npos);
  REQUIRE(assembly_json.find(
              "\"assembly_status\":\"first_of_multiple_primaries\"") !=
          std::string::npos);
  REQUIRE(
      assembly_json.find("\"warnings\":[\"MULTIPLE_PRIMARY_ALIGNMENTS\"]") !=
      std::string::npos);
  REQUIRE(
      assembly_json.find("\"assembly_status\":\"fallback_without_primary\"") !=
      std::string::npos);
  REQUIRE(assembly_json.find("\"warnings\":[\"NO_PRIMARY_ALIGNMENT\"]") !=
          std::string::npos);
  REQUIRE(assembly_json.find("\"role\":\"secondary\"") != std::string::npos);
  std::filesystem::remove(assembly_case_path);

  const auto protocol_case_path = std::filesystem::temp_directory_path() /
                                  "mito_core_protocol_identity.sam";
  {
    std::ofstream out(protocol_case_path);
    out << "@SQ\tSN:NC_012920.1\tLN:16569\n";
    out << "tag-a\t0\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAA\tIII"
           "\tMI:Z:M1\tRX:Z:AAA\tDX:Z:duplex\n";
    out << "tag-b\t0\tNC_012920.1\t1\t55\t3M\t*\t0\t0\tGAA\tIII"
           "\tMI:Z:M1\tRX:Z:AAA\tDX:Z:duplex\n";
    out << "missing\t0\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAT\tIII"
           "\tRX:Z:CCC\tDX:Z:simplex\n";
  }
  mito::AnalysisConfig protocol_config;
  protocol_config.result_schema = mito::ResultSchema::v0_6;
  protocol_config.molecule_id_tag = "MI";
  protocol_config.umi_tag = "RX";
  protocol_config.duplex_tag = "DX";
  const auto protocol_json =
      engine.analyze(protocol_case_path.string(), "", protocol_config);
  REQUIRE(protocol_json.find("\"id\":\"MI:M1\"") != std::string::npos);
  REQUIRE(protocol_json.find("\"identity_policy\":\"sam_tag:MI\"") !=
          std::string::npos);
  REQUIRE(protocol_json.find("\"source_qnames\":[\"tag-a\",\"tag-b\"]") !=
          std::string::npos);
  REQUIRE(protocol_json.find("\"assembly_status\":"
                             "\"explicit_tag_multi_qname_group\"") !=
          std::string::npos);
  REQUIRE(protocol_json.find("\"MISSING_MOLECULE_ID_TAG\"") !=
          std::string::npos);
  REQUIRE(protocol_json.find("\"evidence_eligible_molecules\":1") !=
          std::string::npos);
  std::filesystem::remove(protocol_case_path);

  const auto invalid_sam_path =
      std::filesystem::temp_directory_path() / "mito_core_invalid_cigar.sam";
  {
    std::ofstream out(invalid_sam_path);
    out << "@SQ\tSN:NC_012920.1\tLN:16569\n";
    out << "bad-cigar\t0\tNC_012920.1\t1\t60\t5M\t*\t0\t0\tACGT\tIIII\n";
  }
  bool rejected_invalid_sam = false;
  try {
    static_cast<void>(
        engine.analyze(invalid_sam_path.string(), "", production_config));
  } catch (const mito::AnalysisError &error) {
    rejected_invalid_sam =
        error.code() == mito::AnalysisErrorCode::input_parse_failed;
  }
  REQUIRE(rejected_invalid_sam);
  std::filesystem::remove(invalid_sam_path);

  const auto invalid_sa_path =
      std::filesystem::temp_directory_path() / "mito_core_invalid_sa.sam";
  {
    std::ofstream out(invalid_sa_path);
    out << "@SQ\tSN:NC_012920.1\tLN:16569\n";
    out << "bad-sa\t0\tNC_012920.1\t100\t60\t10M10S\t*\t0\t0\t"
           "AAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\t"
           "SA:Z:NC_012920.1,300,+,10S10M,60;\n";
  }
  bool rejected_invalid_sa = false;
  try {
    static_cast<void>(
        engine.analyze(invalid_sa_path.string(), "", production_config));
  } catch (const mito::AnalysisError &error) {
    rejected_invalid_sa =
        error.code() == mito::AnalysisErrorCode::input_parse_failed &&
        mito::analysis_error_code_name(error.code()) == "MITO-E1103";
  }
  REQUIRE(rejected_invalid_sa);
  std::filesystem::remove(invalid_sa_path);

  const auto invalid_fastq_path =
      std::filesystem::temp_directory_path() / "mito_core_invalid.fastq";
  {
    std::ofstream out(invalid_fastq_path);
    out << "@invalid\nACGT\nnot-plus\nIII\n";
  }
  bool rejected_invalid_fastq = false;
  try {
    static_cast<void>(engine.analyze(invalid_fastq_path.string(), "", config));
  } catch (const mito::AnalysisError &error) {
    rejected_invalid_fastq =
        error.code() == mito::AnalysisErrorCode::input_parse_failed;
  }
  REQUIRE(rejected_invalid_fastq);
  std::filesystem::remove(invalid_fastq_path);

  mito::AnalysisConfig invalid_config;
  invalid_config.numt_threshold = 1.01;
  bool rejected_invalid_config = false;
  try {
    static_cast<void>(engine.analyze(sam_path.string(), "", invalid_config));
  } catch (const mito::AnalysisError &error) {
    rejected_invalid_config =
        error.code() == mito::AnalysisErrorCode::invalid_configuration &&
        mito::analysis_error_code_name(error.code()) == "MITO-E1001";
  }
  REQUIRE(rejected_invalid_config);
  invalid_config = {};
  invalid_config.result_schema = mito::ResultSchema::v0_6;
  invalid_config.evidence_page_size = 0U;
  rejected_invalid_config = false;
  try {
    static_cast<void>(engine.analyze(sam_path.string(), "", invalid_config));
  } catch (const mito::AnalysisError &error) {
    rejected_invalid_config =
        error.code() == mito::AnalysisErrorCode::invalid_configuration;
  }
  REQUIRE(rejected_invalid_config);
  invalid_config = {};
  invalid_config.result_schema = mito::ResultSchema::v0_6;
  invalid_config.molecule_id_tag = "M!";
  rejected_invalid_config = false;
  try {
    static_cast<void>(engine.analyze(sam_path.string(), "", invalid_config));
  } catch (const mito::AnalysisError &error) {
    rejected_invalid_config =
        error.code() == mito::AnalysisErrorCode::invalid_configuration;
  }
  REQUIRE(rejected_invalid_config);
  return 0;
}
