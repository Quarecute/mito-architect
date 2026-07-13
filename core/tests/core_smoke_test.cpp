#include "mito/analysis_engine.hpp"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#define REQUIRE(condition)                                                                         \
  do {                                                                                             \
    if (!(condition)) {                                                                            \
      throw std::runtime_error("test assertion failed: " #condition);                             \
    }                                                                                              \
  } while (false)

int main() {
  const auto path = std::filesystem::temp_directory_path() / "mito_core_smoke.fastq";
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
  config.sample_name = "smoke";
  config.allow_development_tags = true;
  const std::string json = engine.analyze(path.string(), "", config);

  REQUIRE(json.find("\"schema_version\":\"0.4\"") != std::string::npos);
  REQUIRE(json.find("\"clusters\"") != std::string::npos);
  REQUIRE(json.find("\"numt_filtered_reads\":1") != std::string::npos);
  REQUIRE(json.find("del:8483-13459") != std::string::npos);
  REQUIRE(json.find("\"known_event\":true") != std::string::npos);
  REQUIRE(json.find("\"pathogenicity\":\"pathogenic\"") != std::string::npos);
  REQUIRE(json.find("\"protein\":\"MT-ATP6\"") != std::string::npos);
  REQUIRE(json.find("\"structure_id\":\"8H9S\"") != std::string::npos);
  REQUIRE(json.find("\"chain\":\"N\"") != std::string::npos);
  REQUIRE(json.find("\"name\":\"phylotree-rcrs\"") != std::string::npos);
  REQUIRE(json.find("\"version\":\"17.3\"") != std::string::npos);

  std::filesystem::remove(path);

  const auto sam_path = std::filesystem::temp_directory_path() / "mito_core_smoke.sam";
  {
    std::ofstream out(sam_path);
    out << "@SQ\tSN:MT\tLN:16569\n";
    out << "@SQ\tSN:chr1\tLN:248956422\n";
    out << "aligned-snp\t0\tMT\t1\t60\t10M\t*\t0\t0\tGAGCACAGGT\tIIIIIIIIII\tNM:i:1\tMD:Z:2T7\n";
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
    out << "secondary-alt\t256\tMT\t1\t60\t10M\t*\t0\t0\tGAGCACAGGT\tIIIIIIIIII\n";
    out << "nuclear-primary\t0\tchr1\t1\t60\t10M\t*\t0\t0\tGAGCACAGGT\tIIIIIIIIII\n";
    out << "nuclear-SA\t0\tMT\t1\t60\t10M\t*\t0\t0\tGAGCACAGGT\tIIIIIIIIII\t"
        << "SA:Z:chr1,100,+,10M,60,0;\n";
  }
  const std::string quality_json = engine.analyze(quality_path.string(), "", config);
  REQUIRE(quality_json.find("\"alt_depth\":1") != std::string::npos);
  REQUIRE(quality_json.find("\"ref_depth\":1") != std::string::npos);
  REQUIRE(quality_json.find("\"callable_depth\":2") != std::string::npos);
  REQUIRE(quality_json.find("\"heteroplasmy\":0.500000") != std::string::npos);
  REQUIRE(quality_json.find("\"supporting_reads\":[\"high-alt\"]") != std::string::npos);
  REQUIRE(quality_json.find("\"numt_filtered_reads\":2") != std::string::npos);
  REQUIRE(quality_json.find("\"mode\":\"competitive_alignment\"") != std::string::npos);
  REQUIRE(quality_json.find("primary_nuclear_alignment") != std::string::npos);
  REQUIRE(quality_json.find("supplementary_nuclear_alignment") != std::string::npos);
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
        << std::string(100, 'A') << "\t" << std::string(100, 'I') << "\tNM:i:0\n";
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
  const std::string split_json = engine.analyze(split_path.string(), "", config);
  REQUIRE(split_json.find("\"input_alignment_records\":8") != std::string::npos);
  REQUIRE(split_json.find("\"input_molecules\":7") != std::string::npos);
  REQUIRE(split_json.find("deletion:150-299") != std::string::npos);
  REQUIRE(split_json.find("inversion:549-849") != std::string::npos);
  REQUIRE(split_json.find("duplication:1020-1049") != std::string::npos);
  REQUIRE(split_json.find("circular_origin:16549-20") != std::string::npos);
  REQUIRE(split_json.find("\"evidence_sources\":[\"cigar\",\"split_alignment\"]") !=
          std::string::npos);
  REQUIRE(split_json.find("\"orientations\":[\"+/+\",\"-/-\"]") != std::string::npos);
  REQUIRE(split_json.find("\"supporting_reads\":[\"molecule-cigar-equivalent\","
                          "\"molecule-reverse-deletion\",\"molecule-split\"]") !=
          std::string::npos);
  REQUIRE(split_json.find("\"segment_count\":2") != std::string::npos);
  std::filesystem::remove(split_path);

  const auto haplogroup_path =
      std::filesystem::temp_directory_path() / "mito_core_haplogroup.fastq";
  {
    std::ofstream out(haplogroup_path);
    out << "@haplogroup-test snp=15314:G>A\nACGT\n+\nIIII\n";
  }
  const std::string haplogroup_json = engine.analyze(haplogroup_path.string(), "", config);
  REQUIRE(haplogroup_json.find("\"haplogroup\":\"H2a2a1a\"") != std::string::npos);
  REQUIRE(haplogroup_json.find("\"resource\":\"phylotree-rcrs@17.3\"") != std::string::npos);
  REQUIRE(haplogroup_json.find("\"quality\":100.000000") != std::string::npos);
  REQUIRE(haplogroup_json.find("\"matched\":[\"15314A\"]") != std::string::npos);
  std::filesystem::remove(haplogroup_path);

  const auto duplicate_sv_path =
      std::filesystem::temp_directory_path() / "mito_core_duplicate_sv.fastq";
  {
    std::ofstream out(duplicate_sv_path);
    out << "@read-a sv=del:100-120\nACGT\n+\nIIII\n";
    out << "@read-b sv=del:100-120\nTGCA\n+\nIIII\n";
  }
  const std::string duplicate_sv_json = engine.analyze(duplicate_sv_path.string(), "", config);
  REQUIRE(duplicate_sv_json.find("\"supporting_reads\":[\"read-a sv=del:100-120\",\"read-b "
                                 "sv=del:100-120\"]") != std::string::npos);
  std::filesystem::remove(duplicate_sv_path);

  const auto development_tag_path =
      std::filesystem::temp_directory_path() / "mito_core_development_tags.fastq";
  {
    std::ofstream out(development_tag_path);
    out << "@tagged snp=8993:T>G sv=del:100-120\nACGT\n+\nIIII\n";
  }
  mito::AnalysisConfig production_config;
  const std::string production_json =
      engine.analyze(development_tag_path.string(), "", production_config);
  REQUIRE(production_json.find("\"development_tags_enabled\":false") != std::string::npos);
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
  REQUIRE(permissive_numt_json.find("\"numt_filtered_reads\":0") != std::string::npos);
  REQUIRE(permissive_numt_json.find("\"numt_threshold\":0.900000") != std::string::npos);
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
  const std::string edge_case_json = engine.analyze(edge_case_path.string(), "", production_config);
  REQUIRE(edge_case_json.find("\"position\":3,\"ref\":\"T\",\"alt\":\"C\","
                              "\"alt_depth\":1,\"ref_depth\":1,\"other_depth\":1,"
                              "\"callable_depth\":3") != std::string::npos);
  REQUIRE(edge_case_json.find("\"position\":3,\"ref\":\"T\",\"alt\":\"G\","
                              "\"alt_depth\":1,\"ref_depth\":1,\"other_depth\":1,"
                              "\"callable_depth\":3") != std::string::npos);
  REQUIRE(edge_case_json.find("\"position\":1,\"ref\":\"G\",\"alt\":\"A\"") !=
          std::string::npos);
  std::filesystem::remove(edge_case_path);

  const auto molecule_case_path =
      std::filesystem::temp_directory_path() / "mito_core_molecule_evidence.sam";
  {
    std::ofstream out(molecule_case_path);
    out << "@SQ\tSN:NC_012920.1\tLN:16569\n";
    out << "conflict\t0\tNC_012920.1\t1\t60\t3M16566D3M\t*\t0\t0\tGAAGAT\tIIIIII\n";
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
  REQUIRE(molecule_case_json.find("\"position\":3,\"ref\":\"T\",\"alt\":\"A\","
                                  "\"alt_depth\":2,\"ref_depth\":1,\"other_depth\":0,"
                                  "\"callable_depth\":3") != std::string::npos);
  REQUIRE(molecule_case_json.find("\"supporting_reads\":[\"alt\",\"reverse\"]") !=
          std::string::npos);
  std::filesystem::remove(molecule_case_path);

  const auto invalid_sam_path =
      std::filesystem::temp_directory_path() / "mito_core_invalid_cigar.sam";
  {
    std::ofstream out(invalid_sam_path);
    out << "@SQ\tSN:NC_012920.1\tLN:16569\n";
    out << "bad-cigar\t0\tNC_012920.1\t1\t60\t5M\t*\t0\t0\tACGT\tIIII\n";
  }
  bool rejected_invalid_sam = false;
  try {
    static_cast<void>(engine.analyze(invalid_sam_path.string(), "", production_config));
  } catch (const mito::AnalysisError& error) {
    rejected_invalid_sam = error.code() == mito::AnalysisErrorCode::input_parse_failed;
  }
  REQUIRE(rejected_invalid_sam);
  std::filesystem::remove(invalid_sam_path);

  const auto invalid_fastq_path =
      std::filesystem::temp_directory_path() / "mito_core_invalid.fastq";
  {
    std::ofstream out(invalid_fastq_path);
    out << "@invalid\nACGT\nnot-plus\nIII\n";
  }
  bool rejected_invalid_fastq = false;
  try {
    static_cast<void>(engine.analyze(invalid_fastq_path.string(), "", config));
  } catch (const mito::AnalysisError& error) {
    rejected_invalid_fastq = error.code() == mito::AnalysisErrorCode::input_parse_failed;
  }
  REQUIRE(rejected_invalid_fastq);
  std::filesystem::remove(invalid_fastq_path);

  mito::AnalysisConfig invalid_config;
  invalid_config.numt_threshold = 1.01;
  bool rejected_invalid_config = false;
  try {
    static_cast<void>(engine.analyze(sam_path.string(), "", invalid_config));
  } catch (const mito::AnalysisError& error) {
    rejected_invalid_config =
        error.code() == mito::AnalysisErrorCode::invalid_configuration &&
        mito::analysis_error_code_name(error.code()) == "MITO-E1001";
  }
  REQUIRE(rejected_invalid_config);
  return 0;
}
