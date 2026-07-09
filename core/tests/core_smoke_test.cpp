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
  const std::string json = engine.analyze(path.string(), "", config);

  REQUIRE(json.find("\"schema_version\":\"0.2\"") != std::string::npos);
  REQUIRE(json.find("\"clusters\"") != std::string::npos);
  REQUIRE(json.find("\"numt_filtered_reads\":1") != std::string::npos);
  REQUIRE(json.find("del:8483-13459") != std::string::npos);
  REQUIRE(json.find("\"known_event\":true") != std::string::npos);
  REQUIRE(json.find("\"pathogenicity\":\"pathogenic\"") != std::string::npos);
  REQUIRE(json.find("\"protein\":\"MT-ATP6\"") != std::string::npos);

  std::filesystem::remove(path);

  const auto sam_path = std::filesystem::temp_directory_path() / "mito_core_smoke.sam";
  {
    std::ofstream out(sam_path);
    out << "@SQ\tSN:MT\tLN:16569\n";
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
  REQUIRE(sam_json.find("insertion:132+15") != std::string::npos);
  REQUIRE(sam_json.find("soft_clip_right:142+13") != std::string::npos);
  REQUIRE(sam_json.find("\"max_depth\":1") != std::string::npos);

  std::filesystem::remove(sam_path);

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

  const auto invalid_fastq_path =
      std::filesystem::temp_directory_path() / "mito_core_invalid.fastq";
  {
    std::ofstream out(invalid_fastq_path);
    out << "@invalid\nACGT\nnot-plus\nIII\n";
  }
  bool rejected_invalid_fastq = false;
  try {
    static_cast<void>(engine.analyze(invalid_fastq_path.string(), "", config));
  } catch (const std::runtime_error&) {
    rejected_invalid_fastq = true;
  }
  REQUIRE(rejected_invalid_fastq);
  std::filesystem::remove(invalid_fastq_path);
  return 0;
}
