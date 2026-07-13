mod report;

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};
use indicatif::{ProgressBar, ProgressStyle};
use mito_ffi::{AnalyzeOptions, MitoEngine};
use serde_json::Value;
use std::collections::BTreeMap;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command as ProcessCommand, Stdio};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

#[derive(Debug, Parser)]
#[command(name = "mito-cli")]
#[command(about = "Long-read mtDNA analysis and offline report generator")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    Analyze(AnalyzeArgs),
    Doctor,
    UpdateClinical(UpdateClinicalArgs),
    ValidateFixture(ValidateFixtureArgs),
    ValidateSvFixture(ValidateSvFixtureArgs),
}

#[derive(Debug, Args)]
struct AnalyzeArgs {
    #[arg(short = 'i', long = "input")]
    input: PathBuf,

    #[arg(short = 'o', long = "output", default_value = "output.html")]
    output: PathBuf,

    #[arg(short = 'r', long = "reference")]
    reference: Option<PathBuf>,

    #[arg(long = "json")]
    json: bool,

    #[arg(long = "vcf")]
    vcf: Option<PathBuf>,

    #[arg(long = "filter-numt", default_value_t = true)]
    filter_numt: bool,

    #[arg(long = "threads", default_value_t = 1)]
    threads: usize,

    #[arg(long = "min-mapq", default_value_t = 20)]
    min_mapping_quality: u8,

    #[arg(long = "min-base-quality", default_value_t = 10)]
    min_base_quality: u8,

    #[arg(long = "excluded-snp-flags", default_value_t = 3840)]
    excluded_snp_flags: u16,

    #[arg(long = "numt-threshold", default_value_t = 0.30)]
    numt_threshold: f64,

    #[arg(long = "allow-development-tags", default_value_t = false)]
    allow_development_tags: bool,

    #[arg(long = "update-clinical")]
    update_clinical: bool,
}

#[derive(Debug, Args)]
struct UpdateClinicalArgs {
    #[arg(long = "source")]
    source: Option<PathBuf>,

    #[arg(long = "output")]
    output: Option<PathBuf>,

    #[arg(long = "clinvar-live")]
    clinvar_live: bool,

    #[arg(long = "clinvar-gz")]
    clinvar_gz: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct ValidateFixtureArgs {
    #[arg(long = "input")]
    input: PathBuf,

    #[arg(long = "expected-vcf")]
    expected_vcf: PathBuf,

    #[arg(long = "expected-passed")]
    expected_passed: u64,

    #[arg(long = "expected-numt")]
    expected_numt: u64,

    #[arg(long = "expected-snp", required = true, action = clap::ArgAction::Append)]
    expected_snp: Vec<String>,

    #[arg(long = "expected-sv", action = clap::ArgAction::Append)]
    expected_sv: Vec<String>,

    #[arg(long = "expected-mapq", action = clap::ArgAction::Append)]
    expected_mapq: Vec<u8>,

    #[arg(long = "expected-aux", action = clap::ArgAction::Append)]
    expected_aux: Vec<String>,
}

#[derive(Debug, Args)]
struct ValidateSvFixtureArgs {
    #[arg(long = "input")]
    input: PathBuf,

    #[arg(long = "expected-json")]
    expected_json: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Analyze(args) => analyze(args),
        Command::Doctor => doctor(),
        Command::UpdateClinical(args) => update_clinical(args).map(|path| {
            println!("clinical annotation cache: {}", path.display());
            println!("set MITO_CLINICAL_ANNOTATIONS={} to use it", path.display());
        }),
        Command::ValidateFixture(args) => validate_fixture(args),
        Command::ValidateSvFixture(args) => validate_sv_fixture(args),
    }
}

fn doctor() -> Result<()> {
    let capabilities = MitoEngine::capabilities();
    println!(
        "mito-engine {} (result schema {}, error schema {})",
        capabilities.engine_version, capabilities.schema_version, capabilities.error_schema_version
    );
    println!(
        "BAM/CRAM reader: {}",
        if capabilities.htslib {
            "enabled"
        } else {
            "disabled"
        }
    );

    let mut missing = Vec::new();
    for (program, arguments, required) in [
        ("samtools", &["--version"][..], true),
        ("minimap2", &["--version"][..], true),
        ("bcftools", &["--version"][..], false),
        ("fasterq-dump", &["--version"][..], false),
    ] {
        match ProcessCommand::new(program).args(arguments).output() {
            Ok(output) if output.status.success() => {
                let version = String::from_utf8_lossy(&output.stdout)
                    .lines()
                    .find(|line| !line.trim().is_empty())
                    .unwrap_or("version unavailable")
                    .trim()
                    .to_owned();
                println!("{program}: {version}");
            }
            _ => {
                println!("{program}: missing");
                if required {
                    missing.push(program);
                }
            }
        }
    }

    if !capabilities.htslib {
        anyhow::bail!(
            "native core was built without htslib; clean and rebuild after installing htslib"
        );
    }
    if !missing.is_empty() {
        anyhow::bail!("missing required tools: {}", missing.join(", "));
    }
    println!("required native analysis capabilities are available");
    Ok(())
}

fn analyze(args: AnalyzeArgs) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("input file does not exist: {}", args.input.display());
    }
    if let Some(reference) = &args.reference {
        if !reference.exists() {
            anyhow::bail!("reference file does not exist: {}", reference.display());
        }
    }
    if args.update_clinical {
        let cache_path = update_clinical(UpdateClinicalArgs {
            source: None,
            output: None,
            clinvar_live: false,
            clinvar_gz: None,
        })?;
        // SAFETY: the CLI is single-threaded before the C++ engine is created.
        unsafe {
            std::env::set_var("MITO_CLINICAL_ANNOTATIONS", &cache_path);
        }
    }

    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::with_template("{spinner:.cyan} {msg}")
            .unwrap_or_else(|_| ProgressStyle::default_spinner()),
    );
    progress.enable_steady_tick(Duration::from_millis(80));
    progress.set_message("initializing C++ analysis engine");

    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    progress.set_message("analyzing reads");
    let result = engine
        .analyze_with_options(
            &args.input,
            args.reference.as_deref(),
            AnalyzeOptions {
                filter_numt: args.filter_numt,
                threads: args.threads,
                min_mapping_quality: args.min_mapping_quality,
                min_base_quality: args.min_base_quality,
                excluded_snp_flags: args.excluded_snp_flags,
                numt_threshold: args.numt_threshold,
                allow_development_tags: args.allow_development_tags,
            },
        )
        .context("analysis failed")?;

    if let Some(vcf_path) = &args.vcf {
        progress.set_message("writing VCF export");
        let vcf = render_vcf(&result)?;
        fs::write(vcf_path, vcf)
            .with_context(|| format!("failed to write {}", vcf_path.display()))?;
    }

    if args.json {
        progress.finish_and_clear();
        println!("{result}");
        return Ok(());
    }

    progress.set_message("rendering standalone report");
    let html = report::render_report(&result)?;
    fs::write(&args.output, html)
        .with_context(|| format!("failed to write {}", args.output.display()))?;
    progress.finish_with_message(format!("wrote {}", args.output.display()));

    Ok(())
}

fn update_clinical(args: UpdateClinicalArgs) -> Result<PathBuf> {
    if args.clinvar_live {
        return update_clinvar_cache(args.output, args.clinvar_gz);
    }

    let source = args.source.unwrap_or_else(default_bundled_clinical_tsv);
    if !source.exists() {
        anyhow::bail!(
            "clinical annotation source does not exist: {}",
            source.display()
        );
    }

    let explicit_output = args.output.is_some();
    let mut output = args.output.unwrap_or_else(|| {
        default_cache_dir()
            .join("mito-architect")
            .join("clinical_annotations.tsv")
    });
    if let Some(parent) = output.parent() {
        if let Err(error) = fs::create_dir_all(parent) {
            if explicit_output {
                return Err(error)
                    .with_context(|| format!("failed to create {}", parent.display()));
            }
            output = std::env::temp_dir()
                .join("mito-architect")
                .join("clinical_annotations.tsv");
            if let Some(fallback_parent) = output.parent() {
                fs::create_dir_all(fallback_parent)
                    .with_context(|| format!("failed to create {}", fallback_parent.display()))?;
            }
        }
    }
    fs::copy(&source, &output).with_context(|| {
        format!(
            "failed to copy clinical annotations from {} to {}",
            source.display(),
            output.display()
        )
    })?;
    Ok(output)
}

fn validate_fixture(args: ValidateFixtureArgs) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("fixture input does not exist: {}", args.input.display());
    }
    if !args.expected_vcf.exists() {
        anyhow::bail!(
            "expected VCF does not exist: {}",
            args.expected_vcf.display()
        );
    }

    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    let json = engine
        .analyze_with_options(
            &args.input,
            None,
            AnalyzeOptions {
                filter_numt: true,
                threads: 1,
                ..AnalyzeOptions::default()
            },
        )
        .context("analysis failed")?;
    let data: Value = serde_json::from_str(&json).context("analysis returned invalid JSON")?;
    assert_json_u64(
        &data,
        "/filter_stats/passed_reads",
        args.expected_passed,
        "passed read count",
    )?;
    assert_json_u64(
        &data,
        "/filter_stats/numt_filtered_reads",
        args.expected_numt,
        "NUMT filtered read count",
    )?;
    for expected_snp in &args.expected_snp {
        assert_snp_present(&data, expected_snp)?;
    }
    for expected_sv in &args.expected_sv {
        assert_sv_present(&data, expected_sv)?;
    }
    for expected_mapq in args.expected_mapq {
        assert_mapq_present(&data, expected_mapq)?;
    }
    for expected_aux in &args.expected_aux {
        assert_aux_present(&data, expected_aux)?;
    }

    let actual_vcf = normalize_text(&render_vcf(&json)?);
    let expected_vcf = normalize_text(
        &fs::read_to_string(&args.expected_vcf)
            .with_context(|| format!("failed to read {}", args.expected_vcf.display()))?,
    );
    if actual_vcf != expected_vcf {
        anyhow::bail!(
            "VCF mismatch for {}\nexpected:\n{}\nactual:\n{}",
            args.input.display(),
            expected_vcf,
            actual_vcf
        );
    }

    println!("fixture validation passed: {}", args.input.display());
    Ok(())
}

fn validate_sv_fixture(args: ValidateSvFixtureArgs) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("fixture input does not exist: {}", args.input.display());
    }
    if !args.expected_json.exists() {
        anyhow::bail!(
            "expected SV JSON does not exist: {}",
            args.expected_json.display()
        );
    }

    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    let json = engine
        .analyze_with_options(
            &args.input,
            None,
            AnalyzeOptions {
                filter_numt: true,
                threads: 1,
                ..AnalyzeOptions::default()
            },
        )
        .context("analysis failed")?;
    let data: Value = serde_json::from_str(&json).context("analysis returned invalid JSON")?;
    let actual = serde_json::json!({
        "schema_version": data.pointer("/metadata/sv_event_schema_version"),
        "input_alignment_records": data.pointer("/filter_stats/input_alignment_records"),
        "input_molecules": data.pointer("/filter_stats/input_molecules"),
        "passed_reads": data.pointer("/filter_stats/passed_reads"),
        "numt_filtered_reads": data.pointer("/filter_stats/numt_filtered_reads"),
        "svs": data.get("svs"),
    });
    let expected: Value = serde_json::from_str(
        &fs::read_to_string(&args.expected_json)
            .with_context(|| format!("failed to read {}", args.expected_json.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.expected_json.display()))?;
    if actual != expected {
        anyhow::bail!(
            "SV golden mismatch for {}\nexpected:\n{}\nactual:\n{}",
            args.input.display(),
            serde_json::to_string_pretty(&expected)?,
            serde_json::to_string_pretty(&actual)?
        );
    }

    println!("SV fixture validation passed: {}", args.input.display());
    Ok(())
}

fn assert_json_u64(data: &Value, pointer: &str, expected: u64, label: &str) -> Result<()> {
    let actual = data
        .pointer(pointer)
        .and_then(Value::as_u64)
        .with_context(|| format!("missing numeric field {pointer}"))?;
    if actual != expected {
        anyhow::bail!("{label} mismatch: expected {expected}, got {actual}");
    }
    Ok(())
}

fn assert_snp_present(data: &Value, expected: &str) -> Result<()> {
    let mut parts = expected.split(':');
    let position = parts
        .next()
        .context("expected SNP must be formatted as position:ref:alt")?
        .parse::<u64>()
        .context("expected SNP position is not numeric")?;
    let reference = parts
        .next()
        .context("expected SNP must include reference allele")?;
    let alternate = parts
        .next()
        .context("expected SNP must include alternate allele")?;
    if parts.next().is_some() {
        anyhow::bail!("expected SNP must be formatted as position:ref:alt");
    }

    for read in data
        .get("reads")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
    {
        for snp in read
            .get("snps")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
        {
            if snp.get("position").and_then(Value::as_u64) == Some(position)
                && snp.get("ref").and_then(Value::as_str) == Some(reference)
                && snp.get("alt").and_then(Value::as_str) == Some(alternate)
            {
                return Ok(());
            }
        }
    }
    anyhow::bail!("expected SNP not found: {expected}");
}

fn assert_sv_present(data: &Value, expected_id: &str) -> Result<()> {
    let found = data
        .get("svs")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .any(|sv| sv.get("id").and_then(Value::as_str) == Some(expected_id));
    if found {
        Ok(())
    } else {
        anyhow::bail!("expected SV not found: {expected_id}");
    }
}

fn assert_mapq_present(data: &Value, expected_mapq: u8) -> Result<()> {
    let found = data
        .pointer("/coverage_metrics/mapping_quality_histogram")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .any(|entry| {
            entry.get("mapq").and_then(Value::as_u64) == Some(u64::from(expected_mapq))
                && entry.get("count").and_then(Value::as_u64).unwrap_or(0) > 0
        });
    if found {
        Ok(())
    } else {
        anyhow::bail!("expected MAPQ not found in histogram: {expected_mapq}");
    }
}

fn assert_aux_present(data: &Value, expected: &str) -> Result<()> {
    let (key, value) = expected
        .split_once('=')
        .context("expected aux tag must be formatted as TAG=value")?;
    if key.len() != 2 || !key.chars().all(|c| c.is_ascii_alphanumeric()) {
        anyhow::bail!("expected aux tag key must be a two-character SAM tag");
    }

    for read in data
        .get("reads")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
    {
        if read
            .get("aux_tags")
            .and_then(Value::as_object)
            .and_then(|tags| tags.get(key))
            .and_then(Value::as_str)
            == Some(value)
        {
            return Ok(());
        }
    }
    anyhow::bail!("expected aux tag not found: {expected}");
}

fn normalize_text(value: &str) -> String {
    value.replace("\r\n", "\n").trim_end().to_string()
}

fn update_clinvar_cache(output: Option<PathBuf>, clinvar_gz: Option<PathBuf>) -> Result<PathBuf> {
    let explicit_output = output.is_some();
    let mut output = output.unwrap_or_else(|| {
        default_cache_dir()
            .join("mito-architect")
            .join("clinical_annotations.tsv")
    });

    if let Err(error) = try_prepare_parent(&output) {
        if explicit_output {
            return Err(error);
        }
        output = std::env::temp_dir()
            .join("mito-architect")
            .join("clinical_annotations.tsv");
        try_prepare_parent(&output)?;
    }

    let (source, remove_source) = match clinvar_gz {
        Some(path) => (path, false),
        None => (download_clinvar_summary()?, true),
    };
    let output_name = output
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("clinical_annotations.tsv");
    let temporary_output =
        output.with_file_name(format!(".{output_name}.tmp-{}", std::process::id()));
    let write_result = write_clinvar_mtdna_tsv(&source, &temporary_output).and_then(|()| {
        fs::rename(&temporary_output, &output).with_context(|| {
            format!(
                "failed to atomically replace {} with refreshed clinical cache",
                output.display()
            )
        })
    });
    if write_result.is_err() {
        let _ = fs::remove_file(&temporary_output);
    }
    let metadata_result = if write_result.is_ok() {
        write_clinical_cache_metadata(&output, &source, remove_source)
    } else {
        Ok(())
    };
    if remove_source {
        let _ = fs::remove_file(&source);
    }
    write_result?;
    metadata_result?;
    Ok(output)
}

fn write_clinical_cache_metadata(output: &Path, source: &Path, downloaded: bool) -> Result<()> {
    let retrieved_at_unix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system time is before the Unix epoch")?
        .as_secs();
    let metadata = serde_json::json!({
        "schema_version": 1,
        "resource": "ClinVar variant_summary mitochondrial SNVs",
        "source": if downloaded {
            "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz".to_string()
        } else {
            source.display().to_string()
        },
        "retrieved_at_unix": retrieved_at_unix,
        "source_sha256": sha256_file(source)?,
        "cache_sha256": sha256_file(output)?,
        "normalization": "GRCh38 MT/M single-nucleotide records; PositionVCF and VCF alleles; bundled curated rows appended and merged field-wise",
    });
    let file_name = output
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("clinical_annotations.tsv");
    let metadata_path = output.with_file_name(format!("{file_name}.metadata.json"));
    let temporary_path = metadata_path.with_extension(format!("json.tmp-{}", std::process::id()));
    fs::write(&temporary_path, serde_json::to_vec_pretty(&metadata)?)
        .with_context(|| format!("failed to write {}", temporary_path.display()))?;
    fs::rename(&temporary_path, &metadata_path).with_context(|| {
        format!(
            "failed to atomically replace clinical metadata {}",
            metadata_path.display()
        )
    })?;
    Ok(())
}

fn sha256_file(path: &Path) -> Result<String> {
    let output = ProcessCommand::new("sha256sum")
        .arg(path)
        .output()
        .with_context(|| format!("failed to start sha256sum for {}", path.display()))?;
    if !output.status.success() {
        anyhow::bail!("sha256sum failed for {}", path.display());
    }
    let text = String::from_utf8(output.stdout).context("sha256sum returned non-UTF-8 output")?;
    text.split_whitespace()
        .next()
        .filter(|value| value.len() == 64 && value.chars().all(|c| c.is_ascii_hexdigit()))
        .map(str::to_string)
        .context("sha256sum returned an invalid digest")
}

fn try_prepare_parent(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    Ok(())
}

fn download_clinvar_summary() -> Result<PathBuf> {
    let output = std::env::temp_dir().join(format!(
        "clinvar_variant_summary_{}.txt.gz",
        std::process::id()
    ));
    let status = ProcessCommand::new("curl")
        .args([
            "--fail",
            "--silent",
            "--show-error",
            "-L",
            "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
            "-o",
        ])
        .arg(&output)
        .status()
        .context("failed to start curl for ClinVar download")?;
    if !status.success() {
        anyhow::bail!("curl failed while downloading ClinVar variant_summary.txt.gz");
    }
    Ok(output)
}

fn write_clinvar_mtdna_tsv(source_gz: &Path, output: &Path) -> Result<()> {
    let mut gzip = ProcessCommand::new("gzip")
        .arg("-dc")
        .arg(source_gz)
        .stdout(Stdio::piped())
        .spawn()
        .with_context(|| format!("failed to decompress {}", source_gz.display()))?;
    let stdout = gzip.stdout.take().context("gzip did not provide stdout")?;
    let mut reader = BufReader::new(stdout);
    let mut header = String::new();
    reader
        .read_line(&mut header)
        .context("ClinVar summary is empty")?;
    let columns = header
        .trim_start_matches('#')
        .trim_end()
        .split('\t')
        .enumerate()
        .map(|(index, name)| (name.to_string(), index))
        .collect::<std::collections::HashMap<_, _>>();
    let mut out = fs::File::create(output)
        .with_context(|| format!("failed to create {}", output.display()))?;
    writeln!(
        out,
        "position\tref\talt\tgene\tconsequence\tprotein\tresidue\tphenotype\tpathogenicity\treferences\tsources\tstructure_id\tstructure_chain\tstructure_residue\tstructure_complex\tclinvar_allele_id\tmitomap_url"
    )?;

    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        {
            let fields = line.trim_end().split('\t').collect::<Vec<_>>();
            let get = |name: &str| -> &str {
                columns
                    .get(name)
                    .and_then(|index| fields.get(*index))
                    .copied()
                    .unwrap_or_default()
            };
            if get("Assembly") == "GRCh38" && matches!(get("Chromosome"), "MT" | "M") {
                let position = get("PositionVCF");
                let reference = get("ReferenceAlleleVCF");
                let alternate = get("AlternateAlleleVCF");
                if !position.is_empty()
                    && reference.len() == 1
                    && alternate.len() == 1
                    && reference != "na"
                    && alternate != "na"
                {
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}\t\t\t\t{}\t{}\t{}\tClinVar\t\t\t\t\t{}\t",
                        clean_tsv(position),
                        clean_tsv(reference),
                        clean_tsv(alternate),
                        clean_tsv(get("GeneSymbol")),
                        clean_tsv(get("PhenotypeList")),
                        clean_tsv(get("ClinicalSignificance")),
                        clean_tsv(&get("RCVaccession").replace('|', ";")),
                        clean_tsv(get("AlleleID"))
                    )?;
                }
            }
        }
        line.clear();
    }

    let curated = fs::read_to_string(default_bundled_clinical_tsv())
        .context("failed to read bundled curated clinical annotations")?;
    for row in curated.lines().skip(1) {
        if !row.trim().is_empty() {
            writeln!(out, "{row}")?;
        }
    }

    let status = gzip.wait().context("failed to wait for gzip")?;
    if !status.success() {
        anyhow::bail!("gzip failed while reading {}", source_gz.display());
    }
    Ok(())
}

fn clean_tsv(value: &str) -> String {
    value.replace(['\t', '\n', '\r'], " ").trim().to_string()
}

#[derive(Debug)]
struct VcfVariant {
    chrom: String,
    position: u64,
    reference: String,
    alternate: String,
    support: u64,
    depth: u64,
    gene: Option<String>,
    clinical_significance: Option<String>,
    ci95_low: Option<f64>,
    ci95_high: Option<f64>,
}

impl VcfVariant {
    fn heteroplasmy(&self) -> f64 {
        if self.depth == 0 {
            0.0
        } else {
            self.support as f64 / self.depth as f64
        }
    }
}

fn render_vcf(json: &str) -> Result<String> {
    let data: Value = serde_json::from_str(json).context("analysis returned invalid JSON")?;
    validate_result_contract(&data)?;
    let sample = data
        .pointer("/metadata/sample")
        .and_then(Value::as_str)
        .unwrap_or("sample");
    let chrom = data
        .pointer("/metadata/reference_accession")
        .and_then(Value::as_str)
        .filter(|value| !value.is_empty() && *value != "custom")
        .unwrap_or("MT");
    let reference_length = data
        .pointer("/metadata/reference_length")
        .and_then(Value::as_u64)
        .context("result contract violation: /metadata/reference_length must be an integer")?;
    let mut variants = BTreeMap::<(u64, String, String), VcfVariant>::new();
    let aggregates = data
        .get("variants")
        .and_then(Value::as_array)
        .context("result contract violation: /variants must be an array")?;
    for variant in aggregates {
        let Some(position) = variant.get("position").and_then(Value::as_u64) else {
            continue;
        };
        let Some(reference) = variant.get("ref").and_then(Value::as_str) else {
            continue;
        };
        let Some(alternate) = variant.get("alt").and_then(Value::as_str) else {
            continue;
        };
        let support = variant
            .get("alt_depth")
            .and_then(Value::as_u64)
            .unwrap_or(0);
        let depth = variant
            .get("callable_depth")
            .and_then(Value::as_u64)
            .unwrap_or(0);
        if reference.len() != 1 || alternate.len() != 1 || depth == 0 {
            continue;
        }
        variants.insert(
            (position, reference.to_string(), alternate.to_string()),
            VcfVariant {
                chrom: chrom.to_string(),
                position,
                reference: reference.to_string(),
                alternate: alternate.to_string(),
                support,
                depth,
                gene: variant
                    .get("gene")
                    .and_then(Value::as_str)
                    .map(str::to_string),
                clinical_significance: variant
                    .pointer("/annotation/pathogenicity")
                    .and_then(Value::as_str)
                    .map(str::to_string),
                ci95_low: variant.get("ci95_low").and_then(Value::as_f64),
                ci95_high: variant.get("ci95_high").and_then(Value::as_f64),
            },
        );
    }

    let mut out = String::new();
    out.push_str("##fileformat=VCFv4.3\n");
    out.push_str("##source=MitoArchitect\n");
    out.push_str(&format!(
        "##contig=<ID={},length={}>\n",
        sanitize_vcf_token(chrom),
        reference_length
    ));
    out.push_str("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate read support\">\n");
    out.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Locus-callable A/C/G/T molecule depth after quality filters\">\n");
    out.push_str("##INFO=<ID=HF,Number=A,Type=Float,Description=\"Heteroplasmy fraction estimated from alternate depth / locus-callable depth\">\n");
    out.push_str("##INFO=<ID=HF_CI95,Number=2,Type=Float,Description=\"Wilson 95% confidence interval for heteroplasmy fraction\">\n");
    out.push_str(
        "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Annotated mitochondrial gene\">\n",
    );
    out.push_str("##INFO=<ID=CLNSIG,Number=1,Type=String,Description=\"Clinical significance from local annotation cache\">\n");
    out.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype placeholder for haploid mtDNA export\">\n");
    out.push_str("##FORMAT=<ID=HF,Number=1,Type=Float,Description=\"Heteroplasmy fraction\">\n");
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
    out.push_str(&sanitize_vcf_token(sample));
    out.push('\n');

    for variant in variants.values() {
        let mut info = vec![
            format!("AC={}", variant.support),
            format!("DP={}", variant.depth),
            format!("HF={:.6}", variant.heteroplasmy()),
        ];
        if let Some(gene) = &variant.gene {
            info.push(format!("GENE={}", sanitize_vcf_token(gene)));
        }
        if let Some(clinical_significance) = &variant.clinical_significance {
            info.push(format!(
                "CLNSIG={}",
                sanitize_vcf_token(clinical_significance)
            ));
        }
        if let (Some(low), Some(high)) = (variant.ci95_low, variant.ci95_high) {
            info.push(format!("HF_CI95={low:.6},{high:.6}"));
        }
        out.push_str(&format!(
            "{}\t{}\t.\t{}\t{}\t.\tPASS\t{}\tGT:HF\t0/1:{:.6}\n",
            sanitize_vcf_token(&variant.chrom),
            variant.position,
            sanitize_vcf_token(&variant.reference),
            sanitize_vcf_token(&variant.alternate),
            info.join(";"),
            variant.heteroplasmy()
        ));
    }

    Ok(out)
}

fn validate_result_contract(data: &Value) -> Result<()> {
    require_object(data, "/metadata")?;
    require_object(data, "/filter_stats")?;
    require_string(data, "/metadata/schema_version")?;
    require_string(data, "/metadata/engine_version")?;
    require_u64(data, "/metadata/reference_length")?;
    require_array(data, "/reads")?;
    require_array(data, "/variants")?;
    require_array(data, "/svs")?;
    require_array(data, "/clusters")?;
    require_array(data, "/metadata/resources")?;
    require_object(data, "/filter_stats/numt_assessment")?;
    require_u64(data, "/filter_stats/passed_reads")?;
    Ok(())
}

fn require_object<'a>(
    data: &'a Value,
    pointer: &str,
) -> Result<&'a serde_json::Map<String, Value>> {
    data.pointer(pointer)
        .and_then(Value::as_object)
        .with_context(|| format!("result contract violation: {pointer} must be an object"))
}

fn require_array<'a>(data: &'a Value, pointer: &str) -> Result<&'a Vec<Value>> {
    data.pointer(pointer)
        .and_then(Value::as_array)
        .with_context(|| format!("result contract violation: {pointer} must be an array"))
}

fn require_string<'a>(data: &'a Value, pointer: &str) -> Result<&'a str> {
    data.pointer(pointer)
        .and_then(Value::as_str)
        .with_context(|| format!("result contract violation: {pointer} must be a string"))
}

fn require_u64(data: &Value, pointer: &str) -> Result<u64> {
    data.pointer(pointer)
        .and_then(Value::as_u64)
        .with_context(|| {
            format!("result contract violation: {pointer} must be an unsigned integer")
        })
}

fn sanitize_vcf_token(value: &str) -> String {
    value
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.' | ':' | '>') {
                c
            } else {
                '_'
            }
        })
        .collect()
}

fn default_bundled_clinical_tsv() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../core/data/clinical_annotations.tsv")
        .components()
        .collect::<PathBuf>()
}

fn default_cache_dir() -> PathBuf {
    if let Some(path) = std::env::var_os("XDG_CACHE_HOME") {
        return PathBuf::from(path);
    }
    std::env::var_os("HOME")
        .map(PathBuf::from)
        .map(|home| home.join(".cache"))
        .unwrap_or_else(|| Path::new(".").join(".cache"))
}

#[cfg(test)]
mod tests {
    use super::{assert_aux_present, assert_mapq_present, assert_sv_present, render_vcf};

    #[test]
    fn renders_vcf_with_locus_callable_heteroplasmy() {
        let json = r#"{
          "metadata": {
            "schema_version": "0.4",
            "engine_version": "0.4.0",
            "reference_length": 16569,
            "sample": "sample A",
            "reference_accession": "NC_012920.1",
            "resources": []
          },
          "filter_stats": {
            "passed_reads": 40,
            "numt_assessment": { "mode": "competitive_alignment" }
          },
          "variants": [
            {
              "position": 3243,
              "ref": "A",
              "alt": "G",
              "alt_depth": 2,
              "callable_depth": 3,
              "ci95_low": 0.207660,
              "ci95_high": 0.938508,
              "gene": "MT-TL1",
              "annotation": { "pathogenicity": "pathogenic" }
            }
          ],
          "reads": [],
          "svs": [],
          "clusters": []
        }"#;

        let vcf = render_vcf(json).expect("VCF should render");
        assert!(vcf.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A"));
        assert!(vcf.contains("NC_012920.1\t3243\t.\tA\tG\t.\tPASS\tAC=2;DP=3;HF=0.666667;GENE=MT-TL1;CLNSIG=pathogenic;HF_CI95=0.207660,0.938508\tGT:HF\t0/1:0.666667"));
    }

    #[test]
    fn validates_fixture_fields_beyond_snp() {
        let json: serde_json::Value = serde_json::from_str(
            r#"{
              "coverage_metrics": {
                "mapping_quality_histogram": [
                  { "mapq": 42, "count": 1 }
                ]
              },
              "svs": [
                { "id": "deletion:110-121" }
              ],
              "reads": [
                { "aux_tags": { "NM": "1", "MD": "2T7" } }
              ]
            }"#,
        )
        .expect("fixture JSON should parse");

        assert_sv_present(&json, "deletion:110-121").expect("SV should validate");
        assert_mapq_present(&json, 42).expect("MAPQ should validate");
        assert_aux_present(&json, "NM=1").expect("aux tag should validate");
        assert_aux_present(&json, "MD=2T7").expect("aux tag should validate");
    }
}
