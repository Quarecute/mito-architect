mod report;

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use indicatif::{ProgressBar, ProgressStyle};
use mito_ffi::{AnalyzeOptions, MitoEngine};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command as ProcessCommand, Stdio};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

const CLINICAL_TSV_HEADER: &str = "position\tref\talt\tgene\tconsequence\tprotein\tresidue\tstructure_id\tstructure_chain\tstructure_residue\tstructure_complex\tsource\tassertion_id\tallele_id\tdisease\tclinical_significance\treview_status\tassertion_date\tsource_url\treferences\tresource_version\tretrieved_at";

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
    CompareMtDnaServer2(CompareMtDnaServer2Args),
    Doctor,
    UpdateClinical(UpdateClinicalArgs),
    ClinicalSnapshot(ClinicalSnapshotArgs),
    ValidateFixture(ValidateFixtureArgs),
    ValidateClinicalFixture(ValidateClinicalFixtureArgs),
    ValidateClinicalManifest(ValidateClinicalManifestArgs),
    ValidateEvidenceFixture(ValidateEvidenceFixtureArgs),
    ValidateEvidenceGraphFixture(ValidateEvidenceGraphFixtureArgs),
    ValidateErrorManifest(ValidateErrorManifestArgs),
    ValidateHaplogroupManifest(ValidateHaplogroupManifestArgs),
    ValidateSvFixture(ValidateSvFixtureArgs),
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum ComparatorFilterPolicy {
    PassOnly,
    All,
}

#[derive(Debug, Args)]
struct CompareMtDnaServer2Args {
    /// Authoritative Mito-Architect schema 0.6 JSON result.
    #[arg(long = "result")]
    result: PathBuf,

    /// mtDNA-Server 2 variants.annotated.txt file.
    #[arg(long = "comparator")]
    comparator: PathBuf,

    /// Exact mtDNA-Server 2 sample ID. Required for multi-sample files.
    #[arg(long = "sample")]
    sample: Option<String>,

    /// Pinned upstream release whose variants.annotated.txt format is parsed.
    #[arg(long = "comparator-version", default_value = "2.1.16")]
    comparator_version: String,

    /// Include only PASS comparator records or every filter state.
    #[arg(long = "filter-policy", value_enum, default_value_t = ComparatorFilterPolicy::PassOnly)]
    filter_policy: ComparatorFilterPolicy,

    /// Deterministic machine-readable differential report.
    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    /// Optional development gate; comparator concordance is not truth accuracy.
    #[arg(long = "min-call-concordance")]
    min_call_concordance: Option<f64>,

    /// Optional development gate over variants called by both tools.
    #[arg(long = "max-mean-hf-delta")]
    max_mean_hf_delta: Option<f64>,
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

    /// Write the unified schema 0.6 SNV/small-indel projection as TSV.
    #[arg(long = "tsv")]
    tsv: Option<PathBuf>,

    /// Write bgzip-compressed VCF and a tabix index (requires bgzip and tabix).
    #[arg(long = "bgzip-vcf")]
    bgzip_vcf: Option<PathBuf>,

    /// Write a deterministic per-analysis provenance manifest.
    #[arg(long = "provenance-manifest")]
    provenance_manifest: Option<PathBuf>,

    /// Optional directory for deterministic schema 0.6 observation-page sidecars.
    #[arg(long = "evidence-pages-dir")]
    evidence_pages_dir: Option<PathBuf>,

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

    /// Emit the opt-in schema 0.6 fragment/molecule/event evidence graph.
    #[arg(long = "evidence-graph", default_value_t = false)]
    emit_evidence_graph: bool,

    /// Hard cap for in-memory schema 0.6 sparse observations.
    #[arg(long = "max-evidence-observations", default_value_t = 5_000_000)]
    max_evidence_observations: usize,

    /// Hard cap for callable-aware event-pair phase projections.
    #[arg(long = "max-phase-links", default_value_t = 1_000_000)]
    max_phase_links: usize,

    /// Maximum observations per schema 0.6 columnar page.
    #[arg(long = "evidence-page-size", default_value_t = 4096)]
    evidence_page_size: usize,

    /// Explicit SAM tag used as the physical-molecule identifier (for example MI).
    #[arg(long = "molecule-id-tag", default_value = "")]
    molecule_id_tag: String,

    /// Optional SAM tag carrying UMI metadata (for example RX).
    #[arg(long = "umi-tag", default_value = "")]
    umi_tag: String,

    /// Optional SAM tag carrying duplex metadata.
    #[arg(long = "duplex-tag", default_value = "")]
    duplex_tag: String,

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

#[derive(Clone, Debug, ValueEnum)]
enum ClinicalSnapshotAction {
    Stage,
    Activate,
    Rollback,
    Verify,
    Status,
}

#[derive(Debug, Args)]
struct ClinicalSnapshotArgs {
    #[arg(value_enum)]
    action: ClinicalSnapshotAction,

    #[arg(long = "store")]
    store: Option<PathBuf>,

    #[arg(long = "source")]
    source: Option<PathBuf>,

    #[arg(long = "snapshot-id")]
    snapshot_id: Option<String>,

    #[arg(long = "license-id")]
    license_id: Option<String>,

    #[arg(long = "source-policy")]
    source_policy: Option<String>,

    #[arg(long = "activate")]
    activate: bool,

    #[arg(long = "max-age-days")]
    max_age_days: Option<u64>,
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

#[derive(Debug, Args)]
struct ValidateEvidenceFixtureArgs {
    #[arg(long = "input")]
    input: PathBuf,

    #[arg(long = "expected-json")]
    expected_json: PathBuf,

    /// Write the deterministic projection before comparison. Intended for
    /// reviewed fixture regeneration after an explicit contract change.
    #[arg(long = "write-projection")]
    write_projection: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct ValidateEvidenceGraphFixtureArgs {
    #[arg(long = "input")]
    input: PathBuf,

    #[arg(long = "expected-json")]
    expected_json: PathBuf,

    #[arg(long = "evidence-page-size", default_value_t = 3)]
    evidence_page_size: usize,

    #[arg(long = "molecule-id-tag", default_value = "")]
    molecule_id_tag: String,

    #[arg(long = "umi-tag", default_value = "")]
    umi_tag: String,

    #[arg(long = "duplex-tag", default_value = "")]
    duplex_tag: String,

    /// Write the deterministic projection before comparison. Intended for
    /// reviewed fixture regeneration, not normal analysis output.
    #[arg(long = "write-projection")]
    write_projection: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct ValidateClinicalFixtureArgs {
    #[arg(long = "input")]
    input: PathBuf,

    #[arg(long = "annotations")]
    annotations: PathBuf,

    #[arg(long = "expected-json")]
    expected_json: PathBuf,
}

#[derive(Debug, Args)]
struct ValidateClinicalManifestArgs {
    #[arg(long = "manifest")]
    manifest: PathBuf,
}

#[derive(Debug, Args)]
struct ValidateErrorManifestArgs {
    #[arg(long = "manifest")]
    manifest: PathBuf,
}

#[derive(Debug, Args)]
struct ValidateHaplogroupManifestArgs {
    #[arg(long = "manifest")]
    manifest: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Analyze(args) => analyze(args),
        Command::CompareMtDnaServer2(args) => compare_mtdna_server2(args),
        Command::Doctor => doctor(),
        Command::UpdateClinical(args) => update_clinical(args).map(|path| {
            println!("clinical annotation cache: {}", path.display());
            println!("set MITO_CLINICAL_ANNOTATIONS={} to use it", path.display());
        }),
        Command::ClinicalSnapshot(args) => manage_clinical_snapshot(args),
        Command::ValidateFixture(args) => validate_fixture(args),
        Command::ValidateClinicalFixture(args) => validate_clinical_fixture(args),
        Command::ValidateClinicalManifest(args) => validate_clinical_manifest(args),
        Command::ValidateEvidenceFixture(args) => validate_evidence_fixture(args),
        Command::ValidateEvidenceGraphFixture(args) => validate_evidence_graph_fixture(args),
        Command::ValidateErrorManifest(args) => validate_error_manifest(args),
        Command::ValidateHaplogroupManifest(args) => validate_haplogroup_manifest(args),
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
        ("bgzip", &["--version"][..], false),
        ("tabix", &["--version"][..], false),
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
    let emit_evidence_graph = args.emit_evidence_graph
        || args.tsv.is_some()
        || args.bgzip_vcf.is_some()
        || args.evidence_pages_dir.is_some();
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
                emit_evidence_graph,
                max_evidence_observations: args.max_evidence_observations,
                max_phase_links: args.max_phase_links,
                evidence_page_size: args.evidence_page_size,
                molecule_id_tag: args.molecule_id_tag.clone(),
                umi_tag: args.umi_tag.clone(),
                duplex_tag: args.duplex_tag.clone(),
            },
        )
        .context("analysis failed")?;

    let mut exported = Vec::<ExportArtifact>::new();

    if let Some(directory) = &args.evidence_pages_dir {
        progress.set_message("writing evidence page resources");
        export_evidence_pages(&result, directory)?;
        let manifest = directory.join("manifest.json");
        exported.push(ExportArtifact::from_path("evidence_manifest", &manifest)?);
    }

    let vcf = if args.vcf.is_some() || args.bgzip_vcf.is_some() {
        progress.set_message("writing VCF export");
        Some(render_vcf(&result)?)
    } else {
        None
    };
    if let (Some(vcf_path), Some(vcf)) = (&args.vcf, vcf.as_deref()) {
        write_bytes_atomic(vcf_path, vcf.as_bytes())?;
        exported.push(ExportArtifact::from_path("vcf", vcf_path)?);
    }
    if let (Some(vcf_path), Some(vcf)) = (&args.bgzip_vcf, vcf.as_deref()) {
        progress.set_message("writing bgzip VCF and tabix index");
        let index_path = write_bgzip_vcf_with_tabix(vcf_path, vcf)?;
        exported.push(ExportArtifact::from_path("vcf_bgzip", vcf_path)?);
        exported.push(ExportArtifact::from_path("vcf_tabix", &index_path)?);
    }

    if let Some(tsv_path) = &args.tsv {
        progress.set_message("writing unified variant TSV");
        let tsv = render_variant_tsv(&result)?;
        write_bytes_atomic(tsv_path, tsv.as_bytes())?;
        exported.push(ExportArtifact::from_path("variant_tsv", tsv_path)?);
    }

    if !args.json {
        progress.set_message("rendering standalone report");
        let html = report::render_report(&result)?;
        write_bytes_atomic(&args.output, html.as_bytes())?;
        exported.push(ExportArtifact::from_path("html", &args.output)?);
    }

    if let Some(manifest_path) = &args.provenance_manifest {
        progress.set_message("writing provenance manifest");
        let manifest = build_provenance_manifest(&result, &args, &exported)?;
        write_json_atomic(manifest_path, &manifest)?;
    }

    if args.json {
        progress.finish_and_clear();
        println!("{result}");
    } else {
        progress.finish_with_message(format!("wrote {}", args.output.display()));
    }

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
    validate_clinical_tsv_schema(&source)?;

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
    let output_name = output
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("clinical_annotations.tsv");
    let temporary_output =
        output.with_file_name(format!(".{output_name}.tmp-{}", std::process::id()));
    fs::copy(&source, &temporary_output).with_context(|| {
        format!(
            "failed to copy clinical annotations from {} to {}",
            source.display(),
            temporary_output.display()
        )
    })?;
    if let Err(error) = fs::rename(&temporary_output, &output) {
        let _ = fs::remove_file(&temporary_output);
        return Err(error).with_context(|| {
            format!(
                "failed to atomically replace clinical cache {}",
                output.display()
            )
        });
    }
    write_clinical_cache_metadata(&output, &source, false)?;
    Ok(output)
}

fn validate_clinical_tsv_schema(path: &Path) -> Result<()> {
    let input = fs::File::open(path)
        .with_context(|| format!("failed to open clinical TSV {}", path.display()))?;
    let mut lines = BufReader::new(input).lines();
    let header = lines.next().transpose()?.context("clinical TSV is empty")?;
    if header.trim_end_matches('\r') != CLINICAL_TSV_HEADER {
        anyhow::bail!("clinical TSV header does not match assertion schema 1.0");
    }

    let mut row_count = 0usize;
    for (index, line) in lines.enumerate() {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 22 {
            anyhow::bail!(
                "clinical TSV row {} has {} fields; expected 22",
                index + 2,
                fields.len()
            );
        }
        let position = fields[0]
            .parse::<u64>()
            .with_context(|| format!("invalid clinical position at row {}", index + 2))?;
        let valid_allele = |value: &str| {
            value.len() == 1
                && matches!(
                    value.as_bytes()[0].to_ascii_uppercase(),
                    b'A' | b'C' | b'G' | b'T'
                )
        };
        if !(1..=16_569).contains(&position)
            || !valid_allele(fields[1])
            || !valid_allele(fields[2])
            || fields[1].eq_ignore_ascii_case(fields[2])
            || fields[11].trim().is_empty()
        {
            anyhow::bail!(
                "invalid clinical variant/assertion key at row {}",
                index + 2
            );
        }
        if !fields[18].is_empty()
            && !fields[18].starts_with("https://")
            && !fields[18].starts_with("http://")
        {
            anyhow::bail!("clinical source URL must be HTTP(S) at row {}", index + 2);
        }
        row_count += 1;
    }
    if row_count == 0 {
        anyhow::bail!("clinical TSV contains no assertions");
    }
    Ok(())
}

fn manage_clinical_snapshot(args: ClinicalSnapshotArgs) -> Result<()> {
    let store = args.store.unwrap_or_else(default_clinical_snapshot_store);
    let result = match args.action {
        ClinicalSnapshotAction::Stage => stage_clinical_snapshot(
            &store,
            args.source.as_deref(),
            args.snapshot_id.as_deref(),
            args.license_id.as_deref(),
            args.source_policy.as_deref(),
            args.activate,
        )?,
        ClinicalSnapshotAction::Activate => {
            let snapshot_id = args
                .snapshot_id
                .as_deref()
                .context("clinical snapshot activate requires --snapshot-id")?;
            activate_clinical_snapshot(&store, snapshot_id, args.max_age_days)?
        }
        ClinicalSnapshotAction::Rollback => rollback_clinical_snapshot(&store, args.max_age_days)?,
        ClinicalSnapshotAction::Verify => {
            let snapshot_id = match args.snapshot_id.as_deref() {
                Some(snapshot_id) => snapshot_id.to_owned(),
                None => active_snapshot_id(&store)?.context(
                    "clinical snapshot verify requires --snapshot-id or an active snapshot",
                )?,
            };
            verify_clinical_snapshot(&store, &snapshot_id, args.max_age_days)?
        }
        ClinicalSnapshotAction::Status => clinical_snapshot_status(&store, args.max_age_days)?,
    };
    println!("{}", serde_json::to_string_pretty(&result)?);
    Ok(())
}

fn default_clinical_snapshot_store() -> PathBuf {
    default_cache_dir()
        .join("mito-architect")
        .join("clinical-snapshots")
}

fn validate_snapshot_id(snapshot_id: &str) -> Result<()> {
    if snapshot_id.is_empty()
        || snapshot_id.len() > 128
        || snapshot_id.starts_with('.')
        || snapshot_id.ends_with('.')
        || !snapshot_id
            .chars()
            .all(|value| value.is_ascii_alphanumeric() || matches!(value, '.' | '_' | '-'))
    {
        anyhow::bail!(
            "snapshot ID must be 1-128 ASCII letters, digits, '.', '_', or '-', without edge dots"
        );
    }
    Ok(())
}

fn require_governance_label<'a>(name: &str, value: Option<&'a str>) -> Result<&'a str> {
    let value = value
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .with_context(|| format!("clinical snapshot stage requires --{name}"))?;
    if value.len() > 256
        || matches!(
            value.to_ascii_lowercase().as_str(),
            "unknown" | "unrecorded" | "not-recorded" | "none"
        )
        || value.chars().any(char::is_control)
    {
        anyhow::bail!("clinical snapshot --{name} is not a valid governed value");
    }
    Ok(value)
}

fn snapshot_dir(store: &Path, snapshot_id: &str) -> PathBuf {
    store.join("snapshots").join(snapshot_id)
}

fn snapshot_data_path(store: &Path, snapshot_id: &str) -> PathBuf {
    snapshot_dir(store, snapshot_id).join("clinical_annotations.tsv")
}

fn snapshot_manifest_path(store: &Path, snapshot_id: &str) -> PathBuf {
    snapshot_dir(store, snapshot_id).join("manifest.json")
}

fn clinical_tsv_summary(path: &Path) -> Result<Value> {
    validate_clinical_tsv_schema(path)?;
    let input = fs::File::open(path)
        .with_context(|| format!("failed to open clinical TSV {}", path.display()))?;
    let mut rows = 0u64;
    let mut variants = BTreeSet::new();
    let mut sources = BTreeSet::new();
    for (index, line) in BufReader::new(input).lines().enumerate() {
        let line = line?;
        if index == 0 || line.trim().is_empty() {
            continue;
        }
        let fields = line.trim_end_matches('\r').split('\t').collect::<Vec<_>>();
        rows = rows
            .checked_add(1)
            .context("clinical assertion count overflow")?;
        variants.insert(format!("{}:{}:{}", fields[0], fields[1], fields[2]));
        sources.insert(fields[11].to_owned());
    }
    Ok(serde_json::json!({
        "assertion_count": rows,
        "variant_count": variants.len(),
        "sources": sources,
    }))
}

fn stage_clinical_snapshot(
    store: &Path,
    source: Option<&Path>,
    snapshot_id: Option<&str>,
    license_id: Option<&str>,
    source_policy: Option<&str>,
    activate: bool,
) -> Result<Value> {
    let source = source.context("clinical snapshot stage requires --source")?;
    if !source.is_file() {
        anyhow::bail!(
            "clinical snapshot source does not exist: {}",
            source.display()
        );
    }
    let snapshot_id = snapshot_id.context("clinical snapshot stage requires --snapshot-id")?;
    validate_snapshot_id(snapshot_id)?;
    let license_id = require_governance_label("license-id", license_id)?;
    let source_policy = require_governance_label("source-policy", source_policy)?;
    let summary = clinical_tsv_summary(source)?;

    let snapshots = store.join("snapshots");
    fs::create_dir_all(&snapshots)
        .with_context(|| format!("failed to create {}", snapshots.display()))?;
    let destination = snapshot_dir(store, snapshot_id);
    if destination.exists() {
        anyhow::bail!("clinical snapshot already exists: {snapshot_id}");
    }
    let temporary = snapshots.join(format!(
        ".{snapshot_id}.tmp-{}-{}",
        std::process::id(),
        now_unix_secs()?
    ));
    if temporary.exists() {
        anyhow::bail!(
            "clinical snapshot staging path already exists: {}",
            temporary.display()
        );
    }
    fs::create_dir(&temporary)
        .with_context(|| format!("failed to create {}", temporary.display()))?;

    let stage_result = (|| -> Result<Value> {
        let staged_data = temporary.join("clinical_annotations.tsv");
        fs::copy(source, &staged_data).with_context(|| {
            format!(
                "failed to stage clinical snapshot from {} to {}",
                source.display(),
                staged_data.display()
            )
        })?;
        let staged_summary = clinical_tsv_summary(&staged_data)?;
        if staged_summary != summary {
            anyhow::bail!("clinical snapshot summary changed during staging");
        }
        let source_sha256 = sha256_file(source)?;
        let data_sha256 = sha256_file(&staged_data)?;
        if source_sha256 != data_sha256 {
            anyhow::bail!("clinical snapshot bytes changed during staging");
        }
        let staged_at_unix = now_unix_secs()?;
        let manifest = serde_json::json!({
            "schema_version": "1.0",
            "clinical_annotation_schema_version": "1.0",
            "snapshot_id": snapshot_id,
            "status": "staged",
            "source": source.display().to_string(),
            "source_sha256": source_sha256,
            "data_sha256": data_sha256,
            "staged_at_unix": staged_at_unix,
            "license_id": license_id,
            "source_policy": source_policy,
            "summary": summary,
        });
        fs::write(
            temporary.join("manifest.json"),
            serde_json::to_vec_pretty(&manifest)?,
        )
        .with_context(|| format!("failed to write manifest for snapshot {snapshot_id}"))?;
        sync_file(&staged_data)?;
        sync_file(&temporary.join("manifest.json"))?;
        sync_directory(&temporary)?;
        fs::rename(&temporary, &destination).with_context(|| {
            format!(
                "failed to atomically publish clinical snapshot {}",
                destination.display()
            )
        })?;
        sync_directory(&snapshots)?;
        Ok(manifest)
    })();
    if stage_result.is_err() {
        let _ = fs::remove_dir_all(&temporary);
    }
    let manifest = stage_result?;

    if activate {
        return activate_clinical_snapshot(store, snapshot_id, None);
    }
    Ok(serde_json::json!({
        "operation": "stage",
        "store": store,
        "snapshot": manifest,
        "data_path": snapshot_data_path(store, snapshot_id),
    }))
}

fn now_unix_secs() -> Result<u64> {
    Ok(SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system time is before the Unix epoch")?
        .as_secs())
}

fn read_json_file(path: &Path, label: &str) -> Result<Value> {
    serde_json::from_str(
        &fs::read_to_string(path)
            .with_context(|| format!("failed to read {label} {}", path.display()))?,
    )
    .with_context(|| format!("invalid JSON in {label} {}", path.display()))
}

fn read_snapshot_manifest(store: &Path, snapshot_id: &str) -> Result<Value> {
    validate_snapshot_id(snapshot_id)?;
    let path = snapshot_manifest_path(store, snapshot_id);
    let manifest = read_json_file(&path, "clinical snapshot manifest")?;
    if manifest.get("schema_version").and_then(Value::as_str) != Some("1.0")
        || manifest
            .get("clinical_annotation_schema_version")
            .and_then(Value::as_str)
            != Some("1.0")
        || manifest.get("snapshot_id").and_then(Value::as_str) != Some(snapshot_id)
        || manifest.get("status").and_then(Value::as_str) != Some("staged")
    {
        anyhow::bail!("clinical snapshot manifest contract mismatch: {snapshot_id}");
    }
    for field in ["source_sha256", "data_sha256"] {
        let digest = manifest
            .get(field)
            .and_then(Value::as_str)
            .with_context(|| format!("clinical snapshot manifest requires {field}"))?;
        if digest.len() != 64 || !digest.chars().all(|value| value.is_ascii_hexdigit()) {
            anyhow::bail!("clinical snapshot manifest has invalid {field}");
        }
    }
    for field in ["license_id", "source_policy"] {
        require_governance_label(
            &field.replace('_', "-"),
            manifest.get(field).and_then(Value::as_str),
        )?;
    }
    manifest
        .get("staged_at_unix")
        .and_then(Value::as_u64)
        .context("clinical snapshot manifest requires staged_at_unix")?;
    manifest
        .get("summary")
        .and_then(Value::as_object)
        .context("clinical snapshot manifest requires summary")?;
    Ok(manifest)
}

fn verify_clinical_snapshot(
    store: &Path,
    snapshot_id: &str,
    max_age_days: Option<u64>,
) -> Result<Value> {
    let manifest = read_snapshot_manifest(store, snapshot_id)?;
    let data_path = snapshot_data_path(store, snapshot_id);
    let summary = clinical_tsv_summary(&data_path)?;
    if manifest.get("summary") != Some(&summary) {
        anyhow::bail!("clinical snapshot summary mismatch: {snapshot_id}");
    }
    let actual_sha256 = sha256_file(&data_path)?;
    let expected_sha256 = manifest
        .get("data_sha256")
        .and_then(Value::as_str)
        .context("clinical snapshot manifest requires data_sha256")?;
    if actual_sha256 != expected_sha256 {
        anyhow::bail!(
            "clinical snapshot checksum mismatch for {snapshot_id}: expected {expected_sha256}, got {actual_sha256}"
        );
    }
    let staged_at = manifest
        .get("staged_at_unix")
        .and_then(Value::as_u64)
        .context("clinical snapshot manifest requires staged_at_unix")?;
    let now = now_unix_secs()?;
    if staged_at > now.saturating_add(300) {
        anyhow::bail!("clinical snapshot staging time is in the future: {snapshot_id}");
    }
    let age_seconds = now.saturating_sub(staged_at);
    if let Some(max_age_days) = max_age_days {
        let max_age_seconds = max_age_days
            .checked_mul(24 * 60 * 60)
            .context("clinical snapshot maximum age overflow")?;
        if age_seconds > max_age_seconds {
            anyhow::bail!(
                "clinical snapshot {snapshot_id} is stale: age {age_seconds}s exceeds {max_age_seconds}s"
            );
        }
    }
    Ok(serde_json::json!({
        "operation": "verify",
        "verified": true,
        "snapshot_id": snapshot_id,
        "data_path": data_path,
        "data_sha256": actual_sha256,
        "age_seconds": age_seconds,
        "max_age_days": max_age_days,
        "manifest": manifest,
    }))
}

fn snapshot_state_path(store: &Path) -> PathBuf {
    store.join("state.json")
}

fn read_snapshot_state(store: &Path) -> Result<Option<Value>> {
    let path = snapshot_state_path(store);
    if !path.exists() {
        return Ok(None);
    }
    if !path.is_file() {
        anyhow::bail!("clinical snapshot state is not a file: {}", path.display());
    }
    let state = read_json_file(&path, "clinical snapshot state")?;
    if state.get("schema_version").and_then(Value::as_str) != Some("1.0") {
        anyhow::bail!("clinical snapshot state requires schema_version 1.0");
    }
    if let Some(active) = state.get("active_snapshot").and_then(Value::as_str) {
        validate_snapshot_id(active)?;
    }
    if let Some(previous) = state.get("previous_snapshot").and_then(Value::as_str) {
        validate_snapshot_id(previous)?;
    }
    Ok(Some(state))
}

fn active_snapshot_id(store: &Path) -> Result<Option<String>> {
    Ok(read_snapshot_state(store)?.and_then(|state| {
        state
            .get("active_snapshot")
            .and_then(Value::as_str)
            .map(str::to_owned)
    }))
}

fn write_json_atomic(path: &Path, value: &Value) -> Result<()> {
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("state.json");
    let temporary = path.with_file_name(format!(
        ".{name}.tmp-{}-{}",
        std::process::id(),
        now_unix_secs()?
    ));
    fs::write(&temporary, serde_json::to_vec_pretty(value)?)
        .with_context(|| format!("failed to write {}", temporary.display()))?;
    sync_file(&temporary)?;
    if let Err(error) = fs::rename(&temporary, path) {
        let _ = fs::remove_file(&temporary);
        return Err(error)
            .with_context(|| format!("failed to atomically replace {}", path.display()));
    }
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        sync_directory(parent)?;
    }
    Ok(())
}

fn write_bytes_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("export");
    let temporary = path.with_file_name(format!(
        ".{name}.tmp-{}-{}",
        std::process::id(),
        now_unix_secs()?
    ));
    fs::write(&temporary, bytes)
        .with_context(|| format!("failed to write {}", temporary.display()))?;
    sync_file(&temporary)?;
    if let Err(error) = fs::rename(&temporary, path) {
        let _ = fs::remove_file(&temporary);
        return Err(error)
            .with_context(|| format!("failed to atomically replace {}", path.display()));
    }
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        sync_directory(parent)?;
    }
    Ok(())
}

fn path_with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut value = path.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

fn write_bgzip_vcf_with_tabix(path: &Path, vcf: &str) -> Result<PathBuf> {
    if path.extension().and_then(|value| value.to_str()) != Some("gz") {
        anyhow::bail!("--bgzip-vcf output must end in .gz: {}", path.display());
    }
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("variants.vcf.gz");
    let temporary = path.with_file_name(format!(
        ".{name}.tmp-{}-{}.gz",
        std::process::id(),
        now_unix_secs()?
    ));
    let temporary_index = path_with_suffix(&temporary, ".tbi");
    let final_index = path_with_suffix(path, ".tbi");

    let mut child = ProcessCommand::new("bgzip")
        .arg("-c")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .context("failed to start bgzip; install htslib/bgzip")?;
    child
        .stdin
        .as_mut()
        .context("bgzip stdin is unavailable")?
        .write_all(vcf.as_bytes())
        .context("failed to stream VCF into bgzip")?;
    let compressed = child
        .wait_with_output()
        .context("failed to wait for bgzip")?;
    if !compressed.status.success() {
        anyhow::bail!(
            "bgzip failed: {}",
            String::from_utf8_lossy(&compressed.stderr).trim()
        );
    }
    fs::write(&temporary, compressed.stdout)
        .with_context(|| format!("failed to write {}", temporary.display()))?;
    sync_file(&temporary)?;

    let tabix = ProcessCommand::new("tabix")
        .args(["-f", "-p", "vcf"])
        .arg(&temporary)
        .output()
        .context("failed to start tabix; install htslib/tabix")?;
    if !tabix.status.success() {
        let _ = fs::remove_file(&temporary);
        let _ = fs::remove_file(&temporary_index);
        anyhow::bail!(
            "tabix failed: {}",
            String::from_utf8_lossy(&tabix.stderr).trim()
        );
    }
    sync_file(&temporary_index)?;
    if let Err(error) = fs::rename(&temporary, path) {
        let _ = fs::remove_file(&temporary);
        let _ = fs::remove_file(&temporary_index);
        return Err(error).with_context(|| format!("failed to publish {}", path.display()));
    }
    if let Err(error) = fs::rename(&temporary_index, &final_index) {
        let _ = fs::remove_file(path);
        let _ = fs::remove_file(&temporary_index);
        return Err(error).with_context(|| format!("failed to publish {}", final_index.display()));
    }
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        sync_directory(parent)?;
    }
    Ok(final_index)
}

#[derive(Debug)]
struct ExportArtifact {
    kind: &'static str,
    path: PathBuf,
    bytes: u64,
    sha256: String,
}

impl ExportArtifact {
    fn from_path(kind: &'static str, path: &Path) -> Result<Self> {
        Ok(Self {
            kind,
            path: path.to_path_buf(),
            bytes: fs::metadata(path)
                .with_context(|| format!("failed to stat {}", path.display()))?
                .len(),
            sha256: sha256_file(path)?,
        })
    }
}

fn sha256_bytes(bytes: &[u8]) -> Result<String> {
    let mut child = ProcessCommand::new("sha256sum")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .context("failed to start sha256sum")?;
    child
        .stdin
        .as_mut()
        .context("sha256sum stdin is unavailable")?
        .write_all(bytes)
        .context("failed to stream bytes into sha256sum")?;
    let output = child
        .wait_with_output()
        .context("failed to wait for sha256sum")?;
    if !output.status.success() {
        anyhow::bail!("sha256sum failed while hashing analysis result");
    }
    String::from_utf8(output.stdout)
        .context("sha256sum returned non-UTF-8 output")?
        .split_whitespace()
        .next()
        .map(str::to_owned)
        .filter(|digest| digest.len() == 64 && digest.bytes().all(|byte| byte.is_ascii_hexdigit()))
        .context("sha256sum returned an invalid digest")
}

fn current_git_commit() -> Option<String> {
    let output = ProcessCommand::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let value = String::from_utf8(output.stdout).ok()?.trim().to_owned();
    (!value.is_empty()).then_some(value)
}

fn build_provenance_manifest(
    result: &str,
    args: &AnalyzeArgs,
    exported: &[ExportArtifact],
) -> Result<Value> {
    let data: Value = serde_json::from_str(result).context("analysis returned invalid JSON")?;
    validate_result_contract(&data)?;
    let reference_path = data
        .pointer("/metadata/reference_path")
        .and_then(Value::as_str)
        .context("result contract violation: /metadata/reference_path must be a string")?;
    let reference_path = Path::new(reference_path);
    let output_records = exported
        .iter()
        .map(|artifact| {
            serde_json::json!({
                "kind": artifact.kind,
                "path": artifact.path.display().to_string(),
                "bytes": artifact.bytes,
                "sha256": artifact.sha256,
            })
        })
        .collect::<Vec<_>>();
    Ok(serde_json::json!({
        "schema_version": "1.0",
        "timestamp_policy": "omitted_for_determinism",
        "software": {
            "name": "Mito-Architect",
            "engine_version": data.pointer("/metadata/engine_version"),
            "git_commit": current_git_commit(),
            "result_schema_version": data.pointer("/metadata/schema_version"),
            "sv_event_schema_version": data.pointer("/metadata/sv_event_schema_version"),
            "complex_sv_event_schema_version": data.pointer("/metadata/complex_sv_event_schema_version"),
            "clinical_annotation_schema_version": data.pointer("/metadata/clinical_annotation_schema_version"),
        },
        "determinism": {
            "deterministic_algorithms": true,
            "random_seed": Value::Null,
            "thread_count": data.pointer("/metadata/threads"),
        },
        "input": {
            "path": args.input.display().to_string(),
            "bytes": fs::metadata(&args.input)?.len(),
            "sha256": sha256_file(&args.input)?,
        },
        "reference": {
            "path": reference_path.display().to_string(),
            "accession": data.pointer("/metadata/reference_accession"),
            "length": data.pointer("/metadata/reference_length"),
            "sha256": sha256_file(reference_path)?,
        },
        "calling_parameters": data.pointer("/metadata/calling_parameters"),
        "resources": data.pointer("/metadata/resources"),
        "command_line": std::env::args().collect::<Vec<_>>(),
        "authoritative_result": {
            "format": "application/json",
            "bytes": result.len(),
            "sha256": sha256_bytes(result.as_bytes())?,
        },
        "exports": output_records,
    }))
}

fn export_evidence_pages(result: &str, directory: &Path) -> Result<()> {
    let data: Value = serde_json::from_str(result).context("analysis returned invalid JSON")?;
    if data
        .pointer("/metadata/schema_version")
        .and_then(Value::as_str)
        != Some("0.6")
    {
        anyhow::bail!("--evidence-pages-dir requires --evidence-graph");
    }
    validate_result_contract(&data)?;
    fs::create_dir_all(directory)
        .with_context(|| format!("failed to create {}", directory.display()))?;
    let pages = require_array(&data, "/observation_pages")?;
    let mut manifest_pages = Vec::with_capacity(pages.len());
    for (page_index, page) in pages.iter().enumerate() {
        let file_name = format!("observations-{page_index:06}.json");
        let path = directory.join(&file_name);
        write_json_atomic(&path, page)?;
        let bytes = fs::metadata(&path)
            .with_context(|| format!("failed to stat {}", path.display()))?
            .len();
        manifest_pages.push(serde_json::json!({
            "index": page_index,
            "offset": page.get("offset"),
            "count": page.get("count"),
            "path": file_name,
            "bytes": bytes,
            "sha256": sha256_file(&path)?,
        }));
    }
    let manifest = serde_json::json!({
        "schema_version": "1.0",
        "result_schema_version": "0.6",
        "sample": data.pointer("/metadata/sample"),
        "reference_accession": data.pointer("/metadata/reference_accession"),
        "reference_length": data.pointer("/metadata/reference_length"),
        "evidence_encoding": data.get("evidence_encoding"),
        "pages": manifest_pages,
    });
    write_json_atomic(&directory.join("manifest.json"), &manifest)
}

fn sync_file(path: &Path) -> Result<()> {
    fs::File::open(path)
        .with_context(|| format!("failed to open {} for synchronization", path.display()))?
        .sync_all()
        .with_context(|| format!("failed to synchronize {}", path.display()))
}

fn sync_directory(path: &Path) -> Result<()> {
    fs::File::open(path)
        .with_context(|| {
            format!(
                "failed to open directory {} for synchronization",
                path.display()
            )
        })?
        .sync_all()
        .with_context(|| format!("failed to synchronize directory {}", path.display()))
}

fn activate_clinical_snapshot(
    store: &Path,
    snapshot_id: &str,
    max_age_days: Option<u64>,
) -> Result<Value> {
    let verified = verify_clinical_snapshot(store, snapshot_id, max_age_days)?;
    let current_state = read_snapshot_state(store)?;
    let current = current_state
        .as_ref()
        .and_then(|state| state.get("active_snapshot"))
        .and_then(Value::as_str);
    let previous = current
        .filter(|current| *current != snapshot_id)
        .map(str::to_owned)
        .or_else(|| {
            current_state.as_ref().and_then(|state| {
                state
                    .get("previous_snapshot")
                    .and_then(Value::as_str)
                    .map(str::to_owned)
            })
        });
    let state = serde_json::json!({
        "schema_version": "1.0",
        "active_snapshot": snapshot_id,
        "previous_snapshot": previous,
        "updated_at_unix": now_unix_secs()?,
    });
    write_json_atomic(&snapshot_state_path(store), &state)?;
    Ok(serde_json::json!({
        "operation": "activate",
        "state": state,
        "active_path": snapshot_data_path(store, snapshot_id),
        "verification": verified,
    }))
}

fn rollback_clinical_snapshot(store: &Path, max_age_days: Option<u64>) -> Result<Value> {
    let state = read_snapshot_state(store)?.context("no clinical snapshot state to roll back")?;
    let active = state
        .get("active_snapshot")
        .and_then(Value::as_str)
        .context("clinical snapshot state has no active snapshot")?;
    let previous = state
        .get("previous_snapshot")
        .and_then(Value::as_str)
        .context("clinical snapshot state has no previous snapshot")?;
    let verified = verify_clinical_snapshot(store, previous, max_age_days)?;
    let next_state = serde_json::json!({
        "schema_version": "1.0",
        "active_snapshot": previous,
        "previous_snapshot": active,
        "updated_at_unix": now_unix_secs()?,
    });
    write_json_atomic(&snapshot_state_path(store), &next_state)?;
    Ok(serde_json::json!({
        "operation": "rollback",
        "state": next_state,
        "active_path": snapshot_data_path(store, previous),
        "verification": verified,
    }))
}

fn clinical_snapshot_status(store: &Path, max_age_days: Option<u64>) -> Result<Value> {
    let Some(state) = read_snapshot_state(store)? else {
        return Ok(serde_json::json!({
            "operation": "status",
            "store": store,
            "active_snapshot": null,
            "active_path": null,
            "verified": false,
        }));
    };
    let active = state
        .get("active_snapshot")
        .and_then(Value::as_str)
        .context("clinical snapshot state has no active snapshot")?;
    let verification = verify_clinical_snapshot(store, active, max_age_days)?;
    Ok(serde_json::json!({
        "operation": "status",
        "store": store,
        "active_snapshot": active,
        "active_path": snapshot_data_path(store, active),
        "state": state,
        "verified": true,
        "verification": verification,
    }))
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

#[derive(Clone, Debug)]
struct ComparatorVariant {
    sample_id: String,
    filter: String,
    position: u64,
    reference: String,
    alternate: String,
    heteroplasmy: f64,
    coverage: u64,
    coverage_forward: u64,
    coverage_reverse: u64,
    mutation: String,
}

#[derive(Clone, Debug)]
struct MitoVariantProjection {
    event_id: String,
    position: u64,
    reference: String,
    alternate: String,
    heteroplasmy: f64,
    callable_depth: u64,
    alt_depth: u64,
    filter_status: String,
}

fn compare_mtdna_server2(args: CompareMtDnaServer2Args) -> Result<()> {
    if args.comparator_version != "2.1.16" {
        anyhow::bail!(
            "unsupported mtDNA-Server 2 format version {}; only 2.1.16 is pinned",
            args.comparator_version
        );
    }
    validate_unit_interval(args.min_call_concordance, "--min-call-concordance")?;
    validate_unit_interval(args.max_mean_hf_delta, "--max-mean-hf-delta")?;
    if args
        .sample
        .as_deref()
        .is_some_and(|sample| sample.is_empty() || sample.trim() != sample)
    {
        anyhow::bail!("--sample must be a non-empty exact ID without surrounding whitespace");
    }

    let result_bytes = fs::read(&args.result)
        .with_context(|| format!("failed to read {}", args.result.display()))?;
    let result: Value = serde_json::from_slice(&result_bytes)
        .with_context(|| format!("{} is not valid JSON", args.result.display()))?;
    validate_result_contract(&result)
        .with_context(|| format!("{} is not a valid analysis result", args.result.display()))?;
    if result
        .pointer("/metadata/schema_version")
        .and_then(Value::as_str)
        != Some("0.6")
    {
        anyhow::bail!("mtDNA-Server 2 comparison requires a schema 0.6 result");
    }
    if result
        .pointer("/metadata/reference_length")
        .and_then(Value::as_u64)
        != Some(16_569)
    {
        anyhow::bail!(
            "mtDNA-Server 2 comparison requires the 16,569 bp human rCRS coordinate system"
        );
    }
    if let Some(accession) = result
        .pointer("/metadata/reference_accession")
        .and_then(Value::as_str)
    {
        if accession != "NC_012920.1" {
            anyhow::bail!(
                "mtDNA-Server 2 comparison requires rCRS accession NC_012920.1, found {accession}"
            );
        }
    }

    let comparator_rows = parse_mtdna_server2_variants(&args.comparator)?;
    let available_samples = comparator_rows
        .iter()
        .map(|row| row.sample_id.as_str())
        .collect::<BTreeSet<_>>();
    let selected_sample =
        match args.sample.as_deref() {
            Some(sample) if available_samples.is_empty() => sample.to_owned(),
            Some(sample) if available_samples.contains(sample) => sample.to_owned(),
            Some(sample) => anyhow::bail!(
                "sample {sample:?} is absent from {}; available samples: {}",
                args.comparator.display(),
                available_samples
                    .iter()
                    .copied()
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
            None if available_samples.len() == 1 => available_samples
                .iter()
                .next()
                .context("comparator sample set is unexpectedly empty")?
                .to_string(),
            None if available_samples.is_empty() => anyhow::bail!(
                "--sample is required when the comparator file contains only a header"
            ),
            None => {
                anyhow::bail!(
            "--sample is required for a multi-sample comparator file; available samples: {}",
            available_samples.iter().copied().collect::<Vec<_>>().join(", ")
        )
            }
        };

    let include_comparator = |row: &&ComparatorVariant| {
        row.sample_id == selected_sample
            && match args.filter_policy {
                ComparatorFilterPolicy::PassOnly => row.filter == "PASS",
                ComparatorFilterPolicy::All => true,
            }
    };
    let mut comparator_by_key = BTreeMap::<String, &ComparatorVariant>::new();
    for row in comparator_rows.iter().filter(include_comparator) {
        let key = variant_comparison_key(row.position, &row.reference, &row.alternate);
        if comparator_by_key.insert(key.clone(), row).is_some() {
            anyhow::bail!("duplicate mtDNA-Server 2 variant {key} for sample {selected_sample}");
        }
    }
    let mito_variants = parse_mito_variant_projection(&result)?;
    let mut mito_by_key = BTreeMap::<String, &MitoVariantProjection>::new();
    for variant in &mito_variants {
        let key = variant_comparison_key(variant.position, &variant.reference, &variant.alternate);
        if mito_by_key.insert(key.clone(), variant).is_some() {
            anyhow::bail!("duplicate Mito-Architect variant comparison key {key}");
        }
    }

    let mut matched = Vec::new();
    let mut mito_only = Vec::new();
    let mut comparator_only = Vec::new();
    let mut hf_delta_sum = 0.0;
    let mut max_hf_delta = 0.0_f64;
    for (key, mito) in &mito_by_key {
        if let Some(comparator) = comparator_by_key.get(key) {
            let hf_delta = (mito.heteroplasmy - comparator.heteroplasmy).abs();
            hf_delta_sum += hf_delta;
            max_hf_delta = max_hf_delta.max(hf_delta);
            matched.push(serde_json::json!({
                "key": key,
                "event_id": mito.event_id,
                "position": mito.position,
                "ref": mito.reference,
                "alt": mito.alternate,
                "mito_architect": {
                    "heteroplasmy": mito.heteroplasmy,
                    "callable_depth": mito.callable_depth,
                    "alt_depth": mito.alt_depth,
                    "filter_status": mito.filter_status,
                },
                "mtdna_server_2": {
                    "heteroplasmy": comparator.heteroplasmy,
                    "coverage": comparator.coverage,
                    "coverage_forward": comparator.coverage_forward,
                    "coverage_reverse": comparator.coverage_reverse,
                    "filter": comparator.filter,
                    "mutation": comparator.mutation,
                },
                "absolute_hf_delta": hf_delta,
            }));
        } else {
            mito_only.push(serde_json::json!({
                "key": key,
                "event_id": mito.event_id,
                "position": mito.position,
                "ref": mito.reference,
                "alt": mito.alternate,
                "heteroplasmy": mito.heteroplasmy,
                "callable_depth": mito.callable_depth,
                "alt_depth": mito.alt_depth,
                "filter_status": mito.filter_status,
            }));
        }
    }
    for (key, comparator) in &comparator_by_key {
        if !mito_by_key.contains_key(key) {
            comparator_only.push(serde_json::json!({
                "key": key,
                "position": comparator.position,
                "ref": comparator.reference,
                "alt": comparator.alternate,
                "heteroplasmy": comparator.heteroplasmy,
                "coverage": comparator.coverage,
                "coverage_forward": comparator.coverage_forward,
                "coverage_reverse": comparator.coverage_reverse,
                "filter": comparator.filter,
                "mutation": comparator.mutation,
            }));
        }
    }

    let union_count = matched.len() + mito_only.len() + comparator_only.len();
    let call_concordance = if union_count == 0 {
        1.0
    } else {
        matched.len() as f64 / union_count as f64
    };
    let mean_hf_delta = if matched.is_empty() {
        None
    } else {
        Some(hf_delta_sum / matched.len() as f64)
    };
    let filter_policy = match args.filter_policy {
        ComparatorFilterPolicy::PassOnly => "pass_only",
        ComparatorFilterPolicy::All => "all",
    };
    let report = serde_json::json!({
        "schema_version": "1.0",
        "report_type": "variant_callset_differential",
        "interpretation": "Comparator concordance is not analytical truth, sensitivity, specificity, or clinical validation.",
        "coordinate_system": {
            "name": "rCRS",
            "accession": "NC_012920.1",
            "length": 16569,
        },
        "mtdna_server_2": {
            "version": args.comparator_version,
            "format": "variants.annotated.txt",
            "format_contract": "mtdna-server-2-v2.1.16-variants-annotated",
            "upstream_release": "https://github.com/genepi/mtdna-server-2/releases/tag/v2.1.16",
            "sample_id": selected_sample,
            "filter_policy": filter_policy,
        },
        "provenance": {
            "mito_architect_result": {
                "path": args.result.to_string_lossy(),
                "sha256": sha256_bytes(&result_bytes)?,
            },
            "comparator_result": {
                "path": args.comparator.to_string_lossy(),
                "sha256": sha256_file(&args.comparator)?,
            },
        },
        "metrics": {
            "matched": matched.len(),
            "mito_architect_only": mito_only.len(),
            "mtdna_server_2_only": comparator_only.len(),
            "union": union_count,
            "call_concordance": call_concordance,
            "mean_absolute_hf_delta": mean_hf_delta,
            "max_absolute_hf_delta": if matched.is_empty() { None } else { Some(max_hf_delta) },
        },
        "matched": matched,
        "mito_architect_only": mito_only,
        "mtdna_server_2_only": comparator_only,
    });
    write_json_atomic(&args.output, &report)?;

    if let Some(minimum) = args.min_call_concordance {
        if call_concordance + f64::EPSILON < minimum {
            anyhow::bail!(
                "call concordance {call_concordance:.9} is below the development gate {minimum:.9}; report written to {}",
                args.output.display()
            );
        }
    }
    if let Some(maximum) = args.max_mean_hf_delta {
        match mean_hf_delta {
            Some(value) if value <= maximum + f64::EPSILON => {}
            Some(value) => anyhow::bail!(
                "mean absolute HF delta {value:.9} exceeds the development gate {maximum:.9}; report written to {}",
                args.output.display()
            ),
            None => anyhow::bail!(
                "mean absolute HF delta is undefined because no variants matched; report written to {}",
                args.output.display()
            ),
        }
    }
    println!(
        "mtDNA-Server 2 differential: matched={} mito_only={} comparator_only={} concordance={:.6}",
        report
            .pointer("/metrics/matched")
            .and_then(Value::as_u64)
            .unwrap_or(0),
        report
            .pointer("/metrics/mito_architect_only")
            .and_then(Value::as_u64)
            .unwrap_or(0),
        report
            .pointer("/metrics/mtdna_server_2_only")
            .and_then(Value::as_u64)
            .unwrap_or(0),
        call_concordance
    );
    Ok(())
}

fn validate_unit_interval(value: Option<f64>, name: &str) -> Result<()> {
    if value.is_some_and(|value| !value.is_finite() || !(0.0..=1.0).contains(&value)) {
        anyhow::bail!("{name} must be a finite number between 0 and 1");
    }
    Ok(())
}

fn parse_mtdna_server2_variants(path: &Path) -> Result<Vec<ComparatorVariant>> {
    const REQUIRED: [&str; 15] = [
        "ID",
        "Filter",
        "Pos",
        "Ref",
        "Variant",
        "VariantLevel",
        "MajorBase",
        "MajorLevel",
        "MinorBase",
        "MinorLevel",
        "Coverage",
        "CoverageFWD",
        "CoverageREV",
        "Type",
        "Mutation",
    ];
    let file =
        fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut lines = BufReader::new(file).lines();
    let header_line = lines
        .next()
        .transpose()
        .with_context(|| format!("failed to read {}", path.display()))?
        .context("mtDNA-Server 2 comparator file is empty")?;
    let header = header_line
        .trim_end_matches('\r')
        .split('\t')
        .collect::<Vec<_>>();
    let mut indices = BTreeMap::<&str, usize>::new();
    for (index, name) in header.iter().copied().enumerate() {
        if name.is_empty() || indices.insert(name, index).is_some() {
            anyhow::bail!("mtDNA-Server 2 header contains an empty or duplicate column");
        }
    }
    for required in REQUIRED {
        if !indices.contains_key(required) {
            anyhow::bail!(
                "mtDNA-Server 2 v2.1.16 variants.annotated.txt is missing required column {required}"
            );
        }
    }

    let field = |fields: &[&str], name: &str| -> Result<String> {
        let index = *indices
            .get(name)
            .with_context(|| format!("required comparator column {name} is unresolved"))?;
        fields
            .get(index)
            .map(|value| value.trim().to_owned())
            .with_context(|| format!("comparator row has no {name} column"))
    };
    let mut rows = Vec::new();
    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2;
        let line =
            line.with_context(|| format!("failed to read {}:{line_number}", path.display()))?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() {
            anyhow::bail!("empty comparator row at {}:{line_number}", path.display());
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != header.len() {
            anyhow::bail!(
                "{}:{line_number} has {} columns; expected {}",
                path.display(),
                fields.len(),
                header.len()
            );
        }
        let sample_id = field(&fields, "ID")?;
        let filter = field(&fields, "Filter")?;
        let position = parse_comparator_u64(&field(&fields, "Pos")?, "Pos", line_number)?;
        if !(1..=16_569).contains(&position) {
            anyhow::bail!("{}:{line_number} Pos is outside rCRS", path.display());
        }
        let reference = field(&fields, "Ref")?.to_ascii_uppercase();
        let alternate = field(&fields, "Variant")?.to_ascii_uppercase();
        if sample_id.is_empty()
            || filter.is_empty()
            || !valid_comparator_allele(&reference)
            || !valid_comparator_allele(&alternate)
        {
            anyhow::bail!(
                "{}:{line_number} has an invalid identity, filter, REF, or ALT",
                path.display()
            );
        }
        let heteroplasmy = parse_comparator_fraction(
            &field(&fields, "VariantLevel")?,
            "VariantLevel",
            line_number,
        )?;
        let major_base = field(&fields, "MajorBase")?.to_ascii_uppercase();
        let minor_base = field(&fields, "MinorBase")?.to_ascii_uppercase();
        if !valid_comparator_allele(&major_base) || !valid_comparator_allele(&minor_base) {
            anyhow::bail!(
                "{}:{line_number} has an invalid MajorBase or MinorBase",
                path.display()
            );
        }
        let _major_level =
            parse_comparator_fraction(&field(&fields, "MajorLevel")?, "MajorLevel", line_number)?;
        let _minor_level =
            parse_comparator_fraction(&field(&fields, "MinorLevel")?, "MinorLevel", line_number)?;
        let coverage = parse_comparator_u64(&field(&fields, "Coverage")?, "Coverage", line_number)?;
        let coverage_forward =
            parse_comparator_u64(&field(&fields, "CoverageFWD")?, "CoverageFWD", line_number)?;
        let coverage_reverse =
            parse_comparator_u64(&field(&fields, "CoverageREV")?, "CoverageREV", line_number)?;
        if coverage_forward.checked_add(coverage_reverse) != Some(coverage) {
            anyhow::bail!(
                "{}:{line_number} CoverageFWD + CoverageREV does not equal Coverage",
                path.display()
            );
        }
        if coverage == 0 {
            anyhow::bail!("{}:{line_number} variant Coverage is zero", path.display());
        }
        let variant_type = field(&fields, "Type")?;
        let mutation = field(&fields, "Mutation")?;
        if variant_type.is_empty() || mutation.is_empty() {
            anyhow::bail!(
                "{}:{line_number} has an empty Type or Mutation",
                path.display()
            );
        }
        rows.push(ComparatorVariant {
            sample_id,
            filter,
            position,
            reference,
            alternate,
            heteroplasmy,
            coverage,
            coverage_forward,
            coverage_reverse,
            mutation,
        });
    }
    Ok(rows)
}

fn parse_comparator_u64(value: &str, field: &str, line: usize) -> Result<u64> {
    value
        .parse::<u64>()
        .with_context(|| format!("comparator line {line} {field} is not an unsigned integer"))
}

fn parse_comparator_fraction(value: &str, field: &str, line: usize) -> Result<f64> {
    let value = value
        .parse::<f64>()
        .with_context(|| format!("comparator line {line} {field} is not numeric"))?;
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        anyhow::bail!("comparator line {line} {field} is outside [0,1]");
    }
    Ok(value)
}

fn valid_comparator_allele(value: &str) -> bool {
    !value.is_empty()
        && value
            .bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N' | b'-'))
}

fn parse_mito_variant_projection(result: &Value) -> Result<Vec<MitoVariantProjection>> {
    require_array(result, "/variants")?
        .iter()
        .enumerate()
        .map(|(index, variant)| {
            let string = |name: &str| {
                variant
                    .get(name)
                    .and_then(Value::as_str)
                    .map(str::to_owned)
                    .with_context(|| format!("/variants/{index}/{name} must be a string"))
            };
            let unsigned = |name: &str| {
                variant
                    .get(name)
                    .and_then(Value::as_u64)
                    .with_context(|| format!("/variants/{index}/{name} must be an integer"))
            };
            let heteroplasmy = variant
                .get("heteroplasmy")
                .and_then(Value::as_f64)
                .with_context(|| format!("/variants/{index}/heteroplasmy must be a number"))?;
            if !heteroplasmy.is_finite() || !(0.0..=1.0).contains(&heteroplasmy) {
                anyhow::bail!("/variants/{index}/heteroplasmy is outside [0,1]");
            }
            Ok(MitoVariantProjection {
                event_id: string("event_id")?,
                position: unsigned("position")?,
                reference: string("ref")?.to_ascii_uppercase(),
                alternate: string("alt")?.to_ascii_uppercase(),
                heteroplasmy,
                callable_depth: unsigned("callable_depth")?,
                alt_depth: unsigned("alt_depth")?,
                filter_status: string("filter_status")?,
            })
        })
        .collect()
}

fn variant_comparison_key(position: u64, reference: &str, alternate: &str) -> String {
    format!("{position}:{reference}:{alternate}")
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
        "complex_events": data.get("complex_events"),
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

fn validate_evidence_fixture(args: ValidateEvidenceFixtureArgs) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("fixture input does not exist: {}", args.input.display());
    }
    if !args.expected_json.exists() {
        anyhow::bail!(
            "expected evidence JSON does not exist: {}",
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
    validate_result_contract(&data)?;
    let actual = project_evidence(&data)?;
    if let Some(path) = &args.write_projection {
        write_json_atomic(path, &actual)
            .with_context(|| format!("failed to write projection {}", path.display()))?;
    }
    let expected: Value = serde_json::from_str(
        &fs::read_to_string(&args.expected_json)
            .with_context(|| format!("failed to read {}", args.expected_json.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.expected_json.display()))?;
    if actual != expected {
        anyhow::bail!(
            "scientific evidence golden mismatch for {}\nexpected:\n{}\nactual:\n{}",
            args.input.display(),
            serde_json::to_string_pretty(&expected)?,
            serde_json::to_string_pretty(&actual)?
        );
    }

    println!(
        "scientific evidence fixture validation passed: {}",
        args.input.display()
    );
    Ok(())
}

fn validate_evidence_graph_fixture(args: ValidateEvidenceGraphFixtureArgs) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("fixture input does not exist: {}", args.input.display());
    }
    if !args.expected_json.exists() {
        anyhow::bail!(
            "expected evidence-graph JSON does not exist: {}",
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
                threads: 4,
                emit_evidence_graph: true,
                evidence_page_size: args.evidence_page_size,
                molecule_id_tag: args.molecule_id_tag,
                umi_tag: args.umi_tag,
                duplex_tag: args.duplex_tag,
                ..AnalyzeOptions::default()
            },
        )
        .context("schema 0.6 analysis failed")?;
    let data: Value = serde_json::from_str(&json).context("analysis returned invalid JSON")?;
    validate_result_contract(&data)?;
    let actual = project_evidence_graph(&data)?;
    if let Some(path) = &args.write_projection {
        write_json_atomic(path, &actual)
            .with_context(|| format!("failed to write projection {}", path.display()))?;
    }
    let expected: Value = serde_json::from_str(
        &fs::read_to_string(&args.expected_json)
            .with_context(|| format!("failed to read {}", args.expected_json.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.expected_json.display()))?;
    if actual != expected {
        anyhow::bail!(
            "schema 0.6 evidence-graph golden mismatch for {}\nexpected:\n{}\nactual:\n{}",
            args.input.display(),
            serde_json::to_string_pretty(&expected)?,
            serde_json::to_string_pretty(&actual)?
        );
    }

    println!(
        "schema 0.6 evidence-graph fixture validation passed: {}",
        args.input.display()
    );
    Ok(())
}

struct ScopedEnvironmentVariable {
    name: &'static str,
    previous: Option<std::ffi::OsString>,
}

impl ScopedEnvironmentVariable {
    fn set(name: &'static str, value: &Path) -> Self {
        let previous = std::env::var_os(name);
        std::env::set_var(name, value);
        Self { name, previous }
    }
}

impl Drop for ScopedEnvironmentVariable {
    fn drop(&mut self) {
        if let Some(previous) = &self.previous {
            std::env::set_var(self.name, previous);
        } else {
            std::env::remove_var(self.name);
        }
    }
}

fn validate_clinical_fixture(args: ValidateClinicalFixtureArgs) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("fixture input does not exist: {}", args.input.display());
    }
    if !args.annotations.exists() {
        anyhow::bail!(
            "clinical fixture does not exist: {}",
            args.annotations.display()
        );
    }
    if !args.expected_json.exists() {
        anyhow::bail!(
            "expected clinical JSON does not exist: {}",
            args.expected_json.display()
        );
    }
    validate_clinical_tsv_schema(&args.annotations)?;
    let _environment =
        ScopedEnvironmentVariable::set("MITO_CLINICAL_ANNOTATIONS", &args.annotations);
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
    validate_result_contract(&data)?;
    let variants = require_array(&data, "/variants")?
        .iter()
        .filter_map(|variant| {
            variant.get("annotation").map(|annotation| {
                serde_json::json!({
                    "position": variant.get("position"),
                    "ref": variant.get("ref"),
                    "alt": variant.get("alt"),
                    "annotation": annotation,
                })
            })
        })
        .collect::<Vec<_>>();
    let actual = serde_json::json!({
        "schema_version": data.pointer("/metadata/clinical_annotation_schema_version"),
        "variants": variants,
    });
    let expected: Value = serde_json::from_str(
        &fs::read_to_string(&args.expected_json)
            .with_context(|| format!("failed to read {}", args.expected_json.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.expected_json.display()))?;
    if actual != expected {
        anyhow::bail!(
            "clinical assertion golden mismatch for {}\nexpected:\n{}\nactual:\n{}",
            args.input.display(),
            serde_json::to_string_pretty(&expected)?,
            serde_json::to_string_pretty(&actual)?
        );
    }
    println!(
        "clinical assertion fixture validation passed: {}",
        args.input.display()
    );
    Ok(())
}

fn validate_clinical_manifest(args: ValidateClinicalManifestArgs) -> Result<()> {
    if !args.manifest.exists() {
        anyhow::bail!(
            "clinical manifest does not exist: {}",
            args.manifest.display()
        );
    }
    let manifest: Value = serde_json::from_str(
        &fs::read_to_string(&args.manifest)
            .with_context(|| format!("failed to read {}", args.manifest.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.manifest.display()))?;
    if manifest.get("schema_version").and_then(Value::as_str) != Some("1.0") {
        anyhow::bail!("clinical manifest requires schema_version 1.0");
    }
    let base = args.manifest.parent().unwrap_or_else(|| Path::new("."));
    let input = manifest
        .get("input")
        .and_then(Value::as_str)
        .context("clinical manifest requires string input")?;
    let input = base.join(input);
    if !input.exists() {
        anyhow::bail!(
            "clinical manifest input does not exist: {}",
            input.display()
        );
    }
    let cases = manifest
        .get("cases")
        .and_then(Value::as_array)
        .context("clinical manifest requires cases array")?;
    if cases.is_empty() {
        anyhow::bail!("clinical manifest contains no cases");
    }
    for case in cases {
        let name = case
            .get("name")
            .and_then(Value::as_str)
            .context("clinical case requires string name")?;
        let annotations = case
            .get("annotations")
            .and_then(Value::as_str)
            .context("clinical case requires string annotations")?;
        let expected_code = case
            .get("expected_error_code")
            .and_then(Value::as_str)
            .context("clinical case requires string expected_error_code")?;
        let annotations = base.join(annotations);
        if !annotations.exists() {
            anyhow::bail!(
                "clinical case '{name}' annotations do not exist: {}",
                annotations.display()
            );
        }
        let result = {
            let _environment =
                ScopedEnvironmentVariable::set("MITO_CLINICAL_ANNOTATIONS", &annotations);
            let engine = MitoEngine::new().context("failed to create analysis engine")?;
            engine.analyze_with_options(
                &input,
                None,
                AnalyzeOptions {
                    filter_numt: true,
                    threads: 1,
                    ..AnalyzeOptions::default()
                },
            )
        };
        match result {
            Ok(_) => anyhow::bail!("clinical negative case '{name}' unexpectedly succeeded"),
            Err(error) if error.code == expected_code => {}
            Err(error) => anyhow::bail!(
                "clinical negative case '{name}' returned {}, expected {}: {}",
                error.code,
                expected_code,
                error.message
            ),
        }
    }
    println!(
        "clinical negative manifest validation passed: {} cases",
        cases.len()
    );
    Ok(())
}

fn project_evidence(data: &Value) -> Result<Value> {
    let variants = require_array(data, "/variants")?
        .iter()
        .enumerate()
        .map(|(index, variant)| {
            project_object_fields(
                variant,
                &[
                    "position",
                    "ref",
                    "alt",
                    "alt_depth",
                    "ref_depth",
                    "other_depth",
                    "callable_depth",
                    "heteroplasmy",
                    "ci95_low",
                    "ci95_high",
                    "molecule_support",
                    "strand_support",
                    "strand_bias_delta",
                    "allele_quality",
                    "read_position",
                    "supporting_reads",
                ],
                &format!("/variants/{index}"),
            )
        })
        .collect::<Result<Vec<_>>>()?;

    let reads = require_array(data, "/reads")?
        .iter()
        .enumerate()
        .map(|(index, read)| {
            let mut projected = project_object_fields(
                read,
                &["id", "flags", "mapping_quality"],
                &format!("/reads/{index}"),
            )?;
            let source_snps = read
                .get("snps")
                .and_then(Value::as_array)
                .with_context(|| format!("result contract violation: /reads/{index}/snps"))?;
            let snps = source_snps
                .iter()
                .enumerate()
                .map(|(snp_index, snp)| {
                    project_object_fields(
                        snp,
                        &["position", "ref", "alt"],
                        &format!("/reads/{index}/snps/{snp_index}"),
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            projected
                .as_object_mut()
                .context("internal evidence projection must be an object")?
                .insert("snps".to_owned(), Value::Array(snps));
            Ok(projected)
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(serde_json::json!({
        "schema_version": data.pointer("/metadata/schema_version"),
        "input_alignment_records": data.pointer("/filter_stats/input_alignment_records"),
        "input_molecules": data.pointer("/filter_stats/input_molecules"),
        "passed_reads": data.pointer("/filter_stats/passed_reads"),
        "numt_filtered_reads": data.pointer("/filter_stats/numt_filtered_reads"),
        "variants": variants,
        "reads": reads,
    }))
}

fn materialize_observations(data: &Value) -> Result<Vec<Value>> {
    if let Some(observations) = data.get("observations").and_then(Value::as_array) {
        return Ok(observations.clone());
    }
    let pages = require_array(data, "/observation_pages")?;
    let column_names = [
        "molecule_id",
        "event_id",
        "alignment_id",
        "state",
        "observed_allele",
        "base_quality",
        "mapping_quality",
        "strand",
        "evidence_source",
        "read_position",
    ];
    let expected_count = data
        .pointer("/evidence_encoding/observation_count")
        .and_then(Value::as_u64)
        .context("schema 0.6 observation_count must be an integer")?;
    let expected_page_count = data
        .pointer("/evidence_encoding/observation_page_count")
        .and_then(Value::as_u64)
        .context("schema 0.6 observation_page_count must be an integer")?;
    let page_size = data
        .pointer("/evidence_encoding/observation_page_size")
        .and_then(Value::as_u64)
        .filter(|value| (1..=1_000_000).contains(value))
        .context("schema 0.6 observation_page_size is invalid")?;
    if pages.len() as u64 != expected_page_count {
        anyhow::bail!("schema 0.6 observation page cardinality does not match metadata");
    }
    let capacity = usize::try_from(expected_count)
        .context("schema 0.6 observation_count exceeds platform size")?;
    let mut observations = Vec::with_capacity(capacity);
    for (page_index, page) in pages.iter().enumerate() {
        let index = page
            .get("index")
            .and_then(Value::as_u64)
            .with_context(|| format!("/observation_pages/{page_index}/index must be an integer"))?;
        let offset = page
            .get("offset")
            .and_then(Value::as_u64)
            .with_context(|| {
                format!("/observation_pages/{page_index}/offset must be an integer")
            })?;
        let count = page
            .get("count")
            .and_then(Value::as_u64)
            .with_context(|| format!("/observation_pages/{page_index}/count must be an integer"))?;
        if index != page_index as u64 || offset != observations.len() as u64 {
            anyhow::bail!("schema 0.6 observation pages are not contiguous");
        }
        let count =
            usize::try_from(count).context("observation page count exceeds platform size")?;
        if count as u64 > page_size {
            anyhow::bail!("schema 0.6 observation page {page_index} exceeds the page-size limit");
        }
        let columns = page
            .get("columns")
            .and_then(Value::as_object)
            .with_context(|| {
                format!("/observation_pages/{page_index}/columns must be an object")
            })?;
        for name in column_names {
            if columns.get(name).and_then(Value::as_array).map(Vec::len) != Some(count) {
                anyhow::bail!(
                    "schema 0.6 observation page {page_index} column {name} has the wrong length"
                );
            }
        }
        for row_index in 0..count {
            let mut row = serde_json::Map::new();
            row.insert(
                "id".to_owned(),
                Value::String(format!("observation:{}", observations.len())),
            );
            for name in column_names {
                row.insert(name.to_owned(), columns[name][row_index].clone());
            }
            observations.push(Value::Object(row));
        }
    }
    if observations.len() != capacity {
        anyhow::bail!("schema 0.6 observation page count does not match observation_count");
    }
    Ok(observations)
}

fn project_evidence_graph(data: &Value) -> Result<Value> {
    let project_array = |pointer: &str, fields: &[&str]| -> Result<Vec<Value>> {
        require_array(data, pointer)?
            .iter()
            .enumerate()
            .map(|(index, value)| {
                project_object_fields(value, fields, &format!("{pointer}/{index}"))
            })
            .collect()
    };

    let alignments = project_array(
        "/alignments",
        &[
            "id",
            "source_record_index",
            "molecule_id",
            "role",
            "selected_representative",
            "flags",
        ],
    )?;
    let molecules = project_array(
        "/molecules",
        &[
            "id",
            "index",
            "identity_policy",
            "assembly_status",
            "primary_candidate_count",
            "ambiguous",
            "analysis_eligible",
            "evidence_eligible",
            "callability_status",
            "callable_bases",
            "callable_fraction",
            "query_length",
            "mean_base_quality",
            "mapping_quality",
            "numt_score",
            "numt_evidence",
            "cluster_id",
            "alternate_event_ids",
            "evidence_state_counts",
            "representative_alignment_id",
            "source_qnames",
            "protocol_metadata",
            "protocol_flags",
            "exclusion_reasons",
            "alignment_ids",
            "warnings",
        ],
    )?;
    let events = project_array(
        "/events",
        &[
            "id",
            "index",
            "type",
            "start",
            "end",
            "length",
            "ref",
            "alt",
            "normalization",
            "source_projection",
            "negative_evidence_rule",
            "assessability",
            "component_event_ids",
            "supporting_molecule_ids",
            "evidence_counts",
        ],
    )?;
    let observations = materialize_observations(data)?
        .iter()
        .enumerate()
        .map(|(index, value)| {
            project_object_fields(
                value,
                &[
                    "id",
                    "molecule_id",
                    "event_id",
                    "alignment_id",
                    "state",
                    "observed_allele",
                    "base_quality",
                    "mapping_quality",
                    "strand",
                    "evidence_source",
                ],
                &format!("/observation_pages/materialized/{index}"),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let variants = project_array(
        "/variants",
        &[
            "event_id",
            "type",
            "position",
            "start",
            "end",
            "length",
            "ref",
            "alt",
            "normalization",
            "negative_evidence_rule",
            "assessability",
            "alt_depth",
            "ref_depth",
            "other_depth",
            "event_absent_depth",
            "low_quality_depth",
            "conflict_depth",
            "callable_depth",
            "heteroplasmy",
            "ci95_low",
            "ci95_high",
            "filter_status",
            "qc_flags",
            "numt_assessability",
            "multi_allelic",
            "homopolymer_context",
            "supporting_molecule_ids",
            "supporting_reads",
        ],
    )?;

    Ok(serde_json::json!({
        "schema_version": data.pointer("/metadata/schema_version"),
        "input_alignment_records": data.pointer("/filter_stats/input_alignment_records"),
        "input_molecules": data.pointer("/filter_stats/input_molecules"),
        "evidence_encoding": data.get("evidence_encoding"),
        "alignments": alignments,
        "molecules": molecules,
        "callability": data.get("callability"),
        "events": events,
        "observations": observations,
        "variants": variants,
        "phase_links": data.get("phase_links"),
        "architectures": data.get("architectures"),
    }))
}

fn project_object_fields(value: &Value, fields: &[&str], path: &str) -> Result<Value> {
    let source = value
        .as_object()
        .with_context(|| format!("result contract violation: {path} must be an object"))?;
    let mut projected = serde_json::Map::new();
    for field in fields {
        projected.insert(
            (*field).to_owned(),
            source
                .get(*field)
                .cloned()
                .with_context(|| format!("result contract violation: {path}/{field} is missing"))?,
        );
    }
    Ok(Value::Object(projected))
}

fn validate_error_manifest(args: ValidateErrorManifestArgs) -> Result<()> {
    if !args.manifest.exists() {
        anyhow::bail!("error manifest does not exist: {}", args.manifest.display());
    }
    let manifest: Value = serde_json::from_str(
        &fs::read_to_string(&args.manifest)
            .with_context(|| format!("failed to read {}", args.manifest.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.manifest.display()))?;
    let expected_schema = manifest
        .get("error_schema_version")
        .and_then(Value::as_str)
        .context("error manifest requires string error_schema_version")?;
    let capabilities = MitoEngine::capabilities();
    if expected_schema != capabilities.error_schema_version {
        anyhow::bail!(
            "error-schema mismatch: manifest {}, engine {}",
            expected_schema,
            capabilities.error_schema_version
        );
    }
    let cases = manifest
        .get("cases")
        .and_then(Value::as_array)
        .context("error manifest requires a cases array")?;
    if cases.is_empty() {
        anyhow::bail!("error manifest cases array must not be empty");
    }

    let base = args.manifest.parent().unwrap_or_else(|| Path::new("."));
    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    let mut case_names = BTreeSet::new();
    for (index, case) in cases.iter().enumerate() {
        let object = case
            .as_object()
            .with_context(|| format!("error manifest case {index} must be an object"))?;
        let name = object
            .get("name")
            .and_then(Value::as_str)
            .with_context(|| format!("error manifest case {index} requires string name"))?;
        if !case_names.insert(name.to_owned()) {
            anyhow::bail!("duplicate error manifest case name: {name}");
        }
        let input = object
            .get("input")
            .and_then(Value::as_str)
            .with_context(|| format!("error manifest case '{name}' requires string input"))?;
        let expected_code = object
            .get("code")
            .and_then(Value::as_str)
            .with_context(|| format!("error manifest case '{name}' requires string code"))?;
        let input_path = {
            let path = Path::new(input);
            if path.is_absolute() {
                path.to_path_buf()
            } else {
                base.join(path)
            }
        };
        if !input_path.is_file() {
            anyhow::bail!(
                "error manifest case '{name}' input does not exist: {}",
                input_path.display()
            );
        }
        let error = match engine.analyze_with_options(&input_path, None, AnalyzeOptions::default())
        {
            Ok(_) => {
                anyhow::bail!("negative golden '{name}' unexpectedly produced a scientific result")
            }
            Err(error) => error,
        };
        if error.code != expected_code {
            anyhow::bail!(
                "negative golden '{name}' expected {expected_code}, got {}: {}",
                error.code,
                error.message
            );
        }
        println!("negative fixture passed: {name} [{expected_code}]");
    }

    println!(
        "error manifest validation passed: {} cases from {}",
        cases.len(),
        args.manifest.display()
    );
    Ok(())
}

#[derive(Debug)]
struct HsdProfile {
    ranges: BTreeSet<String>,
    markers: BTreeSet<String>,
}

fn canonical_hsd_marker(sample: &str, raw: &str) -> Result<String> {
    let marker = raw.trim();
    if marker.is_empty() || marker.chars().any(char::is_whitespace) {
        anyhow::bail!("HSD profile '{sample}' has invalid marker {marker:?}");
    }
    let bytes = marker.as_bytes();
    if bytes.len() > 3 && bytes[bytes.len() - 3..].eq_ignore_ascii_case(b"DEL") {
        let position = marker
            .get(..bytes.len() - 3)
            .with_context(|| format!("HSD profile '{sample}' has invalid deletion {marker:?}"))?;
        let parsed = position
            .parse::<u64>()
            .with_context(|| format!("HSD profile '{sample}' has invalid deletion {marker:?}"))?;
        if parsed == 0 || parsed > 16_569 {
            anyhow::bail!("HSD profile '{sample}' has invalid deletion {marker:?}");
        }
        return Ok(format!("{parsed}d"));
    }
    Ok(marker.to_owned())
}

fn canonical_hsd_ranges(sample: &str, raw: &str) -> Result<BTreeSet<String>> {
    let mut ranges = BTreeSet::new();
    let unquoted = raw.trim().trim_matches('"');
    for raw_range in unquoted.split(';') {
        let compact: String = raw_range
            .chars()
            .filter(|value| !value.is_whitespace())
            .collect();
        if compact.is_empty() {
            continue;
        }
        let (start, end) = compact.split_once('-').unwrap_or((&compact, &compact));
        let start = start
            .parse::<u64>()
            .with_context(|| format!("HSD profile '{sample}' has invalid range {compact:?}"))?;
        let end = end
            .parse::<u64>()
            .with_context(|| format!("HSD profile '{sample}' has invalid range {compact:?}"))?;
        if start == 0 || end == 0 || start > 16_569 || end > 16_569 {
            anyhow::bail!("HSD profile '{sample}' has invalid range {compact:?}");
        }
        if start <= end {
            ranges.insert(format!("{start}-{end}"));
        } else {
            ranges.insert(format!("1-{end}"));
            ranges.insert(format!("{start}-16569"));
        }
    }
    if ranges.is_empty() {
        anyhow::bail!("HSD profile '{sample}' has no tested range");
    }
    Ok(ranges)
}

fn load_hsd_profiles(path: &Path) -> Result<BTreeMap<String, HsdProfile>> {
    let input = fs::File::open(path)
        .with_context(|| format!("failed to open pinned HSD input {}", path.display()))?;
    let mut profiles = BTreeMap::new();
    for (line_index, line) in BufReader::new(input).lines().enumerate() {
        let line = line.with_context(|| {
            format!(
                "failed to read pinned HSD input {} at line {}",
                path.display(),
                line_index + 1
            )
        })?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            anyhow::bail!(
                "pinned HSD input {} line {} requires sample, range, metadata, and markers",
                path.display(),
                line_index + 1
            );
        }
        let sample = fields[0].trim();
        if sample.is_empty() || sample.chars().any(char::is_whitespace) {
            anyhow::bail!(
                "pinned HSD input {} line {} has invalid sample identifier",
                path.display(),
                line_index + 1
            );
        }
        let mut markers = BTreeSet::new();
        for raw_marker in &fields[3..] {
            if raw_marker.trim().is_empty() {
                continue;
            }
            let marker = canonical_hsd_marker(sample, raw_marker)?;
            if !markers.insert(marker.clone()) {
                anyhow::bail!("HSD profile '{sample}' repeats marker {marker}");
            }
        }
        if markers.is_empty() {
            anyhow::bail!("HSD profile '{sample}' has no markers");
        }
        let profile = HsdProfile {
            ranges: canonical_hsd_ranges(sample, fields[1])?,
            markers,
        };
        if profiles.insert(sample.to_owned(), profile).is_some() {
            anyhow::bail!("duplicate HSD profile: {sample}");
        }
    }
    if profiles.is_empty() {
        anyhow::bail!("pinned HSD input is empty: {}", path.display());
    }
    Ok(profiles)
}

fn validate_haplogroup_manifest(args: ValidateHaplogroupManifestArgs) -> Result<()> {
    if !args.manifest.is_file() {
        anyhow::bail!(
            "haplogroup manifest does not exist: {}",
            args.manifest.display()
        );
    }
    let manifest: Value = serde_json::from_str(
        &fs::read_to_string(&args.manifest)
            .with_context(|| format!("failed to read {}", args.manifest.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.manifest.display()))?;
    if manifest.get("schema_version").and_then(Value::as_str) != Some("1.0") {
        anyhow::bail!("haplogroup manifest requires schema_version 1.0");
    }
    let reference = manifest
        .get("reference")
        .and_then(Value::as_object)
        .context("haplogroup manifest requires a reference object")?;
    for field in [
        "tool",
        "version",
        "tree",
        "asset_sha256",
        "tree_sha256",
        "weights_sha256",
        "alignment_rules_sha256",
        "source_hsd",
        "source_hsd_sha256",
    ] {
        if reference.get(field).and_then(Value::as_str).is_none() {
            anyhow::bail!("haplogroup manifest reference requires string {field}");
        }
    }
    if reference.get("tool").and_then(Value::as_str) != Some("HaploGrep 3")
        || reference.get("version").and_then(Value::as_str) != Some("3.3.2")
        || reference.get("tree").and_then(Value::as_str) != Some("phylotree-rcrs@17.3")
    {
        anyhow::bail!("haplogroup manifest must pin HaploGrep 3.3.2 with phylotree-rcrs@17.3");
    }
    for field in [
        "asset_sha256",
        "tree_sha256",
        "weights_sha256",
        "alignment_rules_sha256",
        "source_hsd_sha256",
    ] {
        let digest = reference
            .get(field)
            .and_then(Value::as_str)
            .context("validated reference digest must be present")?;
        if digest.len() != 64 || !digest.chars().all(|value| value.is_ascii_hexdigit()) {
            anyhow::bail!("haplogroup manifest reference {field} is not a SHA-256 digest");
        }
    }
    let base = args.manifest.parent().unwrap_or_else(|| Path::new("."));
    let source_hsd = base.join(
        reference
            .get("source_hsd")
            .and_then(Value::as_str)
            .context("source_hsd must be present")?,
    );
    if !source_hsd.is_file() {
        anyhow::bail!(
            "pinned HaploGrep input is missing: {}",
            source_hsd.display()
        );
    }
    let source_digest = sha256_file(&source_hsd)?;
    let expected_source_digest = reference
        .get("source_hsd_sha256")
        .and_then(Value::as_str)
        .context("source_hsd_sha256 must be present")?;
    if source_digest != expected_source_digest {
        anyhow::bail!(
            "pinned HaploGrep input checksum mismatch: expected {expected_source_digest}, got {source_digest}"
        );
    }
    let hsd_profiles = load_hsd_profiles(&source_hsd)?;

    let cases = manifest
        .get("cases")
        .and_then(Value::as_array)
        .context("haplogroup manifest requires a cases array")?;
    if cases.is_empty() {
        anyhow::bail!("haplogroup manifest cases array must not be empty");
    }
    let coverage = manifest
        .get("coverage")
        .and_then(Value::as_object)
        .context("haplogroup manifest requires a coverage object")?;
    let expected_differential_profiles = coverage
        .get("differential_profiles")
        .and_then(Value::as_u64)
        .context("haplogroup manifest coverage requires differential_profiles")?;
    let expected_alignment_cases = coverage
        .get("alignment_projection_cases")
        .and_then(Value::as_u64)
        .context("haplogroup manifest coverage requires alignment_projection_cases")?;
    let expected_compound_rule_cases = coverage
        .get("compound_rule_cases")
        .and_then(Value::as_u64)
        .context("haplogroup manifest coverage requires compound_rule_cases")?;
    let expected_backmutation_cases = coverage
        .get("backmutation_cases")
        .and_then(Value::as_u64)
        .context("haplogroup manifest coverage requires backmutation_cases")?;
    let expected_lineage_groups = coverage
        .get("lineage_groups")
        .and_then(Value::as_array)
        .context("haplogroup manifest coverage requires lineage_groups")?
        .iter()
        .map(|value| {
            value
                .as_str()
                .map(str::to_owned)
                .context("haplogroup manifest coverage has a non-string lineage group")
        })
        .collect::<Result<BTreeSet<_>>>()?;
    if expected_lineage_groups.is_empty() {
        anyhow::bail!("haplogroup manifest coverage has no lineage groups");
    }

    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    let mut case_names = BTreeSet::new();
    let mut referenced_hsd_profiles = BTreeSet::new();
    let mut differential_cases = 0usize;
    let mut alignment_cases = 0usize;
    let mut compound_rule_cases = 0usize;
    let mut backmutation_cases = 0usize;
    let mut observed_lineage_groups = BTreeSet::new();
    for (index, case) in cases.iter().enumerate() {
        let object = case
            .as_object()
            .with_context(|| format!("haplogroup case {index} must be an object"))?;
        let name = object
            .get("name")
            .and_then(Value::as_str)
            .with_context(|| format!("haplogroup case {index} requires string name"))?;
        if !case_names.insert(name.to_owned()) {
            anyhow::bail!("duplicate haplogroup case name: {name}");
        }
        let source_profile = object
            .get("source_profile")
            .map(|value| {
                value.as_str().with_context(|| {
                    format!("haplogroup case '{name}' has a non-string source_profile")
                })
            })
            .transpose()?;
        let (expected_markers, expected_ranges) = if let Some(source_profile) = source_profile {
            if object.contains_key("markers") || object.contains_key("range") {
                anyhow::bail!(
                    "haplogroup case '{name}' must not duplicate markers/range from source_profile"
                );
            }
            if !referenced_hsd_profiles.insert(source_profile.to_owned()) {
                anyhow::bail!("HSD profile '{source_profile}' is referenced more than once");
            }
            let profile = hsd_profiles.get(source_profile).with_context(|| {
                format!(
                    "haplogroup case '{name}' references missing HSD profile '{source_profile}'"
                )
            })?;
            (profile.markers.clone(), profile.ranges.clone())
        } else {
            let marker_values = object
                .get("markers")
                .and_then(Value::as_array)
                .with_context(|| format!("haplogroup case '{name}' requires markers"))?;
            if marker_values.is_empty() {
                anyhow::bail!("haplogroup case '{name}' has no markers");
            }
            let mut markers = BTreeSet::new();
            for marker in marker_values {
                let marker = marker
                    .as_str()
                    .with_context(|| format!("haplogroup case '{name}' has a non-string marker"))?;
                if marker.is_empty() || marker.chars().any(char::is_whitespace) {
                    anyhow::bail!("haplogroup case '{name}' has an invalid marker: {marker:?}");
                }
                if !markers.insert(marker.to_owned()) {
                    anyhow::bail!("haplogroup case '{name}' repeats marker {marker}");
                }
            }
            let range_values = object
                .get("range")
                .and_then(Value::as_array)
                .with_context(|| format!("haplogroup case '{name}' requires range"))?;
            if range_values.is_empty() {
                anyhow::bail!("haplogroup case '{name}' has no tested range");
            }
            let mut ranges = BTreeSet::new();
            for range in range_values {
                let range = range
                    .as_str()
                    .with_context(|| format!("haplogroup case '{name}' has a non-string range"))?;
                let (start, end) = range.split_once('-').unwrap_or((range, range));
                let start = start.parse::<u64>().with_context(|| {
                    format!("haplogroup case '{name}' has invalid range {range:?}")
                })?;
                let end = end.parse::<u64>().with_context(|| {
                    format!("haplogroup case '{name}' has invalid range {range:?}")
                })?;
                if start == 0 || end < start || end > 16_569 {
                    anyhow::bail!("haplogroup case '{name}' has invalid range {range:?}");
                }
                ranges.insert(format!("{start}-{end}"));
            }
            (markers, ranges)
        };
        let expected_candidates = object
            .get("haplogrep_candidates")
            .map(|value| {
                let values = value.as_array().with_context(|| {
                    format!("haplogroup case '{name}' has invalid haplogrep_candidates")
                })?;
                values
                    .iter()
                    .map(|candidate| {
                        candidate.as_str().map(str::to_owned).with_context(|| {
                            format!("haplogroup case '{name}' has a non-string reference candidate")
                        })
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .transpose()?;
        let backmutation_absent_markers = object
            .get("backmutation_absent_markers")
            .map(|value| {
                value
                    .as_array()
                    .with_context(|| {
                        format!("haplogroup case '{name}' has invalid backmutation_absent_markers")
                    })?
                    .iter()
                    .map(|marker| {
                        marker.as_str().map(str::to_owned).with_context(|| {
                            format!("haplogroup case '{name}' has a non-string backmutation marker")
                        })
                    })
                    .collect::<Result<BTreeSet<_>>>()
            })
            .transpose()?
            .unwrap_or_default();
        if !backmutation_absent_markers.is_empty() {
            if source_profile.is_none() {
                anyhow::bail!("haplogroup backmutation case '{name}' requires a source_profile");
            }
            for marker in &backmutation_absent_markers {
                if expected_markers.contains(marker) {
                    anyhow::bail!(
                        "haplogroup backmutation case '{name}' still contains removed marker {marker}"
                    );
                }
            }
            backmutation_cases += 1;
        }
        if expected_candidates.as_ref().is_some_and(Vec::is_empty) {
            anyhow::bail!("haplogroup case '{name}' has no reference top hit");
        }
        if let Some(candidates) = &expected_candidates {
            if candidates.len() != 3 {
                anyhow::bail!(
                    "haplogroup case '{name}' must pin exactly three HaploGrep candidates"
                );
            }
            let quality = object
                .get("haplogrep_quality")
                .and_then(Value::as_f64)
                .with_context(|| {
                    format!("haplogroup case '{name}' requires numeric haplogrep_quality")
                })?;
            if !quality.is_finite() || quality <= 0.0 {
                anyhow::bail!("haplogroup case '{name}' has invalid haplogrep_quality");
            }
            let lineage_group = object
                .get("lineage_group")
                .and_then(Value::as_str)
                .with_context(|| {
                    format!("haplogroup case '{name}' requires string lineage_group")
                })?;
            if lineage_group.is_empty() || lineage_group.chars().any(char::is_whitespace) {
                anyhow::bail!("haplogroup case '{name}' has invalid lineage_group");
            }
            observed_lineage_groups.insert(lineage_group.to_owned());
        }

        let input = object.get("input").and_then(Value::as_str);
        if source_profile.is_some() && input.is_some() {
            anyhow::bail!("haplogroup case '{name}' cannot combine source_profile with input");
        }
        if source_profile.is_some() && expected_candidates.is_none() {
            anyhow::bail!("haplogroup case '{name}' source_profile requires reference candidates");
        }
        if input.is_none() && expected_candidates.is_none() {
            anyhow::bail!("haplogroup case '{name}' requires either input or haplogrep_candidates");
        }
        let rule_class = object
            .get("rule_class")
            .map(|value| {
                value
                    .as_str()
                    .with_context(|| format!("haplogroup case '{name}' has invalid rule_class"))
            })
            .transpose()?;
        if input.is_some() {
            if !matches!(rule_class, Some("cigar_3prime" | "compound")) {
                anyhow::bail!(
                    "haplogroup alignment case '{name}' requires rule_class cigar_3prime or compound"
                );
            }
            let normalization_rules = object
                .get("normalization_rules")
                .and_then(Value::as_array)
                .with_context(|| {
                    format!("haplogroup alignment case '{name}' requires normalization_rules")
                })?;
            if normalization_rules.is_empty()
                || normalization_rules.iter().any(|value| {
                    value
                        .as_str()
                        .map_or(true, |rule| rule.trim().is_empty() || !rule.contains("->"))
                })
            {
                anyhow::bail!("haplogroup alignment case '{name}' has invalid normalization_rules");
            }
        } else if rule_class.is_some() {
            anyhow::bail!("haplogroup differential case '{name}' must not set rule_class");
        }
        let temporary_fixture = input.is_none();
        let fixture_path = if let Some(input) = input {
            base.join(input)
        } else {
            std::env::temp_dir().join(format!(
                "mito-haplogroup-{}-{index}.fastq",
                std::process::id()
            ))
        };
        if temporary_fixture {
            let mut fastq = format!(
                "@{name} phylo_range={}",
                expected_ranges
                    .iter()
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(",")
            );
            for marker in &expected_markers {
                fastq.push_str(" phylo=");
                fastq.push_str(marker);
            }
            fastq.push_str("\nA\n+\nI\n");
            fs::write(&fixture_path, fastq)
                .with_context(|| format!("failed to write {}", fixture_path.display()))?;
        } else if !fixture_path.is_file() {
            anyhow::bail!(
                "haplogroup case '{name}' input does not exist: {}",
                fixture_path.display()
            );
        }
        let analysis = engine.analyze_with_options(
            &fixture_path,
            None,
            AnalyzeOptions {
                filter_numt: false,
                threads: 1,
                allow_development_tags: true,
                ..AnalyzeOptions::default()
            },
        );
        if temporary_fixture {
            let _ = fs::remove_file(&fixture_path);
        }
        let json = analysis.with_context(|| format!("haplogroup case '{name}' failed"))?;
        let data: Value = serde_json::from_str(&json)
            .with_context(|| format!("haplogroup case '{name}' returned invalid JSON"))?;
        validate_result_contract(&data)?;
        let clusters = require_array(&data, "/clusters")?;
        if clusters.len() != 1 {
            anyhow::bail!(
                "haplogroup case '{name}' expected one cluster, got {}",
                clusters.len()
            );
        }
        let assignment = clusters[0]
            .get("haplogroup_assignment")
            .and_then(Value::as_object)
            .with_context(|| format!("haplogroup case '{name}' has no assignment object"))?;
        let actual_best = clusters[0]
            .get("haplogroup")
            .and_then(Value::as_str)
            .with_context(|| format!("haplogroup case '{name}' has no best assignment"))?;
        if let Some(expected_best) = expected_candidates
            .as_ref()
            .and_then(|candidates| candidates.first())
        {
            if actual_best != expected_best {
                anyhow::bail!(
                    "HaploGrep differential mismatch for '{name}': expected {expected_best}, got {actual_best}"
                );
            }
        }
        if let Some(expected_candidates) = &expected_candidates {
            let candidate_values = assignment
                .get("candidates")
                .and_then(Value::as_array)
                .with_context(|| format!("haplogroup case '{name}' has no candidates array"))?;
            if !backmutation_absent_markers.is_empty() {
                let winning_candidate = candidate_values.first().with_context(|| {
                    format!("haplogroup backmutation case '{name}' has no winning candidate")
                })?;
                for field in ["matched", "missing"] {
                    let values = winning_candidate
                        .get(field)
                        .and_then(Value::as_array)
                        .with_context(|| {
                            format!(
                                "haplogroup backmutation case '{name}' winner has no {field} array"
                            )
                        })?;
                    for marker in values.iter().filter_map(Value::as_str) {
                        if backmutation_absent_markers.contains(marker) {
                            anyhow::bail!(
                                "haplogroup backmutation case '{name}' emitted removed marker {marker} in {field}"
                            );
                        }
                    }
                }
            }
            let actual_candidates = candidate_values
                .iter()
                .map(|candidate| {
                    candidate
                        .get("name")
                        .and_then(Value::as_str)
                        .map(str::to_owned)
                        .with_context(|| {
                            format!("haplogroup case '{name}' emitted a candidate without name")
                        })
                })
                .collect::<Result<Vec<_>>>()?;
            if &actual_candidates != expected_candidates {
                anyhow::bail!(
                    "HaploGrep top-3 mismatch for '{name}': expected {:?}, got {:?}",
                    expected_candidates,
                    actual_candidates
                );
            }
        }
        let actual_markers = assignment
            .get("observed_markers")
            .and_then(Value::as_array)
            .with_context(|| format!("haplogroup case '{name}' has no observed_markers"))?
            .iter()
            .map(|marker| {
                marker.as_str().map(str::to_owned).with_context(|| {
                    format!("haplogroup case '{name}' emitted a non-string marker")
                })
            })
            .collect::<Result<BTreeSet<_>>>()?;
        if actual_markers != expected_markers {
            anyhow::bail!(
                "haplogroup marker projection mismatch for '{name}': expected {:?}, got {:?}",
                expected_markers,
                actual_markers
            );
        }
        let actual_ranges = assignment
            .get("callable_ranges")
            .and_then(Value::as_array)
            .with_context(|| format!("haplogroup case '{name}' has no callable_ranges"))?
            .iter()
            .map(|range| {
                let start = range
                    .get("start")
                    .and_then(Value::as_u64)
                    .with_context(|| {
                        format!("haplogroup case '{name}' emitted a range without start")
                    })?;
                let end = range.get("end").and_then(Value::as_u64).with_context(|| {
                    format!("haplogroup case '{name}' emitted a range without end")
                })?;
                Ok(format!("{start}-{end}"))
            })
            .collect::<Result<BTreeSet<_>>>()?;
        if actual_ranges != expected_ranges {
            anyhow::bail!(
                "haplogroup callable-range mismatch for '{name}': expected {:?}, got {:?}",
                expected_ranges,
                actual_ranges
            );
        }
        if expected_candidates.is_some() {
            differential_cases += 1;
            println!("haplogroup differential passed: {name} -> {actual_best}");
        } else {
            alignment_cases += 1;
            if rule_class == Some("compound") {
                compound_rule_cases += 1;
            }
            println!("haplogroup alignment projection passed: {name}");
        }
    }

    let available_hsd_profiles = hsd_profiles.keys().cloned().collect::<BTreeSet<_>>();
    if referenced_hsd_profiles != available_hsd_profiles {
        let missing = available_hsd_profiles
            .difference(&referenced_hsd_profiles)
            .cloned()
            .collect::<Vec<_>>();
        let unknown = referenced_hsd_profiles
            .difference(&available_hsd_profiles)
            .cloned()
            .collect::<Vec<_>>();
        anyhow::bail!(
            "haplogroup manifest/HSD profile mismatch: unreferenced {:?}, unknown {:?}",
            missing,
            unknown
        );
    }
    if u64::try_from(differential_cases).context("differential case count overflow")?
        != expected_differential_profiles
        || u64::try_from(alignment_cases).context("alignment case count overflow")?
            != expected_alignment_cases
        || u64::try_from(compound_rule_cases).context("compound rule case count overflow")?
            != expected_compound_rule_cases
        || u64::try_from(backmutation_cases).context("backmutation case count overflow")?
            != expected_backmutation_cases
    {
        anyhow::bail!(
            "haplogroup coverage count mismatch: expected {expected_differential_profiles} differential/{expected_alignment_cases} alignment/{expected_compound_rule_cases} compound/{expected_backmutation_cases} backmutation, got {differential_cases}/{alignment_cases}/{compound_rule_cases}/{backmutation_cases}"
        );
    }
    if observed_lineage_groups != expected_lineage_groups {
        anyhow::bail!(
            "haplogroup lineage coverage mismatch: expected {:?}, got {:?}",
            expected_lineage_groups,
            observed_lineage_groups
        );
    }

    println!(
        "haplogroup manifest validation passed: {differential_cases} HaploGrep differential cases, {alignment_cases} alignment-projection cases ({compound_rule_cases} compound), and {backmutation_cases} explicit backmutation cases"
    );
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
    let retrieved_at_unix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system time is before the Unix epoch")?
        .as_secs();
    let write_result = write_clinvar_mtdna_tsv(&source, &temporary_output, retrieved_at_unix)
        .and_then(|()| validate_clinical_tsv_schema(&temporary_output))
        .and_then(|()| {
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
        "schema_version": 2,
        "clinical_annotation_schema_version": "1.0",
        "resource": "ClinVar variant_summary mitochondrial SNV assertions",
        "source": if downloaded {
            "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz".to_string()
        } else {
            source.display().to_string()
        },
        "retrieved_at_unix": retrieved_at_unix,
        "source_sha256": sha256_file(source)?,
        "cache_sha256": sha256_file(output)?,
        "normalization": "GRCh38 MT/M single-nucleotide records; PositionVCF and VCF alleles; one source record per assertion row; bundled curated assertions appended without field-wise conflict loss",
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

fn write_clinvar_mtdna_tsv(source_gz: &Path, output: &Path, retrieved_at_unix: u64) -> Result<()> {
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
    writeln!(out, "{CLINICAL_TSV_HEADER}")?;

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
                    let variation_id = clean_tsv(get("VariationID"));
                    let source_url = if !variation_id.is_empty()
                        && variation_id.chars().all(|value| value.is_ascii_digit())
                    {
                        format!("https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/")
                    } else {
                        String::new()
                    };
                    let row = [
                        clean_tsv(position),
                        clean_tsv(reference),
                        clean_tsv(alternate),
                        clean_tsv(get("GeneSymbol")),
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                        "ClinVar".to_string(),
                        clean_tsv(&get("RCVaccession").replace('|', ";")),
                        clean_tsv(get("AlleleID")),
                        clean_tsv(get("PhenotypeList")),
                        clean_tsv(get("ClinicalSignificance")),
                        clean_tsv(get("ReviewStatus")),
                        clean_tsv(get("LastEvaluated")),
                        source_url,
                        clean_tsv(&get("RCVaccession").replace('|', ";")),
                        "clinvar-variant-summary-grch38".to_string(),
                        format!("unix:{retrieved_at_unix}"),
                    ];
                    writeln!(out, "{}", row.join("\t"))?;
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

fn clean_export_field(value: &str) -> String {
    value
        .chars()
        .map(|character| match character {
            '\t' | '\n' | '\r' => ' ',
            other => other,
        })
        .collect::<String>()
        .trim()
        .to_owned()
}

fn render_variant_tsv(json: &str) -> Result<String> {
    let data: Value = serde_json::from_str(json).context("analysis returned invalid JSON")?;
    validate_result_contract(&data)?;
    let schema = require_string(&data, "/metadata/schema_version")?;
    let variants = require_array(&data, "/variants")?;
    let mut out = String::from(
        "event_id\ttype\tposition\tstart\tend\tref\talt\talt_molecules\tref_molecules\tother_molecules\tcallable_molecules\theteroplasmy\tci95_low\tci95_high\talt_forward\talt_reverse\tref_forward\tref_reverse\talt_mean_mapq\talt_mean_baseq\tmulti_allelic\thomopolymer_run\tnumt_assessability\tfilter_status\tqc_flags\tsupporting_molecule_ids\tgene\tconsequence\tclinical_significance\n",
    );
    for (index, variant) in variants.iter().enumerate() {
        let string = |field: &str| {
            variant
                .get(field)
                .and_then(Value::as_str)
                .map(clean_export_field)
                .unwrap_or_default()
        };
        let integer = |field: &str| {
            variant
                .get(field)
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default()
        };
        let number = |field: &str| {
            variant
                .get(field)
                .and_then(Value::as_f64)
                .map(|value| format!("{value:.8}"))
                .unwrap_or_default()
        };
        let event_id = variant
            .get("event_id")
            .and_then(Value::as_str)
            .map(clean_export_field)
            .unwrap_or_else(|| {
                format!(
                    "legacy:{}:{}:{}",
                    integer("position"),
                    string("ref"),
                    string("alt")
                )
            });
        let supporting = variant
            .get("supporting_molecule_ids")
            .or_else(|| variant.get("supporting_reads"))
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(Value::as_str)
            .map(clean_export_field)
            .collect::<Vec<_>>()
            .join(";");
        let qc_flags = variant
            .get("qc_flags")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(Value::as_str)
            .map(clean_export_field)
            .collect::<Vec<_>>()
            .join(";");
        let clinical = variant
            .pointer("/annotation/consensus_significance")
            .or_else(|| variant.pointer("/annotation/pathogenicity"))
            .and_then(Value::as_str)
            .map(clean_export_field)
            .unwrap_or_default();
        let fields = [
            event_id,
            if schema == "0.6" {
                string("type")
            } else {
                "SNV".to_owned()
            },
            integer("position"),
            integer("start"),
            integer("end"),
            string("ref"),
            string("alt"),
            integer("alt_depth"),
            integer("ref_depth"),
            integer("other_depth"),
            integer("callable_depth"),
            number("heteroplasmy"),
            number("ci95_low"),
            number("ci95_high"),
            variant
                .pointer("/strand_support/alt_forward")
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            variant
                .pointer("/strand_support/alt_reverse")
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            variant
                .pointer("/strand_support/ref_forward")
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            variant
                .pointer("/strand_support/ref_reverse")
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            variant
                .pointer("/mapping_quality/alternate/mean")
                .and_then(Value::as_f64)
                .map(|value| format!("{value:.4}"))
                .unwrap_or_default(),
            variant
                .pointer("/allele_quality/alternate/mean_phred")
                .and_then(Value::as_f64)
                .map(|value| format!("{value:.4}"))
                .unwrap_or_default(),
            variant
                .get("multi_allelic")
                .and_then(Value::as_bool)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            variant
                .pointer("/homopolymer_context/run_length")
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            string("numt_assessability"),
            string("filter_status"),
            qc_flags,
            supporting,
            string("gene"),
            string("consequence"),
            clinical,
        ];
        if fields
            .iter()
            .any(|field| field.contains('\t') || field.contains('\n') || field.contains('\r'))
        {
            anyhow::bail!("variant {index} contains an unsafe TSV field after sanitization");
        }
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }
    Ok(out)
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
    clinical_conflict: bool,
    clinical_sources: Vec<String>,
    ci95_low: Option<f64>,
    ci95_high: Option<f64>,
    event_id: Option<String>,
    event_type: Option<String>,
    reference_support: u64,
    other_support: u64,
    low_quality: u64,
    conflict: u64,
    qc_flags: Vec<String>,
    numt_assessability: Option<String>,
    normalization: Option<String>,
    alt_mapping_quality: Option<f64>,
    alt_base_quality: Option<f64>,
    strand_support: Option<[u64; 4]>,
    schema_0_6: bool,
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
    let schema_0_6 = data
        .pointer("/metadata/schema_version")
        .and_then(Value::as_str)
        == Some("0.6");
    let mut variants = BTreeMap::<(u64, String, String), VcfVariant>::new();
    let aggregates = data
        .get("variants")
        .and_then(Value::as_array)
        .context("result contract violation: /variants must be an array")?;
    for variant in aggregates {
        let Some(position) = variant
            .get(if schema_0_6 {
                "vcf_position"
            } else {
                "position"
            })
            .and_then(Value::as_u64)
        else {
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
        if schema_0_6 && variant.get("vcf_representable").and_then(Value::as_bool) == Some(false) {
            let event_id = variant
                .get("event_id")
                .and_then(Value::as_str)
                .unwrap_or("unknown");
            anyhow::bail!(
                "variant {event_id} crosses the circular origin and has no lossless linear VCF representation; use JSON or TSV"
            );
        }
        if reference.is_empty()
            || alternate.is_empty()
            || (!schema_0_6 && (reference.len() != 1 || alternate.len() != 1))
            || depth == 0
        {
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
                    .pointer("/annotation/consensus_significance")
                    .and_then(Value::as_str)
                    .or_else(|| {
                        variant
                            .pointer("/annotation/pathogenicity")
                            .and_then(Value::as_str)
                    })
                    .map(str::to_string),
                clinical_conflict: variant
                    .pointer("/annotation/conflict_status")
                    .and_then(Value::as_str)
                    == Some("conflicting"),
                clinical_sources: variant
                    .pointer("/annotation/sources")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                    .filter_map(Value::as_str)
                    .map(str::to_string)
                    .collect(),
                ci95_low: variant.get("ci95_low").and_then(Value::as_f64),
                ci95_high: variant.get("ci95_high").and_then(Value::as_f64),
                event_id: variant
                    .get("event_id")
                    .and_then(Value::as_str)
                    .map(str::to_string),
                event_type: variant
                    .get("type")
                    .and_then(Value::as_str)
                    .map(str::to_string)
                    .or_else(|| {
                        variant
                            .get("event_id")
                            .and_then(Value::as_str)
                            .map(|_| "SNV".to_string())
                    }),
                reference_support: variant
                    .get("ref_depth")
                    .and_then(Value::as_u64)
                    .unwrap_or(0),
                other_support: variant
                    .get("other_depth")
                    .and_then(Value::as_u64)
                    .unwrap_or(0),
                low_quality: variant
                    .get("low_quality_depth")
                    .and_then(Value::as_u64)
                    .unwrap_or(0),
                conflict: variant
                    .get("conflict_depth")
                    .and_then(Value::as_u64)
                    .unwrap_or(0),
                qc_flags: variant
                    .get("qc_flags")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                    .filter_map(Value::as_str)
                    .map(str::to_owned)
                    .collect(),
                numt_assessability: variant
                    .get("numt_assessability")
                    .and_then(Value::as_str)
                    .map(str::to_owned),
                normalization: variant
                    .get("normalization")
                    .and_then(Value::as_str)
                    .map(str::to_owned),
                alt_mapping_quality: variant
                    .pointer("/mapping_quality/alternate/mean")
                    .and_then(Value::as_f64),
                alt_base_quality: variant
                    .pointer("/allele_quality/alternate/mean_phred")
                    .and_then(Value::as_f64),
                strand_support: [
                    "/strand_support/alt_forward",
                    "/strand_support/alt_reverse",
                    "/strand_support/ref_forward",
                    "/strand_support/ref_reverse",
                ]
                .map(|pointer| variant.pointer(pointer).and_then(Value::as_u64))
                .into_iter()
                .collect::<Option<Vec<_>>>()
                .and_then(|values| values.try_into().ok()),
                schema_0_6,
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
    if schema_0_6 {
        out.push_str(
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate molecule support\">\n",
        );
        out.push_str("##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Reference and alternate molecule support\">\n");
        out.push_str("##INFO=<ID=ODC,Number=1,Type=Integer,Description=\"Callable molecules supporting another allele at the normalized event\">\n");
        out.push_str(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Event-callable molecule depth\">\n",
        );
        out.push_str("##INFO=<ID=LOWQ,Number=1,Type=Integer,Description=\"Low-quality molecule observations excluded from DP\">\n");
        out.push_str("##INFO=<ID=CONFLICT,Number=1,Type=Integer,Description=\"Conflicting molecule observations excluded from DP\">\n");
        out.push_str("##INFO=<ID=MOLECULE_SUPPORT,Number=1,Type=Integer,Description=\"Physical molecules supporting the alternate event\">\n");
        out.push_str("##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Mean mapping quality of alternate-supporting molecules\">\n");
        out.push_str("##INFO=<ID=BQ,Number=1,Type=Float,Description=\"Mean base quality of alternate-supporting observations when defined\">\n");
        out.push_str("##INFO=<ID=STRAND_SUPPORT,Number=4,Type=Integer,Description=\"ALT forward, ALT reverse, REF forward, REF reverse molecule counts\">\n");
        out.push_str("##INFO=<ID=NUMT_ASSESSABLE,Number=1,Type=String,Description=\"Whether competitive nuclear-plus-mitochondrial evidence makes NUMT specificity assessable\">\n");
        out.push_str("##INFO=<ID=QC_FLAGS,Number=.,Type=String,Description=\"Observed reason-coded QC facts; thresholds are not calibrated filters\">\n");
        out.push_str("##INFO=<ID=NORMALIZATION,Number=1,Type=String,Description=\"Normalized mitochondrial event representation rule\">\n");
    } else {
        out.push_str(
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate read support\">\n",
        );
        out.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Locus-callable A/C/G/T molecule depth after quality filters\">\n");
    }
    out.push_str("##INFO=<ID=HF,Number=A,Type=Float,Description=\"Heteroplasmy fraction estimated from alternate depth / locus-callable depth\">\n");
    out.push_str("##INFO=<ID=HF_CI95,Number=2,Type=Float,Description=\"Wilson 95% confidence interval for heteroplasmy fraction\">\n");
    out.push_str(
        "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Annotated mitochondrial gene\">\n",
    );
    out.push_str("##INFO=<ID=CLNSIG,Number=1,Type=String,Description=\"Deterministic clinical significance summary; inspect source assertions before interpretation\">\n");
    out.push_str("##INFO=<ID=CLNCONFLICT,Number=0,Type=Flag,Description=\"Source assertions contain incompatible or explicitly conflicting significance groups\">\n");
    out.push_str("##INFO=<ID=CLNSRC,Number=.,Type=String,Description=\"Clinical assertion sources preserved in the JSON result\">\n");
    const EVENT_TYPE_HEADER: &str = "##INFO=<ID=EVENT_TYPE,Number=1,Type=String,Description=\"Schema 0.6 normalized event class\">\n";
    if variants
        .values()
        .any(|variant| variant.event_type.is_some())
    {
        out.push_str(EVENT_TYPE_HEADER);
    }
    out.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype placeholder for haploid mtDNA export\">\n");
    if schema_0_6 {
        out.push_str("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Event-callable molecule depth\">\n");
        out.push_str("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Reference and alternate molecule support\">\n");
    }
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
        if variant.schema_0_6 {
            info.push(format!(
                "AD={},{}",
                variant.reference_support, variant.support
            ));
            info.push(format!("ODC={}", variant.other_support));
            info.push(format!("LOWQ={}", variant.low_quality));
            info.push(format!("CONFLICT={}", variant.conflict));
            info.push(format!("MOLECULE_SUPPORT={}", variant.support));
            if let Some(value) = variant.alt_mapping_quality {
                info.push(format!("MQ={value:.4}"));
            }
            if let Some(value) = variant.alt_base_quality {
                info.push(format!("BQ={value:.4}"));
            }
            if let Some(strands) = variant.strand_support {
                info.push(format!(
                    "STRAND_SUPPORT={},{},{},{}",
                    strands[0], strands[1], strands[2], strands[3]
                ));
            }
            if let Some(value) = &variant.numt_assessability {
                info.push(format!("NUMT_ASSESSABLE={}", sanitize_vcf_token(value)));
            }
            if !variant.qc_flags.is_empty() {
                info.push(format!(
                    "QC_FLAGS={}",
                    variant
                        .qc_flags
                        .iter()
                        .map(|value| sanitize_vcf_token(value))
                        .collect::<Vec<_>>()
                        .join(",")
                ));
            }
            if let Some(value) = &variant.normalization {
                info.push(format!("NORMALIZATION={}", sanitize_vcf_token(value)));
            }
        }
        if let Some(gene) = &variant.gene {
            info.push(format!("GENE={}", sanitize_vcf_token(gene)));
        }
        if let Some(clinical_significance) = &variant.clinical_significance {
            info.push(format!(
                "CLNSIG={}",
                sanitize_vcf_token(clinical_significance)
            ));
        }
        if variant.clinical_conflict {
            info.push("CLNCONFLICT".to_string());
        }
        if !variant.clinical_sources.is_empty() {
            info.push(format!(
                "CLNSRC={}",
                variant
                    .clinical_sources
                    .iter()
                    .map(|value| sanitize_vcf_token(value))
                    .collect::<Vec<_>>()
                    .join(",")
            ));
        }
        if let (Some(low), Some(high)) = (variant.ci95_low, variant.ci95_high) {
            info.push(format!("HF_CI95={low:.6},{high:.6}"));
        }
        if let Some(event_type) = &variant.event_type {
            info.push(format!("EVENT_TYPE={}", sanitize_vcf_token(event_type)));
        }
        let event_id = variant
            .event_id
            .as_deref()
            .map(sanitize_vcf_token)
            .unwrap_or_else(|| ".".to_string());
        if variant.schema_0_6 {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t.\t.\t{}\tGT:DP:AD:HF\t.:{}:{},{}:{:.6}\n",
                sanitize_vcf_token(&variant.chrom),
                variant.position,
                event_id,
                sanitize_vcf_token(&variant.reference),
                sanitize_vcf_token(&variant.alternate),
                info.join(";"),
                variant.depth,
                variant.reference_support,
                variant.support,
                variant.heteroplasmy()
            ));
        } else {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t.\tPASS\t{}\tGT:HF\t0/1:{:.6}\n",
                sanitize_vcf_token(&variant.chrom),
                variant.position,
                event_id,
                sanitize_vcf_token(&variant.reference),
                sanitize_vcf_token(&variant.alternate),
                info.join(";"),
                variant.heteroplasmy()
            ));
        }
    }

    Ok(out)
}

fn validate_result_contract(data: &Value) -> Result<()> {
    require_object(data, "/metadata")?;
    require_object(data, "/filter_stats")?;
    let schema_version = require_string(data, "/metadata/schema_version")?;
    if schema_version != "0.5" && schema_version != "0.6" {
        anyhow::bail!("result contract violation: unsupported schema {schema_version}");
    }
    require_string(data, "/metadata/sv_event_schema_version")?;
    require_string(data, "/metadata/complex_sv_event_schema_version")?;
    require_string(data, "/metadata/clinical_annotation_schema_version")?;
    require_string(data, "/metadata/engine_version")?;
    require_u64(data, "/metadata/reference_length")?;
    require_array(data, "/reads")?;
    require_array(data, "/variants")?;
    require_array(data, "/svs")?;
    require_array(data, "/complex_events")?;
    require_array(data, "/clusters")?;
    require_array(data, "/metadata/resources")?;
    require_object(data, "/filter_stats/numt_assessment")?;
    require_u64(data, "/filter_stats/passed_reads")?;
    if schema_version == "0.6" {
        require_object(data, "/evidence_encoding")?;
        require_array(data, "/alignments")?;
        require_array(data, "/molecules")?;
        require_array(data, "/callability")?;
        require_array(data, "/events")?;
        require_array(data, "/observation_pages")?;
        require_array(data, "/phase_links")?;
        require_array(data, "/architectures")?;
        validate_schema_0_6_contract(data)?;
    }
    Ok(())
}

fn validate_schema_0_6_contract(data: &Value) -> Result<()> {
    if data.get("observations").is_some() {
        anyhow::bail!("schema 0.6 row observations are forbidden; use observation_pages");
    }
    let reference_length = require_u64(data, "/metadata/reference_length")?;
    let encoding = require_object(data, "/evidence_encoding")?;
    if encoding.get("layout").and_then(Value::as_str) != Some("paged_columnar_molecule_event")
        || encoding.get("observation_storage").and_then(Value::as_str)
            != Some("embedded_columnar_pages")
        || encoding.get("missing_pair_state").and_then(Value::as_str) != Some("NOT_CALLABLE")
        || encoding
            .get("phase_molecule_policy")
            .and_then(Value::as_str)
            != Some("evidence_eligible_only")
        || encoding
            .get("phase_molecule_reference")
            .and_then(Value::as_str)
            != Some("molecules[].index")
        || encoding.get("phase_null_model").and_then(Value::as_str)
            != Some("independent_marginals_within_jointly_callable")
    {
        anyhow::bail!("schema 0.6 evidence encoding is unsupported");
    }
    let observation_limit = encoding
        .get("observation_limit")
        .and_then(Value::as_u64)
        .context("schema 0.6 observation_limit must be an integer")?;
    let phase_link_limit = encoding
        .get("phase_link_limit")
        .and_then(Value::as_u64)
        .context("schema 0.6 phase_link_limit must be an integer")?;
    let page_size = encoding
        .get("observation_page_size")
        .and_then(Value::as_u64)
        .context("schema 0.6 observation_page_size must be an integer")?;
    if page_size == 0 || page_size > 1_000_000 {
        anyhow::bail!("schema 0.6 observation_page_size is outside the supported envelope");
    }

    let alignments = require_array(data, "/alignments")?;
    let molecules = require_array(data, "/molecules")?;
    let callability = require_array(data, "/callability")?;
    let events = require_array(data, "/events")?;
    let observations = materialize_observations(data)?;
    let phase_links = require_array(data, "/phase_links")?;
    if observations.len() as u64 > observation_limit || phase_links.len() as u64 > phase_link_limit
    {
        anyhow::bail!("schema 0.6 evidence payload exceeds its declared resource limit");
    }
    if molecules.len() != callability.len() {
        anyhow::bail!("schema 0.6 callability/molecule cardinality mismatch");
    }

    let alignment_ids = alignments
        .iter()
        .enumerate()
        .map(|(index, alignment)| {
            alignment
                .get("id")
                .and_then(Value::as_str)
                .map(str::to_owned)
                .with_context(|| format!("/alignments/{index}/id must be a string"))
        })
        .collect::<Result<BTreeSet<_>>>()?;
    if alignment_ids.len() != alignments.len() {
        anyhow::bail!("schema 0.6 alignment IDs are not unique");
    }
    let mut alignment_owners = BTreeMap::new();
    for (index, alignment) in alignments.iter().enumerate() {
        let alignment_id = alignment
            .get("id")
            .and_then(Value::as_str)
            .with_context(|| format!("/alignments/{index}/id must be a string"))?;
        let molecule_id = alignment
            .get("molecule_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/alignments/{index}/molecule_id must be a string"))?;
        alignment_owners.insert(alignment_id, molecule_id);
    }
    let molecule_ids = molecules
        .iter()
        .enumerate()
        .map(|(index, molecule)| {
            molecule
                .get("id")
                .and_then(Value::as_str)
                .map(str::to_owned)
                .with_context(|| format!("/molecules/{index}/id must be a string"))
        })
        .collect::<Result<BTreeSet<_>>>()?;
    if molecule_ids.len() != molecules.len() {
        anyhow::bail!("schema 0.6 molecule IDs are not unique");
    }
    let mut projected_alignments = BTreeSet::new();
    let mut phase_eligible_molecule_ids = BTreeSet::new();
    let mut phase_eligible_molecules = BTreeMap::<usize, String>::new();
    for (index, molecule) in molecules.iter().enumerate() {
        let molecule_id = molecule
            .get("id")
            .and_then(Value::as_str)
            .with_context(|| format!("/molecules/{index}/id must be a string"))?;
        if molecule.get("index").and_then(Value::as_u64) != Some(index as u64) {
            anyhow::bail!("schema 0.6 molecule index is not contiguous at {molecule_id}");
        }
        let representative = molecule
            .get("representative_alignment_id")
            .and_then(Value::as_str)
            .with_context(|| {
                format!("/molecules/{index}/representative_alignment_id must be a string")
            })?;
        let identity_policy = molecule
            .get("identity_policy")
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
            .with_context(|| format!("/molecules/{index}/identity_policy must be non-empty"))?;
        let source_qnames = molecule
            .get("source_qnames")
            .and_then(Value::as_array)
            .with_context(|| format!("/molecules/{index}/source_qnames must be an array"))?;
        let source_qname_set = source_qnames
            .iter()
            .map(|value| value.as_str().context("source QNAME must be a string"))
            .collect::<Result<BTreeSet<_>>>()?;
        if source_qname_set.is_empty() || source_qname_set.len() != source_qnames.len() {
            anyhow::bail!("schema 0.6 molecule {molecule_id} has invalid source QNAMEs");
        }
        let protocol_metadata = molecule
            .get("protocol_metadata")
            .and_then(Value::as_object)
            .with_context(|| format!("/molecules/{index}/protocol_metadata must be an object"))?;
        if protocol_metadata
            .values()
            .any(|value| value.as_str().is_none())
        {
            anyhow::bail!("schema 0.6 molecule {molecule_id} protocol metadata must be strings");
        }
        if identity_policy.starts_with("sam_tag:")
            && protocol_metadata
                .get("molecule_id_tag")
                .and_then(Value::as_str)
                .is_none()
        {
            anyhow::bail!("schema 0.6 tagged molecule {molecule_id} lacks tag provenance");
        }
        let analysis_eligible = molecule
            .get("analysis_eligible")
            .and_then(Value::as_bool)
            .with_context(|| format!("/molecules/{index}/analysis_eligible must be boolean"))?;
        let evidence_eligible = molecule
            .get("evidence_eligible")
            .and_then(Value::as_bool)
            .with_context(|| format!("/molecules/{index}/evidence_eligible must be boolean"))?;
        if evidence_eligible && !analysis_eligible {
            anyhow::bail!(
                "schema 0.6 molecule {molecule_id} cannot be evidence-eligible while analysis-ineligible"
            );
        }
        if evidence_eligible {
            phase_eligible_molecule_ids.insert(molecule_id.to_owned());
            phase_eligible_molecules.insert(index, molecule_id.to_owned());
        }
        let exclusion_reasons = molecule
            .get("exclusion_reasons")
            .and_then(Value::as_array)
            .with_context(|| format!("/molecules/{index}/exclusion_reasons must be an array"))?;
        let exclusion_set = exclusion_reasons
            .iter()
            .map(|value| value.as_str().context("exclusion reason must be a string"))
            .collect::<Result<BTreeSet<_>>>()?;
        if exclusion_set.len() != exclusion_reasons.len()
            || (analysis_eligible && !exclusion_set.is_empty())
            || (!analysis_eligible && exclusion_set.is_empty())
        {
            anyhow::bail!("schema 0.6 molecule {molecule_id} exclusion contract is inconsistent");
        }
        let mut molecule_alignments = BTreeSet::new();
        for alignment_id in molecule
            .get("alignment_ids")
            .and_then(Value::as_array)
            .with_context(|| format!("/molecules/{index}/alignment_ids must be an array"))?
        {
            let alignment_id = alignment_id
                .as_str()
                .context("molecule alignment ID must be a string")?;
            if alignment_owners.get(alignment_id).copied() != Some(molecule_id)
                || !molecule_alignments.insert(alignment_id)
                || !projected_alignments.insert(alignment_id)
            {
                anyhow::bail!(
                    "schema 0.6 molecule/alignment ownership mismatch for {alignment_id}"
                );
            }
        }
        if !molecule_alignments.contains(representative) {
            anyhow::bail!("schema 0.6 representative alignment is unresolved for {molecule_id}");
        }
    }
    if projected_alignments != alignment_ids.iter().map(String::as_str).collect() {
        anyhow::bail!("schema 0.6 not every alignment resolves to exactly one molecule");
    }
    let mut callability_molecules = BTreeSet::new();
    for (index, summary) in callability.iter().enumerate() {
        let molecule_id = summary
            .get("molecule_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/callability/{index}/molecule_id must be a string"))?;
        if !molecule_ids.contains(molecule_id) {
            anyhow::bail!("schema 0.6 callability references unknown molecule {molecule_id}");
        }
        if !callability_molecules.insert(molecule_id) {
            anyhow::bail!("schema 0.6 callability duplicates molecule {molecule_id}");
        }
        let ranges = summary
            .get("ranges")
            .and_then(Value::as_array)
            .with_context(|| format!("/callability/{index}/ranges must be an array"))?;
        let mut previous_end = 0_u64;
        for (range_index, range) in ranges.iter().enumerate() {
            let start = range
                .get("start")
                .and_then(Value::as_u64)
                .with_context(|| {
                    format!("/callability/{index}/ranges/{range_index}/start must be an integer")
                })?;
            let end = range.get("end").and_then(Value::as_u64).with_context(|| {
                format!("/callability/{index}/ranges/{range_index}/end must be an integer")
            })?;
            if start == 0 || end < start || end > reference_length || start <= previous_end {
                anyhow::bail!(
                    "schema 0.6 callability ranges are invalid at molecule {molecule_id}"
                );
            }
            previous_end = end;
        }
        for (alignment_index, alignment) in summary
            .get("alignments")
            .and_then(Value::as_array)
            .with_context(|| format!("/callability/{index}/alignments must be an array"))?
            .iter()
            .enumerate()
        {
            let alignment_id = alignment
                .get("alignment_id")
                .and_then(Value::as_str)
                .with_context(|| format!("/callability/{index}/alignments/{alignment_index}/alignment_id must be a string"))?;
            if alignment_owners.get(alignment_id).copied() != Some(molecule_id) {
                anyhow::bail!(
                    "schema 0.6 callability alignment {alignment_id} belongs to another molecule"
                );
            }
        }
    }

    let mut event_ids = BTreeSet::new();
    let mut event_order = BTreeMap::<String, usize>::new();
    let mut event_complete_callability = BTreeMap::<String, bool>::new();
    let mut variant_events = BTreeMap::<String, (String, u64, u64, String, String)>::new();
    for (index, event) in events.iter().enumerate() {
        let event_id = event
            .get("id")
            .and_then(Value::as_str)
            .with_context(|| format!("/events/{index}/id must be a string"))?;
        if event.get("index").and_then(Value::as_u64) != Some(index as u64) {
            anyhow::bail!("schema 0.6 event index is not contiguous at {event_id}");
        }
        if !event_ids.insert(event_id.to_owned()) {
            anyhow::bail!("schema 0.6 event ID is duplicated: {event_id}");
        }
        event_order.insert(event_id.to_owned(), index);
        let event_type = event
            .get("type")
            .and_then(Value::as_str)
            .with_context(|| format!("/events/{index}/type must be a string"))?;
        if matches!(event_type, "SNV" | "SMALL_INSERTION" | "SMALL_DELETION") {
            variant_events.insert(
                event_id.to_owned(),
                (
                    event_type.to_owned(),
                    event
                        .get("start")
                        .and_then(Value::as_u64)
                        .with_context(|| {
                            format!("/events/{index}/start must be an integer for a variant event")
                        })?,
                    event.get("end").and_then(Value::as_u64).with_context(|| {
                        format!("/events/{index}/end must be an integer for a variant event")
                    })?,
                    event
                        .get("ref")
                        .and_then(Value::as_str)
                        .with_context(|| format!("/events/{index}/ref must be a string"))?
                        .to_owned(),
                    event
                        .get("alt")
                        .and_then(Value::as_str)
                        .with_context(|| format!("/events/{index}/alt must be a string"))?
                        .to_owned(),
                ),
            );
        }
        let assessability = event
            .get("assessability")
            .and_then(Value::as_str)
            .with_context(|| format!("/events/{index}/assessability must be a string"))?;
        let negative_rule = event
            .get("negative_evidence_rule")
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
            .with_context(|| format!("/events/{index}/negative_evidence_rule must be non-empty"))?;
        event_complete_callability.insert(
            event_id.to_owned(),
            assessability == "REFERENCE_AND_ALTERNATE",
        );
        if (assessability == "ALTERNATE_SUPPORT_ONLY")
            != (negative_rule == "support_only_no_negative_inference")
        {
            anyhow::bail!("schema 0.6 event {event_id} has an inconsistent negative-evidence rule");
        }
        for molecule_id in event
            .get("supporting_molecule_ids")
            .and_then(Value::as_array)
            .with_context(|| format!("/events/{index}/supporting_molecule_ids must be an array"))?
        {
            let molecule_id = molecule_id
                .as_str()
                .context("supporting molecule ID must be a string")?;
            if !molecule_ids.contains(molecule_id) {
                anyhow::bail!(
                    "schema 0.6 event {event_id} references unknown molecule {molecule_id}"
                );
            }
        }
    }

    let mut molecule_event_pairs = BTreeSet::new();
    let mut projected_support: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut projected_counts: BTreeMap<String, [u64; 5]> = BTreeMap::new();
    let mut projected_variant_counts: BTreeMap<String, [u64; 6]> = BTreeMap::new();
    let mut projected_molecule_counts: BTreeMap<String, [u64; 5]> = BTreeMap::new();
    let mut projected_molecule_alternates: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut phase_observations = BTreeMap::<String, BTreeMap<String, String>>::new();
    for (index, observation) in observations.iter().enumerate() {
        let expected_id = format!("observation:{index}");
        if observation.get("id").and_then(Value::as_str) != Some(expected_id.as_str()) {
            anyhow::bail!("schema 0.6 observation ID is not contiguous at index {index}");
        }
        let molecule_id = observation
            .get("molecule_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/observations/{index}/molecule_id must be a string"))?;
        let event_id = observation
            .get("event_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/observations/{index}/event_id must be a string"))?;
        let alignment_id = observation
            .get("alignment_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/observations/{index}/alignment_id must be a string"))?;
        let state = observation
            .get("state")
            .and_then(Value::as_str)
            .with_context(|| format!("/observations/{index}/state must be a string"))?;
        if !molecule_ids.contains(molecule_id)
            || !event_ids.contains(event_id)
            || !alignment_ids.contains(alignment_id)
        {
            anyhow::bail!("schema 0.6 observation {index} has an unresolved reference");
        }
        if !phase_eligible_molecule_ids.contains(molecule_id) {
            anyhow::bail!(
                "schema 0.6 observation {index} belongs to a molecule excluded from evidence"
            );
        }
        if !matches!(
            state,
            "REFERENCE" | "ALTERNATE" | "EVENT_ABSENT" | "LOW_QUALITY" | "CONFLICT"
        ) {
            anyhow::bail!("schema 0.6 observation {index} has unsupported state {state}");
        }
        let evidence_source = observation
            .get("evidence_source")
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
            .with_context(|| format!("/observations/{index}/evidence_source must be non-empty"))?;
        if !molecule_event_pairs.insert((molecule_id.to_owned(), event_id.to_owned())) {
            anyhow::bail!("schema 0.6 molecule/event observation is duplicated");
        }
        phase_observations
            .entry(molecule_id.to_owned())
            .or_default()
            .insert(event_id.to_owned(), state.to_owned());
        if state == "ALTERNATE" {
            projected_support
                .entry(event_id.to_owned())
                .or_default()
                .insert(molecule_id.to_owned());
        }
        let count_index = match state {
            "ALTERNATE" => 0,
            "REFERENCE" => 1,
            "EVENT_ABSENT" => 2,
            "LOW_QUALITY" => 3,
            "CONFLICT" => 4,
            _ => unreachable!("observation state was validated above"),
        };
        let count = &mut projected_counts.entry(event_id.to_owned()).or_default()[count_index];
        *count = count
            .checked_add(1)
            .context("schema 0.6 observation count overflow")?;
        let molecule_count = &mut projected_molecule_counts
            .entry(molecule_id.to_owned())
            .or_default()[count_index];
        *molecule_count = molecule_count
            .checked_add(1)
            .context("schema 0.6 molecule observation count overflow")?;
        if state == "ALTERNATE" {
            projected_molecule_alternates
                .entry(molecule_id.to_owned())
                .or_default()
                .insert(event_id.to_owned());
        }
        if let Some((event_type, _, _, _, _)) = variant_events.get(event_id) {
            let counts = projected_variant_counts
                .entry(event_id.to_owned())
                .or_default();
            match state {
                "ALTERNATE" => counts[0] += 1,
                "REFERENCE" => counts[1] += 1,
                "EVENT_ABSENT" => {
                    counts[3] += 1;
                    if event_type == "SNV" || evidence_source == "cigar_alternative_small_indel" {
                        counts[2] += 1;
                    } else {
                        counts[1] += 1;
                    }
                }
                "LOW_QUALITY" => counts[4] += 1,
                "CONFLICT" => counts[5] += 1,
                _ => unreachable!("observation state was validated above"),
            }
        }
    }
    let mut verified_event_support = BTreeMap::<String, BTreeSet<String>>::new();
    for (index, event) in events.iter().enumerate() {
        let event_id = event.get("id").and_then(Value::as_str).unwrap_or_default();
        let declared = event
            .get("supporting_molecule_ids")
            .and_then(Value::as_array)
            .with_context(|| format!("/events/{index}/supporting_molecule_ids must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_str()
                    .map(str::to_owned)
                    .context("supporting molecule ID must be a string")
            })
            .collect::<Result<BTreeSet<_>>>()?;
        if declared != projected_support.remove(event_id).unwrap_or_default() {
            anyhow::bail!("schema 0.6 support projection mismatch for {event_id}");
        }
        verified_event_support.insert(event_id.to_owned(), declared);
        let counts = projected_counts.remove(event_id).unwrap_or_default();
        let declared_counts = event
            .get("evidence_counts")
            .and_then(Value::as_object)
            .with_context(|| format!("/events/{index}/evidence_counts must be an object"))?;
        let expected = [
            ("alternate", counts[0]),
            ("reference", counts[1]),
            ("event_absent", counts[2]),
            ("low_quality", counts[3]),
            ("conflict", counts[4]),
            ("callable", counts[0] + counts[1] + counts[2]),
        ];
        for (field, expected_count) in expected {
            if declared_counts.get(field).and_then(Value::as_u64) != Some(expected_count) {
                anyhow::bail!("schema 0.6 evidence count mismatch for {event_id}/{field}");
            }
        }
    }

    for (index, molecule) in molecules.iter().enumerate() {
        let molecule_id = molecule
            .get("id")
            .and_then(Value::as_str)
            .with_context(|| format!("/molecules/{index}/id must be a string"))?;
        let projected = projected_molecule_counts
            .remove(molecule_id)
            .unwrap_or_default();
        let declared = molecule
            .get("evidence_state_counts")
            .and_then(Value::as_object)
            .with_context(|| {
                format!("/molecules/{index}/evidence_state_counts must be an object")
            })?;
        for (field, expected) in [
            ("alternate", projected[0]),
            ("reference", projected[1]),
            ("event_absent", projected[2]),
            ("low_quality", projected[3]),
            ("conflict", projected[4]),
        ] {
            if declared.get(field).and_then(Value::as_u64) != Some(expected) {
                anyhow::bail!(
                    "schema 0.6 molecule evidence count mismatch for {molecule_id}/{field}"
                );
            }
        }
        let alternate_event_ids = molecule
            .get("alternate_event_ids")
            .and_then(Value::as_array)
            .with_context(|| format!("/molecules/{index}/alternate_event_ids must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_str()
                    .map(str::to_owned)
                    .context("molecule alternate event ID must be a string")
            })
            .collect::<Result<BTreeSet<_>>>()?;
        if alternate_event_ids
            != projected_molecule_alternates
                .remove(molecule_id)
                .unwrap_or_default()
        {
            anyhow::bail!(
                "schema 0.6 molecule alternate-event projection mismatch for {molecule_id}"
            );
        }
        for field in ["query_length", "mapping_quality"] {
            if molecule.get(field).and_then(Value::as_u64).is_none() {
                anyhow::bail!("schema 0.6 molecule {molecule_id}/{field} must be an integer");
            }
        }
        for field in ["mean_base_quality", "numt_score"] {
            if molecule.get(field).and_then(Value::as_f64).is_none() {
                anyhow::bail!("schema 0.6 molecule {molecule_id}/{field} must be a number");
            }
        }
        require_array(molecule, "/numt_evidence")?;
    }

    let variants = require_array(data, "/variants")?;
    if variants.len() != variant_events.len() {
        anyhow::bail!("schema 0.6 unified variant/event cardinality mismatch");
    }
    let mut projected_variant_ids = BTreeSet::new();
    for (index, variant) in variants.iter().enumerate() {
        let event_id = variant
            .get("event_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/variants/{index}/event_id must be a string"))?;
        if !projected_variant_ids.insert(event_id.to_owned()) {
            anyhow::bail!("schema 0.6 variant event ID is duplicated: {event_id}");
        }
        let (event_type, start, end, reference, alternate) = variant_events
            .get(event_id)
            .with_context(|| format!("variant {event_id} has no normalized source event"))?;
        if variant.get("type").and_then(Value::as_str) != Some(event_type.as_str())
            || variant.get("position").and_then(Value::as_u64) != Some(*start)
            || variant.get("start").and_then(Value::as_u64) != Some(*start)
            || variant.get("end").and_then(Value::as_u64) != Some(*end)
            || variant.get("ref").and_then(Value::as_str) != Some(reference.as_str())
            || variant.get("alt").and_then(Value::as_str) != Some(alternate.as_str())
        {
            anyhow::bail!("schema 0.6 variant/event identity mismatch for {event_id}");
        }
        let counts = projected_variant_counts
            .get(event_id)
            .copied()
            .unwrap_or_default();
        let callable = counts[0]
            .checked_add(counts[1])
            .and_then(|value| value.checked_add(counts[2]))
            .context("schema 0.6 variant callable count overflow")?;
        for (field, expected) in [
            ("alt_depth", counts[0]),
            ("ref_depth", counts[1]),
            ("other_depth", counts[2]),
            ("event_absent_depth", counts[3]),
            ("low_quality_depth", counts[4]),
            ("conflict_depth", counts[5]),
            ("callable_depth", callable),
        ] {
            if variant.get(field).and_then(Value::as_u64) != Some(expected) {
                anyhow::bail!("schema 0.6 variant count mismatch for {event_id}/{field}");
            }
        }
        let declared_support = variant
            .get("supporting_molecule_ids")
            .and_then(Value::as_array)
            .with_context(|| format!("/variants/{index}/supporting_molecule_ids must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_str()
                    .map(str::to_owned)
                    .context("variant supporting molecule ID must be a string")
            })
            .collect::<Result<BTreeSet<_>>>()?;
        if verified_event_support.get(event_id) != Some(&declared_support) {
            anyhow::bail!("schema 0.6 variant support mismatch for {event_id}");
        }
        let heteroplasmy = variant
            .get("heteroplasmy")
            .and_then(Value::as_f64)
            .with_context(|| format!("/variants/{index}/heteroplasmy must be a number"))?;
        let expected_hf = if callable == 0 {
            0.0
        } else {
            counts[0] as f64 / callable as f64
        };
        if (heteroplasmy - expected_hf).abs() > 2e-9 {
            anyhow::bail!("schema 0.6 variant HF mismatch for {event_id}");
        }
        require_string(variant, "/filter_status")?;
        require_array(variant, "/qc_flags")?;
        require_string(variant, "/numt_assessability")?;
    }
    if projected_variant_ids != variant_events.keys().cloned().collect() {
        anyhow::bail!("schema 0.6 variant projection does not cover every SNV/small-indel event");
    }
    let mut expected_phase_pairs = BTreeSet::<(String, String)>::new();
    for observations in phase_observations.values() {
        for (alternate_event, state) in observations {
            if state != "ALTERNATE" {
                continue;
            }
            for other_event in observations.keys() {
                if other_event == alternate_event {
                    continue;
                }
                let alternate_order = *event_order
                    .get(alternate_event)
                    .context("phase alternate event order is unresolved")?;
                let other_order = *event_order
                    .get(other_event)
                    .context("phase neighbor event order is unresolved")?;
                let pair = if alternate_order < other_order {
                    (alternate_event.clone(), other_event.clone())
                } else {
                    (other_event.clone(), alternate_event.clone())
                };
                expected_phase_pairs.insert(pair);
            }
        }
    }
    let mut phase_ids = BTreeSet::new();
    for (index, link) in phase_links.iter().enumerate() {
        let phase_id = link
            .get("id")
            .and_then(Value::as_str)
            .with_context(|| format!("/phase_links/{index}/id must be a string"))?;
        if !phase_ids.insert(phase_id) {
            anyhow::bail!("schema 0.6 phase link ID is duplicated: {phase_id}");
        }
        let event_a = link
            .get("event_a_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/phase_links/{index}/event_a_id must be a string"))?;
        let event_b = link
            .get("event_b_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/phase_links/{index}/event_b_id must be a string"))?;
        if !event_ids.contains(event_a) || !event_ids.contains(event_b) || event_a == event_b {
            anyhow::bail!("schema 0.6 phase link {index} has invalid event references");
        }
        if event_order.get(event_a) >= event_order.get(event_b)
            || phase_id != format!("phase:{event_a}|{event_b}")
        {
            anyhow::bail!("schema 0.6 phase link {index} is not in canonical event order");
        }
        if !expected_phase_pairs.remove(&(event_a.to_owned(), event_b.to_owned())) {
            anyhow::bail!(
                "schema 0.6 phase link {index} is not derived from an alternate observation"
            );
        }
        let complete = event_complete_callability.get(event_a) == Some(&true)
            && event_complete_callability.get(event_b) == Some(&true);
        let expected_assessability = if complete {
            "COMPLETE_FOR_BOTH_EVENTS"
        } else {
            "SUPPORT_CONDITIONED"
        };
        if link.get("assessability").and_then(Value::as_str) != Some(expected_assessability) {
            anyhow::bail!("schema 0.6 phase assessability mismatch at link {index}");
        }

        let mut counts = [0_u64; 5];
        let mut supporting_molecule_indices = Vec::new();
        let mut uncertain_molecule_indices = Vec::new();
        for (molecule_index, molecule_id) in &phase_eligible_molecules {
            let Some(observations) = phase_observations.get(molecule_id) else {
                continue;
            };
            let (Some(state_a), Some(state_b)) =
                (observations.get(event_a), observations.get(event_b))
            else {
                continue;
            };
            let callable_a = matches!(state_a.as_str(), "REFERENCE" | "ALTERNATE" | "EVENT_ABSENT");
            let callable_b = matches!(state_b.as_str(), "REFERENCE" | "ALTERNATE" | "EVENT_ABSENT");
            if !callable_a || !callable_b {
                counts[4] = counts[4]
                    .checked_add(1)
                    .context("phase uncertain count overflow")?;
                uncertain_molecule_indices.push(*molecule_index as u64);
                continue;
            }
            match (
                state_a.as_str() == "ALTERNATE",
                state_b.as_str() == "ALTERNATE",
            ) {
                (true, true) => {
                    counts[0] += 1;
                    supporting_molecule_indices.push(*molecule_index as u64);
                }
                (true, false) => counts[1] += 1,
                (false, true) => counts[2] += 1,
                (false, false) => counts[3] += 1,
            }
        }
        let jointly_callable = counts[0]
            .checked_add(counts[1])
            .and_then(|value| value.checked_add(counts[2]))
            .and_then(|value| value.checked_add(counts[3]))
            .context("phase jointly-callable count overflow")?;
        for (field, expected_count) in [
            ("both_alternate", counts[0]),
            ("a_alternate_b_absent", counts[1]),
            ("a_absent_b_alternate", counts[2]),
            ("neither_alternate", counts[3]),
            ("jointly_uncertain", counts[4]),
            ("jointly_callable", jointly_callable),
        ] {
            if link.get(field).and_then(Value::as_u64) != Some(expected_count) {
                anyhow::bail!(
                    "schema 0.6 phase observation projection mismatch at link {index}/{field}"
                );
            }
        }
        let declared_support = phase_u64_array(link, "supporting_molecule_indices", index)?;
        let declared_uncertain = phase_u64_array(link, "uncertain_molecule_indices", index)?;
        if declared_support != supporting_molecule_indices
            || declared_uncertain != uncertain_molecule_indices
        {
            anyhow::bail!("schema 0.6 phase molecule traceability mismatch at link {index}");
        }
        let expected_qc_flags = [
            (!complete).then_some("SUPPORT_CONDITIONED"),
            (counts[4] != 0).then_some("UNCERTAIN_COOCCURRENCE_EXCLUDED"),
        ]
        .into_iter()
        .flatten()
        .map(str::to_owned)
        .collect::<Vec<_>>();
        if phase_string_array(link, "qc_flags", index)? != expected_qc_flags {
            anyhow::bail!("schema 0.6 phase QC facts mismatch at link {index}");
        }

        let observed = if jointly_callable == 0 {
            0.0
        } else {
            counts[0] as f64 / jointly_callable as f64
        };
        let expected = if jointly_callable == 0 {
            0.0
        } else {
            let n = jointly_callable as f64;
            ((counts[0] + counts[1]) as f64 / n) * ((counts[0] + counts[2]) as f64 / n)
        };
        let (ci_low, ci_high) = phase_wilson_interval(counts[0], jointly_callable);
        for (field, expected_value) in [
            ("co_alternate_fraction", observed),
            ("co_alternate_ci95_low", ci_low),
            ("co_alternate_ci95_high", ci_high),
            ("expected_co_alternate_fraction", expected),
            ("linkage_delta", observed - expected),
        ] {
            let actual = link
                .get(field)
                .and_then(Value::as_f64)
                .with_context(|| format!("phase {field} must be a number"))?;
            if (actual - expected_value).abs() > 2e-9 {
                anyhow::bail!("schema 0.6 phase statistic mismatch at link {index}/{field}");
            }
        }
    }
    if !expected_phase_pairs.is_empty() {
        anyhow::bail!("schema 0.6 phase projection omits candidate event pairs");
    }
    Ok(())
}

fn phase_string_array(link: &Value, field: &str, index: usize) -> Result<Vec<String>> {
    link.get(field)
        .and_then(Value::as_array)
        .with_context(|| format!("/phase_links/{index}/{field} must be an array"))?
        .iter()
        .map(|value| {
            value
                .as_str()
                .map(str::to_owned)
                .with_context(|| format!("/phase_links/{index}/{field} values must be strings"))
        })
        .collect()
}

fn phase_u64_array(link: &Value, field: &str, index: usize) -> Result<Vec<u64>> {
    link.get(field)
        .and_then(Value::as_array)
        .with_context(|| format!("/phase_links/{index}/{field} must be an array"))?
        .iter()
        .map(|value| {
            value
                .as_u64()
                .with_context(|| format!("/phase_links/{index}/{field} values must be integers"))
        })
        .collect()
}

fn phase_wilson_interval(successes: u64, total: u64) -> (f64, f64) {
    if total == 0 {
        return (0.0, 0.0);
    }
    const Z: f64 = 1.959_963_984_540_054;
    let n = total as f64;
    let p = successes as f64 / n;
    let z2 = Z * Z;
    let denominator = 1.0 + z2 / n;
    let center = (p + z2 / (2.0 * n)) / denominator;
    let margin = Z * ((p * (1.0 - p) / n) + z2 / (4.0 * n * n)).sqrt() / denominator;
    ((center - margin).max(0.0), (center + margin).min(1.0))
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
    use super::{
        assert_aux_present, assert_mapq_present, assert_sv_present, materialize_observations,
        render_variant_tsv, render_vcf,
    };

    #[test]
    fn renders_vcf_with_locus_callable_heteroplasmy() {
        let json = r#"{
          "metadata": {
            "schema_version": "0.5",
            "sv_event_schema_version": "1.0",
            "complex_sv_event_schema_version": "1.0",
            "clinical_annotation_schema_version": "1.0",
            "engine_version": "0.5.0-dev",
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
          "complex_events": [],
          "clusters": []
        }"#;

        let vcf = render_vcf(json).expect("VCF should render");
        assert!(vcf.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A"));
        assert!(vcf.contains("NC_012920.1\t3243\t.\tA\tG\t.\tPASS\tAC=2;DP=3;HF=0.666667;GENE=MT-TL1;CLNSIG=pathogenic;HF_CI95=0.207660,0.938508\tGT:HF\t0/1:0.666667"));
    }

    #[test]
    fn renders_schema_0_6_small_indel_with_stable_event_id() {
        let json = r#"{
          "metadata": {
            "schema_version": "0.6",
            "sv_event_schema_version": "1.0",
            "complex_sv_event_schema_version": "1.0",
            "clinical_annotation_schema_version": "1.0",
            "engine_version": "0.5.0-dev",
            "reference_length": 16569,
            "sample": "phase",
            "reference_accession": "NC_012920.1",
            "resources": []
          },
          "filter_stats": {
            "passed_reads": 2,
            "numt_assessment": { "mode": "mt_only_or_unknown" }
          },
          "evidence_encoding": {
            "layout": "paged_columnar_molecule_event",
            "observation_storage": "embedded_columnar_pages",
            "missing_pair_state": "NOT_CALLABLE",
            "phase_molecule_policy": "evidence_eligible_only",
            "phase_molecule_reference": "molecules[].index",
            "phase_null_model": "independent_marginals_within_jointly_callable",
            "observation_limit": 10,
            "observation_count": 2,
            "observation_page_size": 2,
            "observation_page_count": 1,
            "phase_link_limit": 10
          },
          "alignments": [
            {"id":"alignment:0","molecule_id":"m1"},
            {"id":"alignment:1","molecule_id":"m2"}
          ],
          "molecules": [
            {"id":"m1","index":0,"identity_policy":"sam_qname","source_qnames":["m1"],"protocol_metadata":{},"analysis_eligible":true,"evidence_eligible":true,"exclusion_reasons":[],"representative_alignment_id":"alignment:0","alignment_ids":["alignment:0"],"query_length":12,"mean_base_quality":40,"mapping_quality":60,"numt_score":0,"numt_evidence":[],"alternate_event_ids":["indel:insertion:5:G"],"evidence_state_counts":{"alternate":1,"reference":0,"event_absent":0,"low_quality":0,"conflict":0}},
            {"id":"m2","index":1,"identity_policy":"sam_qname","source_qnames":["m2"],"protocol_metadata":{},"analysis_eligible":true,"evidence_eligible":true,"exclusion_reasons":[],"representative_alignment_id":"alignment:1","alignment_ids":["alignment:1"],"query_length":12,"mean_base_quality":40,"mapping_quality":60,"numt_score":0,"numt_evidence":[],"alternate_event_ids":[],"evidence_state_counts":{"alternate":0,"reference":0,"event_absent":1,"low_quality":0,"conflict":0}}
          ],
          "callability": [
            {"molecule_id":"m1","ranges":[{"start":1,"end":12}],"alignments":[{"alignment_id":"alignment:0"}]},
            {"molecule_id":"m2","ranges":[{"start":1,"end":12}],"alignments":[{"alignment_id":"alignment:1"}]}
          ],
          "events": [{
            "id":"indel:insertion:5:G",
            "index":0,
            "type":"SMALL_INSERTION",
            "start":5,
            "end":5,
            "ref":"A",
            "alt":"AG",
            "negative_evidence_rule":"same_fragment_callable_reference_adjacency",
            "assessability":"REFERENCE_AND_ALTERNATE",
            "supporting_molecule_ids":["m1"],
            "evidence_counts":{"alternate":1,"reference":0,"event_absent":1,"callable":2,"low_quality":0,"conflict":0}
          }],
          "observation_pages": [
            {
              "index": 0,
              "offset": 0,
              "count": 2,
              "columns": {
                "molecule_id": ["m1", "m2"],
                "event_id": ["indel:insertion:5:G", "indel:insertion:5:G"],
                "alignment_id": ["alignment:0", "alignment:1"],
                "state": ["ALTERNATE", "EVENT_ABSENT"],
                "observed_allele": ["G", null],
                "base_quality": [40, null],
                "mapping_quality": [60, 60],
                "strand": ["+", "+"],
                "evidence_source": ["cigar_small_indel", "callable_reference_path"],
                "read_position": [0, null]
              }
            }
          ],
          "phase_links": [],
          "architectures": [],
          "variants": [{
            "event_id":"indel:insertion:5:G",
            "type":"SMALL_INSERTION",
            "position":5,
            "start":5,
            "end":5,
            "ref":"A",
            "alt":"AG",
            "normalization":"rcrs_3prime_small_indel_v1",
            "vcf_position":5,
            "vcf_representable":true,
            "alt_depth":1,
            "ref_depth":1,
            "other_depth":0,
            "event_absent_depth":1,
            "low_quality_depth":0,
            "conflict_depth":0,
            "callable_depth":2,
            "heteroplasmy":0.5,
            "ci95_low":0.094531,
            "ci95_high":0.905469,
            "filter_status":"NOT_CALIBRATED",
            "qc_flags":["NUMT_NOT_ASSESSABLE","SINGLE_STRAND_ALT_SUPPORT"],
            "numt_assessability":"NOT_ASSESSABLE",
            "supporting_molecule_ids":["m1"],
            "mapping_quality":{"alternate":{"mean":60}},
            "allele_quality":{"alternate":{"mean_phred":40}},
            "strand_support":{"alt_forward":1,"alt_reverse":0,"ref_forward":1,"ref_reverse":0}
          }],
          "reads": [],
          "svs": [],
          "complex_events": [],
          "clusters": []
        }"#;

        let vcf = render_vcf(json).expect("schema 0.6 VCF should render");
        assert!(vcf.contains("##INFO=<ID=EVENT_TYPE"));
        assert!(vcf.contains(
            "NC_012920.1\t5\tindel:insertion:5:G\tA\tAG\t.\t.\tAC=1;DP=2;HF=0.500000;AD=1,1;ODC=0;LOWQ=0;CONFLICT=0;MOLECULE_SUPPORT=1;MQ=60.0000;BQ=40.0000;STRAND_SUPPORT=1,0,1,0;NUMT_ASSESSABLE=NOT_ASSESSABLE;QC_FLAGS=NUMT_NOT_ASSESSABLE,SINGLE_STRAND_ALT_SUPPORT;NORMALIZATION=rcrs_3prime_small_indel_v1;HF_CI95=0.094531,0.905469;EVENT_TYPE=SMALL_INSERTION\tGT:DP:AD:HF\t.:2:1,1:0.500000"
        ));
        let tsv = render_variant_tsv(json).expect("schema 0.6 TSV should render");
        assert!(tsv.contains(
            "indel:insertion:5:G\tSMALL_INSERTION\t5\t5\t5\tA\tAG\t1\t1\t0\t2\t0.50000000"
        ));
    }

    #[test]
    fn materializes_columnar_observation_pages_and_rejects_cardinality_drift() {
        let mut value = serde_json::json!({
            "evidence_encoding": {
                "observation_count": 2,
                "observation_page_size": 2,
                "observation_page_count": 1
            },
            "observation_pages": [{
                "index": 0,
                "offset": 0,
                "count": 2,
                "columns": {
                    "molecule_id": ["m1", "m2"],
                    "event_id": ["e1", "e2"],
                    "alignment_id": ["alignment:0", "alignment:1"],
                    "state": ["ALTERNATE", "REFERENCE"],
                    "observed_allele": ["G", "A"],
                    "base_quality": [40, 39],
                    "mapping_quality": [60, 50],
                    "strand": ["+", "-"],
                    "evidence_source": ["aligned_base", "aligned_base"],
                    "read_position": [0, 1]
                }
            }]
        });
        let rows = materialize_observations(&value).expect("page should materialize");
        assert_eq!(rows.len(), 2);
        assert_eq!(
            rows[1].get("id").and_then(serde_json::Value::as_str),
            Some("observation:1")
        );

        value["observation_pages"][0]["columns"]["strand"] = serde_json::json!(["+"]);
        let error = materialize_observations(&value).expect_err("column drift must fail closed");
        assert!(error
            .to_string()
            .contains("column strand has the wrong length"));
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
