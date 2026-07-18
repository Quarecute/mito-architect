# Mito-Architect

Mito-Architect is a long-read platform for reconstructing the molecular
architecture of human mitochondrial DNA. It combines a C++ analysis core, a
Rust CLI/backend, a TypeScript visualization library, and a React web UI around
one result contract. The target instrument derives linked variant-, molecule-,
phasing-, and architecture-centric projections from one callable-aware evidence
graph.

The current `0.5.0-dev` build is an end-to-end research-use mtDNA analysis workbench: it
accepts FASTQ/SAM inputs directly, can use htslib for SAM/BAM/CRAM when htslib
is available at build time, uses a bundled rCRS reference (`NC_012920.1`),
reconstructs CIGAR and split-alignment events, clusters molecules by DBSCAN over
SNP/SV feature tokens, assigns ranked haplogroups from PhyloTree rCRS 17.3,
annotates known variants from a versioned local MITOMAP/ClinVar-compatible
cache, and renders interactive reports.

This repository is not yet an analytically validated clinical or operated
production bioinformatics service. The 0.3.x and 0.4 software milestones are
implemented, but blinded performance evidence, durable service infrastructure,
security controls, and laboratory sign-off remain mandatory for 1.0.

The scientific product identity and defensible novelty relative to
variant-centric and existing molecule-level tools are defined in
[`docs/SCIENTIFIC_POSITIONING.md`](docs/SCIENTIFIC_POSITIONING.md). The
executable production programme and release gates are defined in
[`docs/PRODUCTION_ACCEPTANCE.md`](docs/PRODUCTION_ACCEPTANCE.md). A release must
not be labelled production-ready until its external analytical and operational
evidence is complete. The ordered implementation milestones are in
[`docs/ROADMAP.md`](docs/ROADMAP.md).

Version 1.0 execution is active at **0.5 / RC1**. The live evidence ledger is
[`docs/RELEASE_STATUS.md`](docs/RELEASE_STATUS.md); the normative scientific
draft is [`docs/SCIENTIFIC_SPECIFICATION.md`](docs/SCIENTIFIC_SPECIFICATION.md),
and stable machine-readable failures are documented in
[`docs/ERROR_CONTRACT.md`](docs/ERROR_CONTRACT.md). Status is gate-based rather
than an unverifiable percentage.

For a map of user guides, scientific contracts, and verification evidence,
start with [`docs/README.md`](docs/README.md).
FASTQ ingestion does not perform read alignment. Alignment-derived SNPs,
CIGAR SVs, MAPQ, and coverage require SAM/BAM/CRAM input; raw FASTQ currently
contributes read-level quality/NUMT heuristics and explicit test tags only.

The current schema 0.5 `clusters` are deterministic development groups over
observed SNP/SV/complex-path tokens. They are not yet validated molecular
architectures. The opt-in schema 0.6 RC2 draft now introduces explicit alignment
fragments, assembly decisions, per-molecule callable ranges, normalized
SNV/small-indel/SV/complex-path events, bounded paged-columnar observations,
and callable-aware pairwise phase links with aggregate/evidence invariants.
Its draft contract is [`docs/SCHEMA_0.6_DRAFT.md`](docs/SCHEMA_0.6_DRAFT.md);
explicit SAM-tag molecule identity, audited UMI/duplex provenance, evidence
page resources, exact SNV/indel/SV phase traceability, a pinned mtDNA-Server 2
differential, remote-page virtualization, and RSS/storage/cancellation gates are
implemented development baselines. Wet-lab protocol review, laboratory UAT,
architecture inference, and independent scientific validation remain open
before freeze.
The exact local gate for the completed RC3 engineering slice is recorded in
[`docs/VERIFICATION_0.7_RC3_DEV.md`](docs/VERIFICATION_0.7_RC3_DEV.md).

## Highlights

- Bundled rCRS reference, length-checked at 16,569 bp.
- FASTA read input is rejected; FASTA is reference-only.
- FASTQ and SAM records are validated for required separators, field types,
  sequence/quality lengths, and CIGAR/query-length consistency. Every `SA:Z`
  entry is validated before filtering or parallel analysis; malformed split-
  alignment evidence fails closed with `MITO-E1103`.
- Optional htslib-backed SAM/BAM/CRAM reader with CIGAR, MAPQ, flags, and aux
  tag preservation.
- NUMT filtering hook with explicit read accounting.
- Competitive-alignment NUMT evidence and explicit warnings when specificity
  cannot be assessed.
- Molecule-level deletion, insertion, duplication, inversion, and circular-
  origin reconstruction from primary and supplementary alignments.
- Canonical SV event schema 1.0 merges equivalent CIGAR, forward-split, and
  reverse-split evidence with exact external JSON regression data.
- Exact external evidence goldens pin circular SNPs, multi-allelic depth,
  reverse-strand observations, missing quality, excluded flags, and molecule-
  conflict no-calls; a versioned negative manifest pins malformed FASTQ/CIGAR/
  `SA` failures.
- Versioned PhyloTree rCRS 17.3 haplogroup ranking with fail-closed
  substitution/deletion/insertion/wildcard/backmutation grammar, molecule-
  callable ranges, strict-majority cluster markers, alternatives, missing
  sites, quality, contamination warnings, and checksums. Twenty-eight pinned
  HaploGrep 3.3.2 profiles require exact top-three agreement across 23 lineage
  groups, including three inherited-marker removals by backmutation. The
  checksum-pinned official alignment-rule table is applied deterministically
  to compound cluster marker sets; five projection goldens cover 3-prime CIGAR
  normalization, motif rotation, SNP-to-deletion equivalence, a mixed
  insertion/SNP/deletion rewrite, and explicit reference restoration.
- CIGAR-aware coverage excludes unmapped and NUMT-filtered reads and reports
  exact site-level mean, maximum, and >20x statistics alongside display bins.
- Cluster drill-down: click a cluster to inspect its SNPs, SVs, multi-junction
  paths, support, frequencies, molecule IDs, and genome distribution.
- Clinical assertion schema 1.0 with per-source provenance, deterministic
  conflict/consensus, strict resource validation, JSON/VCF/UI drill-down, and a
  checksum-pinned MITOMAP/ClinVar-style development cache. Governed external
  snapshots support validate-before-publish staging, SHA-256 manifests,
  durable atomic activation, freshness checks, status, and verified rollback.
- Molecule-level alternate/reference/other support, forward/reverse evidence,
  strand/read-position observations, and per-allele Phred summaries. These QC
  fields are reported as evidence and do not apply uncalibrated hard filters.
- Opt-in schema 0.6 evidence graph with all alignment fragments, conservative
  base callability, primary/supplementary SNV evidence merging, normalized small
  indels and structural events, explicit support-only event states, bounded
  callable-aware phase partitions, exact co-ALT/uncertain molecule traceability,
  and a linked Molecules & phase UI.
- Pinned mtDNA-Server 2 v2.1.16 `variants.annotated.txt` differential with
  strict sample/coordinate/coverage validation, deterministic checksums and an
  explicit no-truth/no-clinical-validation interpretation boundary.
- Compact authenticated result summaries and pre-serialized immutable evidence
  pages; an interned server-side index provides exact bounded molecule/event/
  state search over all pages, while the web client virtualizes only the
  returned rows and keeps full JSON as an explicit authoritative download.
- A linked Rearrangements inspector traces ordered complex junction paths and
  canonical SVs through supporting molecules to retained alignment fragments;
  missing junctions, circular closure, and biological architecture are never
  inferred from the visualization.
- Lazy-loaded NGL WebGL viewer for curated MT-ATP6/MT-ND4 residue inspection,
  with bundled local-first RCSB experimental coordinates, mapped-chain/full-
  complex views, network fallbacks, B-factor semantics, 5 Å molecular
  context, checksum verification, and provenance links.
- Offline HTML report and React web UI consume the same JSON schema.
- Stable native/C/Rust error categories and a versioned HTTP failure envelope;
  callers never need to parse diagnostic prose. Route-level HTTP goldens cover
  authentication, malformed IDs/multipart, job states, JSON results, and HTML
  downloads through the real Axum middleware stack.

## Repository

```text
core/               C++20 analysis core
ffi/                C ABI wrapper used from Rust
cli/                Rust CLI and standalone HTML report generator
server/             Rust axum backend
visualization-lib/  TypeScript circular mtDNA plot library
web/                React/Vite frontend
fixtures/           Deterministic FASTQ/SAM truth and expected outputs
docs/               User, scientific, roadmap, and release evidence
scripts/            Bootstrap, fixture, resource, and verification tools
```

## Build And Test

For a user-facing setup guide, including native dependencies and Docker-first
installation, see [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md).

Audit the complete native toolchain without changing the host, or install it on
a supported Debian/Ubuntu/Arch host and run the release gate:

```bash
bash scripts/bootstrap.sh --check
bash scripts/bootstrap.sh --install --verify
```

The installer deliberately uses the host package manager for native tools and
the committed Cargo/npm lockfiles for project packages. It does not silently
install Docker unless `--include-docker` is supplied, and it cannot install
laboratory validation evidence, clinical-data licences, or deployment services.

Run the full local verification gate:

```bash
bash scripts/verify.sh
```

```bash
npm ci
cmake -S . -B build -DMITO_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build
cargo check --workspace
cargo test --workspace
npm run typecheck
npm --workspace visualization-lib run build
npm --workspace web run build
cargo run -p mito-cli -- validate-fixture \
  --input fixtures/truth_snp.sam \
  --expected-vcf fixtures/truth_snp.expected.vcf \
  --expected-passed 3 \
  --expected-numt 1 \
  --expected-snp 3243:A:G
```

If htslib headers are installed (`pkg-config --exists htslib`), CMake enables
BAM/CRAM support automatically. Without htslib, FASTQ/SAM still work and
BAM/CRAM inputs return a clear rebuild error.

Local htslib build used on Arch without a package:

```bash
curl -L https://github.com/samtools/htslib/releases/download/1.23.1/htslib-1.23.1.tar.bz2 \
  -o .deps/htslib-1.23.1.tar.bz2
tar -xjf .deps/htslib-1.23.1.tar.bz2 -C .deps
cd .deps/htslib-1.23.1
./configure --prefix="$PWD/../../.local/htslib" --disable-libcurl
make
make install
```

Build with local htslib and optional HDBSCAN-C++:

```bash
PKG_CONFIG_PATH="$PWD/.local/htslib/lib/pkgconfig" \
cmake -S . -B build-hft-deps -DMITO_BUILD_TESTS=ON \
  -DMITO_HDBSCAN_CPP_ROOT="$PWD/.deps/hdbscan-cpp-master/HDBSCAN-CPP"
cmake --build build-hft-deps -j2
LD_LIBRARY_PATH="$PWD/.local/htslib/lib" ctest --test-dir build-hft-deps --output-on-failure
```

## CLI

Check compiled native capabilities and external tools:

```bash
cargo run -p mito-cli -- doctor
```

The doctor output includes engine, result-schema, and error-schema versions.
Analysis failures are emitted as `[MITO-E....] message`; automation must branch
on the code described in [`docs/ERROR_CONTRACT.md`](docs/ERROR_CONTRACT.md).

Emit JSON:

```bash
cargo run -p mito-cli -- analyze -i fixtures/tiny.fastq --json
```

Use multiple C++ workers for per-read feature extraction:

```bash
cargo run -p mito-cli -- analyze -i sample.bam --json --threads 8
```

Tune the explicit evidence thresholds when a validated protocol requires it:

```bash
cargo run -p mito-cli -- analyze -i sample.bam --json --threads 8 \
  --min-mapq 30 --min-base-quality 20 --numt-threshold 0.30
```

Emit the opt-in RC2 evidence graph while keeping schema 0.5 as the default:

```bash
cargo run -p mito-cli -- analyze -i sample.bam --json --threads 8 \
  --evidence-graph --max-evidence-observations 5000000 --max-phase-links 1000000 \
  --evidence-page-size 4096 --evidence-pages-dir analysis-evidence-pages
```

When the protocol provides an explicit physical-molecule identifier, declare
its SAM tag rather than assuming QNAME independence:

```bash
cargo run -p mito-cli -- analyze -i sample.bam --json \
  --evidence-graph --molecule-id-tag MI --umi-tag RX --duplex-tag DX
```

Only the configured molecule-ID tag changes grouping. Missing or conflicting
required IDs remain auditable and are excluded from scientific evidence.

Schema 0.6 is a review draft. Missing sparse molecule/event pairs mean
`NOT_CALLABLE`; they must never be treated as reference observations. Phase
links marked `SUPPORT_CONDITIONED` do not have a whole-sample denominator.

Synthetic FASTQ header controls are test-only and disabled by default. The
`--allow-development-tags` flag exists solely for committed development
fixtures and must not be used for biological analysis.

Generate an offline report:

```bash
cargo run -p mito-cli -- analyze -i fixtures/tiny.fastq -o output.html
```

Export the schema-0.6 unified SNV/small-indel projection, indexed VCF, and
deterministic provenance:

```bash
cargo run -p mito-cli -- analyze -i sample.bam --evidence-graph --json \
  --vcf sample.vcf --tsv sample.variants.tsv \
  --bgzip-vcf sample.vcf.gz \
  --provenance-manifest sample.provenance.json
```

The VCF is a compatibility projection; JSON/paged evidence remains
authoritative. Schema-0.6 `FILTER=.` means thresholds are not scientifically
calibrated. Circular-origin indels that cannot be represented losslessly in a
linear VCF fail with an explicit error and remain available in JSON/TSV.

Compare the schema-0.6 callset with a pinned mtDNA-Server 2 v2.1.16
`variants.annotated.txt` output:

```bash
cargo run -p mito-cli -- compare-mt-dna-server2 \
  --result sample.json --comparator variants.annotated.txt \
  --sample sample-id --output comparator.differential.json
```

The resulting concordance and heteroplasmy deltas are reproducible comparator
evidence, not analytical truth, sensitivity/specificity, LoD, or clinical
validation.

Refresh the local clinical annotation cache:

```bash
cargo run -p mito-cli -- update-clinical
```

Refresh from the official ClinVar tab-delimited export:

```bash
cargo run -p mito-cli -- update-clinical --clinvar-live
```

Use a custom preprocessed MITOMAP/ClinVar TSV:

```bash
MITO_CLINICAL_ANNOTATIONS=/path/to/clinical_annotations.tsv \
  cargo run -p mito-cli -- analyze -i sample.sam --json
```

Clinical TSVs must follow assertion schema 1.0. The engine preserves each
source record, computes an explicit conflict state, and fails closed on malformed
or contradictory resource rows. See
[`docs/CLINICAL_ANNOTATION_SCHEMA_1.0.md`](docs/CLINICAL_ANNOTATION_SCHEMA_1.0.md).
An explicit override is labeled unpinned/unverified in the result provenance;
the bundled development subset remains checksum-pinned but is not complete.

Stage and activate an immutable governed snapshot:

```bash
cargo run -p mito-cli -- clinical-snapshot stage \
  --store /srv/mito/clinical-snapshots \
  --source /approved/clinical_annotations.tsv \
  --snapshot-id 2026-07-approved \
  --license-id LAB-CLINICAL-RESOURCE-2026 \
  --source-policy laboratory-approved-v1 \
  --activate
```

Verify freshness and obtain the exact path for the engine:

```bash
cargo run -p mito-cli -- clinical-snapshot status \
  --store /srv/mito/clinical-snapshots --max-age-days 30
export MITO_CLINICAL_ANNOTATIONS="$(cargo run --quiet -p mito-cli -- \
  clinical-snapshot status --store /srv/mito/clinical-snapshots \
  --max-age-days 30 | jq -r .active_path)"
```

Rollback swaps the active and previous snapshots only after the rollback target
passes schema, checksum, governance-label, and freshness validation:

```bash
cargo run -p mito-cli -- clinical-snapshot rollback \
  --store /srv/mito/clinical-snapshots --max-age-days 30
```

MITOMAP does not provide a stable public live API suitable for direct querying
from the app; MITOMAP rows are expected through the local TSV cache.

## Web Version

Start the backend:

```bash
cargo run -p mito-server
```

Start the frontend:

```bash
npm --workspace web run dev
```

Open:

```text
http://127.0.0.1:5173/
```

The demo result is available at:

```text
http://127.0.0.1:5173/result/demo
```

Cancel a backend job:

```bash
curl -X POST http://127.0.0.1:8080/cancel/<job_id>
```

Cancellation is cooperative: the server sets a per-job cancellation flag and the
native analysis core checks it while parsing reads, extracting features,
clustering, and calculating coverage. Third-party HDBSCAN execution remains a
single non-interruptible call.

The backend admits one analysis job at a time by default because each job can use
all available CPU threads. Tune this explicitly when deploying on a larger host:

```bash
MITO_MAX_ACTIVE_JOBS=16 MITO_MAX_CONCURRENT_JOBS=2 MITO_JOB_THREADS=8 \
MITO_MAX_UPLOAD_BYTES=67108864 MITO_MAX_EVIDENCE_OBSERVATIONS=5000000 \
MITO_MAX_PHASE_LINKS=1000000 cargo run -p mito-server
```

Uploads default to 64 MiB because the current result contract retains per-read
features and can expand substantially in memory. Raise the byte limit only with
a measured memory envelope. The active queued/running job cap defaults to 16.

Browser origins default to `http://127.0.0.1:5173` and
`http://localhost:5173`. Set a comma-separated deployment allowlist with
`MITO_CORS_ORIGINS`; the API no longer enables unrestricted CORS.
The provided Compose ports bind to loopback. The server refuses a non-loopback
bind unless `MITO_API_KEY` is set; protected endpoints then require
`Authorization: Bearer <key>`. A deployment still needs TLS, user-level
authorization, audit logging, and policy enforcement at a trusted ingress.

## Result Contract

Every layer consumes the same JSON document. The current development result
schema version is `0.5`.

- `metadata`: result, SV-edge, and complex-SV schema versions, sample, engine version,
  rCRS/custom reference, thread request, and versioned resource
  paths/checksums/sources/licenses.
- `filter_stats`: alignment-record and molecule counts, passed and NUMT-filtered
  molecules, evidence mode, and whether NUMT specificity is assessable.
- `genes`: rCRS gene intervals.
- `coverage` and `coverage_metrics`: rounded binned depth plus exact site-level
  mean, maximum, and >20x summary statistics for passed mapped reads.
- `coverage_metrics.mapping_quality_histogram`: MAPQ distribution from parsed
  passed, mapped alignment metadata.
- `svs`: canonical structural-event IDs, breakpoints, sorted evidence-source
  and strand-transition sets, segment count, known-event flag, and sorted
  supporting molecules. Equivalent CIGAR/forward-split/reverse-split evidence
  merges under SV event schema `1.0`.
- `complex_events`: exact strand-invariant ordered paths containing two or more
  split-alignment junction IDs and strand transitions, with segment count and
  sorted supporting molecules. These are observed paths, not inferred complete
  circular genomes.
- `clusters`: cluster ID, ranked PhyloTree assignment, alternatives, missing and
  extra sites, contamination warning, size, reads, SV signature, and complex-
  event signature.
- `reads`: per-read SNPs, SV and complex-event IDs, NUMT status, quality,
  cluster assignment, DBSCAN outlier flag, MAPQ, flags, reference name, and
  preserved aux tags.
- SNP `annotation`: local MITOMAP/ClinVar-style phenotype, pathogenicity,
  references, and source labels.
- SNP `structure`: optional protein/residue/structure mapping for 3D display.
- `variants`: authoritative locus-level SNP evidence with quality-filtered
  reference/alternate/other depth, callable DP, HF, Wilson 95% interval, and
  supporting molecule IDs.

## Analysis Modules

The React result page contains module tabs for:

- Haplogroups: cluster-level assignment field and signatures.
- Clinical: deterministic summary, explicit conflict state, and expandable
  source assertions with review/snapshot provenance.
- Coverage: mean depth, high-depth fraction, max depth, histogram preview.
- 3D Protein: NGL molecular viewer with representation modes, residue focus,
  experimental B-factor coloring, mapped-chain or full-complex scope, 5 Å
  residue context, spin, reset, provenance, PNG export, and adjacent read-
  cluster/molecule support. The viewer accepts curated PDB mappings only and
  never fabricates fallback coordinates.

Multi-sample, single-cell, and longitudinal tabs are intentionally not exposed
until their multi-input schemas, statistics, tests, and failure states are
implemented. The production UI does not present placeholder analyses.

## Production Readiness

Mito-Architect is designed as a single workbench for native mtDNA analysis,
interactive inspection, clinical annotation surfaces, and reproducible exports.
Production readiness is gated on the repository's own criteria: real aligned
input handling, transparent dependency setup, VCF/heteroplasmy interoperability,
NUMT validation, reproducible truth fixtures, documented failure modes,
cancellable jobs, and benchmarked performance envelopes.

## Benchmarks

Build the native benchmark target:

```bash
cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DMITO_BUILD_TESTS=ON -DMITO_BUILD_BENCHMARKS=ON
cmake --build build-bench -j2
MITO_BENCH_READS=10000 ./build-bench/core/mito_core_benchmark
bash scripts/verify_benchmark_budget.sh
bash scripts/verify_evidence_benchmark_budget.sh
```

The native benchmark reports phase-separated timings outside the deterministic
scientific JSON. The budget script runs three cold processes per supported
smoke size and fails when their median exceeds the documented same-host
development ceiling; it is not a substitute for the RC5 multi-host envelope.
The evidence-budget companion exercises schema 0.6 on an adverse high-mismatch
workload and bounds latency, JSON size, and native `getrusage` peak RSS so
evidence-graph growth cannot regress silently.

## Container And Verification

Build the combined image:

```bash
docker build -t mito-architect .
```

On Linux, make sure `docker.service` is running and the current login session
can access the Docker socket, usually through the `docker` group or `sudo`.

Run the backend:

```bash
docker compose up server
```

The canonical local verification gate builds/tests C++, Rust, TypeScript, truth
fixtures, the Release benchmark target, and, when requested, the Docker image:

```bash
bash scripts/verify.sh
```

Run the native sanitizer gate separately or include it in the canonical gate:

```bash
bash scripts/verify_sanitizers.sh
MITO_VERIFY_SANITIZERS=1 bash scripts/verify.sh
```

Include Docker in that local gate with:

```bash
MITO_VERIFY_DOCKER=1 bash scripts/verify.sh
```

After running the entire gate into a log on a clean candidate commit, generate
the exact release evidence index. The command refuses a dirty tree or a missing
verification log:

```bash
bash scripts/verify.sh 2>&1 | tee build/rc1-verify.log
bash scripts/generate_rc1_evidence.sh \
  --output build/rc1-evidence.json \
  --verification-log build/rc1-verify.log
```

Release notes for the current version are in
[`docs/RELEASE_NOTES_0.4.1.md`](docs/RELEASE_NOTES_0.4.1.md). The staged path to
1.0 is in [`docs/PLAN_1.0.md`](docs/PLAN_1.0.md), and the latest local evidence
is recorded in
[`docs/VERIFICATION_0.5_RC1_DEV.md`](docs/VERIFICATION_0.5_RC1_DEV.md).

## Threading Model

The C++ core uses a bounded worker pool for independent per-read feature
extraction: NUMT scoring, SNP calling, clinical cache lookup, CIGAR/SV parsing,
and read-level metadata capture run in parallel. Shared SV aggregation and
clustering are reduced deterministically after workers join. The DBSCAN fallback
compresses identical feature profiles and uses an inverted token index to avoid
comparing profiles that cannot be Jaccard neighbors. The CLI exposes workers
through `--threads`; the web backend uses `MITO_JOB_THREADS` and bounds concurrent
jobs with `MITO_MAX_CONCURRENT_JOBS`.

## Data Files

- `core/data/rcrs.fasta`: bundled RefSeq rCRS (`NC_012920.1`).
- `core/data/clinical_annotations.tsv`: checksum-pinned assertion-schema 1.0
  development subset used by tests, CLI, backend, and UI; not a complete
  clinical knowledge base.
- `core/data/phylotree-rcrs-17.3.xml`, `-weights.txt`, and `-rules.csv`:
  versioned haplogroup tree, weights, and official HaploGrep normalization
  rules with the upstream MIT license.
- `core/data/resource_manifest.tsv`: source, license, retrieval date, and
  SHA-256 identity for bundled scientific resources.
- `web/public/structures/`: curated checksum-verified RCSB BinaryCIF models used
  by the current variant-to-protein mappings.
- `fixtures/`: deterministic SAM/FASTQ truth inputs and expected VCF output;
  see [`fixtures/README.md`](fixtures/README.md) for BAM generation and a
  bounded public Oxford Nanopore validation dataset.

## Current Boundaries

The code intentionally does not fake unavailable production evidence or
resources. The HDBSCAN-C++ adapter is build-tested when explicitly configured;
complete clinical databases, durable
storage, service security, and blinded analytical validation remain outside
the current 0.5/RC1 software baseline and are mandatory later gates on the path
to 1.0. DBSCAN is implemented in-core as the current density-clustering baseline.

Algorithmic caveats in the current prototype:

- Plain FASTQ has no alignment coordinates, so production SNP evidence should
  come from aligned SAM/BAM/CRAM or an explicit variant-calling stage.
- NUMT specificity is assessable only when the input was competitively aligned
  to a nuclear-plus-mitochondrial reference. FASTQ and mtDNA-only inputs are
  explicitly marked non-assessable; their heuristic evidence is not sufficient
  for low-frequency interpretation.
- The pinned PhyloTree 17.3 token grammar, representative 28-profile/23-lineage
  HaploGrep top-three differential, three backmutation profiles, tree-order tie
  breaking, repeat-equivalent 3-prime CIGAR normalization, and the official
  callable-marker alignment-rule table are implemented. This is not a
  population accuracy study: broader HSD/IUPAC input grammar, mapper-level
  differential evidence, independent calibration, and review remain RC1/RC5
  work. Ranking quality is an evidence-similarity score, not a probability.
- The bundled clinical annotation table is a development subset, not a complete
  clinical knowledgebase.
- Cancellation is cooperative and is checked between major native analysis
  steps; blocking third-party file reads may only observe cancellation after the
  read call returns.
