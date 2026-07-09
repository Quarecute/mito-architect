# Mito-Architect

Mito-Architect is a long-read mitochondrial DNA analysis workbench. It combines a
C++ analysis core, a Rust CLI/backend, a TypeScript visualization library, and a
React web UI around one JSON result contract.

The current build is an end-to-end prototype mtDNA analysis workbench: it
accepts FASTQ/SAM inputs directly, can use htslib for SAM/BAM/CRAM when htslib
is available at build time, uses a bundled rCRS reference (`NC_012920.1`),
detects tagged/CIGAR SVs including soft-clips, clusters reads by DBSCAN over
SNP/SV feature tokens, annotates known variants from a local MITOMAP/
ClinVar-style cache, and renders interactive reports.

This repository is not yet a validated clinical or production bioinformatics
pipeline. SNP discovery, NUMT filtering, haplogroup assignment, and several
analysis modules are intentionally lightweight extension points so the UI,
backend, FFI, and result contract can be exercised end to end.

FASTQ ingestion does not perform read alignment. Alignment-derived SNPs,
CIGAR SVs, MAPQ, and coverage require SAM/BAM/CRAM input; raw FASTQ currently
contributes read-level quality/NUMT heuristics and explicit test tags only.

## Highlights

- Bundled rCRS reference, length-checked at 16,569 bp.
- FASTA read input is rejected; FASTA is reference-only.
- FASTQ and SAM records are validated for required separators, field types,
  sequence/quality lengths, and CIGAR/query-length consistency.
- Optional htslib-backed SAM/BAM/CRAM reader with CIGAR, MAPQ, flags, and aux
  tag preservation.
- NUMT filtering hook with explicit read accounting.
- Large deletion/insertion calls with supporting reads.
- CIGAR-aware coverage excludes unmapped and NUMT-filtered reads and reports
  exact site-level mean, maximum, and >20x statistics alongside display bins.
- Cluster drill-down: click a cluster to inspect its SNPs, SVs, support,
  frequencies, molecule IDs, and genome distribution.
- Local clinical cache for MITOMAP/ClinVar-style annotations.
- Lazy-loaded NGL WebGL viewer for MT-ATP6/MT-ND4 protein-residue inspection.
- Offline HTML report and React web UI consume the same JSON schema.

## Repository

```text
core/               C++20 analysis core
ffi/                C ABI wrapper used from Rust
cli/                Rust CLI and standalone HTML report generator
server/             Rust axum backend
visualization-lib/  TypeScript circular mtDNA plot library
web/                React/Vite frontend
fixtures/           Small FASTQ smoke data
```

## Build And Test

For a user-facing setup guide, including native dependencies and Docker-first
installation, see [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md).

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

Emit JSON:

```bash
cargo run -p mito-cli -- analyze -i fixtures/tiny.fastq --json
```

Use multiple C++ workers for per-read feature extraction:

```bash
cargo run -p mito-cli -- analyze -i sample.bam --json --threads 8
```

Generate an offline report:

```bash
cargo run -p mito-cli -- analyze -i fixtures/tiny.fastq -o output.html
```

Export SNP support and heteroplasmy estimates as VCF:

```bash
cargo run -p mito-cli -- analyze -i sample.sam --json --vcf sample.vcf
```

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
MITO_MAX_UPLOAD_BYTES=67108864 cargo run -p mito-server
```

Uploads default to 64 MiB because the current result contract retains per-read
features and can expand substantially in memory. Raise the byte limit only with
a measured memory envelope. The active queued/running job cap defaults to 16.

Browser origins default to `http://127.0.0.1:5173` and
`http://localhost:5173`. Set a comma-separated deployment allowlist with
`MITO_CORS_ORIGINS`; the API no longer enables unrestricted CORS.
The provided Compose ports bind to loopback; exposing the API beyond the local
host requires an authenticated reverse proxy because the server itself does not
implement user authentication.

## Result Contract

Every layer consumes the same JSON document. The current result schema version
is `0.2`.

- `metadata`: schema version, sample, engine version, rCRS/custom reference,
  thread request.
- `filter_stats`: input reads, passed reads, NUMT-filtered reads.
- `genes`: rCRS gene intervals.
- `coverage` and `coverage_metrics`: rounded binned depth plus exact site-level
  mean, maximum, and >20x summary statistics for passed mapped reads.
- `coverage_metrics.mapping_quality_histogram`: MAPQ distribution from parsed
  passed, mapped alignment metadata.
- `svs`: structural variants, breakpoints, known-event flag, supporting reads.
- `clusters`: cluster ID, label, haplogroup field, size, reads, SV signature.
- `reads`: per-read SNPs, SV IDs, NUMT status, quality, cluster assignment,
  DBSCAN outlier flag, MAPQ, flags, reference name, and preserved aux tags.
- SNP `annotation`: local MITOMAP/ClinVar-style phenotype, pathogenicity,
  references, and source labels.
- SNP `structure`: optional protein/residue/structure mapping for 3D display.

## Analysis Modules

The React result page contains module tabs for:

- Haplogroups: cluster-level assignment field and signatures.
- Clinical: variant pathogenicity, phenotype, references, and sources.
- Coverage: mean depth, high-depth fraction, max depth, histogram preview.
- Multi-sample: sample x variant matrix surface for merged jobs.
- Single-cell: barcode-aware summary surface for `CB/CR` tagged reads.
- 3D Protein: NGL molecular viewer with representation modes, residue focus,
  spin, reset, and PNG export.
- Timeline: dated-sample surface for future tree and molecular-clock output.

Production-scale versions of these modules need the full local resource bundle:
PhyloTree build 17, full MITOMAP/ClinVar exports, optional MAFFT/FastTree/
IQ-TREE/TreeTime binaries, and curated protein structure mappings.

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
```

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

Include Docker in that local gate with:

```bash
MITO_VERIFY_DOCKER=1 bash scripts/verify.sh
```

Release notes for the current version are in
[`docs/RELEASE_NOTES_0.2.0.md`](docs/RELEASE_NOTES_0.2.0.md).

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
- `core/data/clinical_annotations.tsv`: small local annotation cache used by
  tests, fixtures, CLI, backend, and demo UI.

## Current Boundaries

The code intentionally does not fake unavailable production dependencies.
hg38 NUMT remapping, a native HDBSCAN library adapter, full PhyloTree
classification, and complete licensed clinical databases are extension points
behind the current APIs and JSON schema. DBSCAN is implemented in-core as the
current density-clustering baseline.

Algorithmic caveats in the current prototype:

- Plain FASTQ has no alignment coordinates, so production SNP evidence should
  come from aligned SAM/BAM/CRAM or an explicit variant-calling stage.
- NUMT filtering is a heuristic based on read name, read length, and GC
  fraction until nuclear remapping is implemented.
- Cancellation is cooperative and is checked between major native analysis
  steps; blocking third-party file reads may only observe cancellation after the
  read call returns.
