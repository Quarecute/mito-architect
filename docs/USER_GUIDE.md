# Mito-Architect User Guide

## Dependency Model

Mito-Architect has three dependency layers:

- Rust and JavaScript packages are resolved automatically by `cargo` and `npm ci`.
- The C++ core is built automatically by CMake when building the Rust FFI crate.
- Native system tools and libraries must be installed by the user, or avoided by
  running the Docker image.

For the least fragile setup, use Docker:

```bash
docker build -t mito-architect .
docker compose up server
```

The Docker client must be able to reach the daemon. On Linux this usually means
that `docker.service` is running and the current user is in the `docker` group,
or that commands are run through `sudo`:

```bash
systemctl start docker
sudo usermod -aG docker "$USER"
```

After changing group membership, start a new login session before running
`docker build`.

For a native setup on Debian/Ubuntu:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake libhts-dev nodejs npm pkg-config curl gzip
npm ci
cmake -S . -B build -DMITO_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
cargo test --workspace
npm run typecheck
npm run build
cargo run -p mito-cli -- validate-fixture \
  --input fixtures/truth_snp.sam \
  --expected-vcf fixtures/truth_snp.expected.vcf \
  --expected-passed 3 \
  --expected-numt 1 \
  --expected-snp 3243:A:G
cargo run -p mito-cli -- validate-fixture \
  --input fixtures/truth_mixed.sam \
  --expected-vcf fixtures/truth_mixed.expected.vcf \
  --expected-passed 2 \
  --expected-numt 0 \
  --expected-snp 3:T:G \
  --expected-sv deletion:110-121 \
  --expected-sv insertion:132+15 \
  --expected-sv soft_clip_right:142+13 \
  --expected-mapq 42 \
  --expected-mapq 60 \
  --expected-aux NM=1 \
  --expected-aux MD=2T7
```

After dependencies are installed, the same checks are available through:

```bash
bash scripts/verify.sh
```

Set `MITO_VERIFY_DOCKER=1` to include the Docker image build in that local gate.

On other Linux distributions, install equivalent packages:

- C++20 compiler and linker.
- CMake 3.24 or newer.
- Rust 1.80 or newer.
- Node.js and npm.
- `pkg-config`.
- htslib development headers for BAM/CRAM support.
- `curl` and `gzip` for `mito-cli update-clinical --clinvar-live`.

Without htslib, FASTQ and SAM still work. BAM and CRAM return a clear rebuild
error because those formats need htslib.

FASTQ and text SAM input is strict: malformed separators, invalid numeric SAM
fields, inconsistent sequence/quality lengths, and CIGAR/query-length mismatches
are rejected instead of being silently skipped.

Raw FASTQ is not aligned by the application. Use aligned SAM/BAM/CRAM for
alignment-derived SNP, structural-variant, MAPQ, and coverage results. FASTQ is
currently useful for read-quality/NUMT heuristics and tagged development data.

## Backend Resource Controls

One analysis job runs at a time by default, and that job may use all available
CPU threads. Set explicit limits for shared deployments:

```bash
MITO_MAX_ACTIVE_JOBS=16 \
MITO_MAX_CONCURRENT_JOBS=2 \
MITO_JOB_THREADS=8 \
MITO_MAX_UPLOAD_BYTES=67108864 \
MITO_CORS_ORIGINS=https://mito.example.org \
cargo run -p mito-server
```

`MITO_CORS_ORIGINS` is a comma-separated browser-origin allowlist. Its local
development default allows only `http://127.0.0.1:5173` and
`http://localhost:5173`.

The Compose file binds both ports to loopback. Do not publish the backend on a
public interface without an authenticated reverse proxy; CORS is a browser
policy, not API authentication.

The upload limit defaults to 64 MiB. Per-read JSON can be much larger than its
input, so increase `MITO_MAX_UPLOAD_BYTES` only after measuring peak resident
memory on representative data. The server accepts at most 16 queued/running jobs
by default; excess uploads receive HTTP 429.

## Production Validation Checklist

Before using results as scientific evidence, validate against a representative
truth set:

- Aligned SAM/BAM/CRAM with known SNPs, indels, soft clips, and MAPQ values.
- Expected VCF or comparable truth calls.
- Known NUMT-contaminated samples.
- Haplogroup-labelled samples.
- Large read sets that exercise clustering memory and runtime.

The repository includes smoke tests, type checks, and offline truth fixtures for
SNP, NUMT, SV, MAPQ, auxiliary-tag, and VCF behavior. Production datasets should
extend the same validation pattern with larger aligned samples.

## Exporting VCF

The CLI can write a VCF sidecar from the JSON result:

```bash
cargo run -p mito-cli -- analyze -i sample.sam --json --vcf sample.vcf
```

The current VCF export aggregates SNP support across passed reads and reports
heteroplasmy fraction as `support / passed_reads`. Treat that as an operational
export format, not a replacement for validation against a dedicated variant
caller and truth set.

This VCF denominator is not locus-specific depth. The JSON coverage metrics are
CIGAR-aware, but the current result contract does not yet expose per-locus depth
to the VCF exporter.

## Native Benchmark

Build and run the C++ benchmark target:

```bash
cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DMITO_BUILD_BENCHMARKS=ON
cmake --build build-bench -j2
MITO_BENCH_READS=10000 ./build-bench/core/mito_core_benchmark
```
