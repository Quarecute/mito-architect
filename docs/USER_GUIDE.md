# Mito-Architect User Guide

## Dependency Model

Mito-Architect has four dependency layers:

- Rust and JavaScript packages are resolved automatically by `cargo` and `npm ci`.
- The C++ core is built automatically by CMake when building the Rust FFI crate.
- Native system tools and libraries are installed through the host package
  manager by the repository bootstrap. Application build dependencies are
  already contained in the Docker image when using the container path.
- Large reference bundles, governed clinical snapshots, and independently
  characterized validation material are deliberate operator inputs. They are
  not ordinary software dependencies and are never downloaded implicitly.

Audit without changing the machine:

```bash
bash scripts/bootstrap.sh --check
```

On Debian/Ubuntu or Arch Linux, install the native development/tooling set,
resolve the committed Rust/JavaScript lockfiles, rebuild the htslib bridge, and
run the canonical gate:

```bash
bash scripts/bootstrap.sh --install --verify
```

The script requests `sudo` only for system packages. Omit `--verify` for a
faster dependency-only setup, add `--yes` for a non-interactive package-manager
run, and add `--include-docker` only when a container workflow is required.
Run `bash scripts/bootstrap.sh --help` for the exact contract. On unsupported
Linux distributions, `--check` remains useful and reports every missing item.

For an isolated application build, use Docker:

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
sudo apt-get install -y build-essential cmake libhts-dev samtools minimap2 bcftools \
  sra-toolkit nodejs npm rustc cargo pkg-config curl gzip jq
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
  --expected-sv insertion:131+15 \
  --expected-sv soft_clip_right:141+13 \
  --expected-mapq 42 \
  --expected-mapq 60 \
  --expected-aux NM=1 \
  --expected-aux MD=2T7
cargo run -p mito-cli -- validate-sv-fixture \
  --input fixtures/truth_split.sam \
  --expected-json fixtures/truth_split.expected.svs.json
```

For Arch Linux, the current bioinformatics packages are available through AUR:

```bash
yay -S --needed samtools minimap2 bcftools sra-tools
cargo clean -p mito-ffi  # required only if Cargo previously cached a no-htslib build
cargo run -p mito-cli -- doctor
```

`mito-cli doctor` must report `BAM/CRAM reader: enabled` before BAM or CRAM is
accepted as a verified capability. It also reports result-contract and error-
contract versions so recorded runs can be interpreted against the correct
schemas.

Analysis errors use stable codes such as `MITO-E1103` for malformed input and
`MITO-E1501` for cancellation. The backend returns a JSON error object with
`schema_version`, `code`, `message`, and `retryable`; internal host details are
not exposed to clients. Automation must branch on `code`, never message text.
The complete catalog is in [`ERROR_CONTRACT.md`](ERROR_CONTRACT.md).

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
- samtools and minimap2 for fixture generation and competitive mapping.
- bcftools for VCF interoperability checks.
- SRA Toolkit for bounded public-data acquisition.
- `curl`, `gzip`, and `jq` for resource refresh and deterministic verification.

Dependency roles are intentionally distinct:

| Dependency | Required for | Bundled into application container |
| --- | --- | --- |
| C++20, CMake, pkg-config | native source build | build stage |
| Rust 1.80+, Cargo | CLI/server source build | build stage |
| Node 18+, npm | web source build | build stage |
| htslib | BAM/CRAM engine capability | yes |
| samtools | BAM/CRAM fixture construction and checks | verification image/host only |
| minimap2 | explicit pre-alignment today; first-class pipeline in RC3 | no |
| bcftools | VCF interoperability gate | verification image/host only |
| SRA Toolkit | optional bounded public validation fixture | no |
| Docker | alternate build/deployment path | host capability, opt-in |

“No” in the last column does not mean the tool is silently downloaded at
runtime. The current server consumes aligned SAM/BAM/CRAM; RC3 will package and
version the minimap2 competitive-mapping stage explicitly.

Without htslib, FASTQ and SAM still work. BAM and CRAM return a clear rebuild
error because those formats need htslib.

FASTQ and text SAM input is strict: malformed separators, invalid numeric SAM
fields, inconsistent sequence/quality lengths, and CIGAR/query-length mismatches
are rejected instead of being silently skipped.

Raw FASTQ is not aligned by the application. Use aligned SAM/BAM/CRAM for
alignment-derived SNP, structural-variant, MAPQ, and coverage results. FASTQ is
currently useful for parser/read-quality checks. Synthetic `snp=` and `sv=`
header controls are disabled by default and require the explicit
`--allow-development-tags` CLI option; never enable it for biological input.

## Validation Data

The committed SAM fixtures provide deterministic SNP/SV/MAPQ/aux-tag truth.
Build the same mixed fixture as coordinate-sorted BAM and reference-backed CRAM:

```bash
bash scripts/build_bam_fixture.sh
bash scripts/verify_determinism.sh
```

For realistic Oxford Nanopore read lengths and alignment behavior, download a
bounded subset of public run `SRR18110025` (BioProject `PRJNA809571`) and align
it to the bundled rCRS reference:

```bash
bash scripts/fetch_public_fixture.sh
```

The latter requires SRA Toolkit, minimap2, and samtools, writes only under
`.data/`, and records checksums and tool versions. It is a realistic public
sample, not a calibrated low-frequency truth material. See
[`fixtures/README.md`](../fixtures/README.md) for scope and controls.

## Protein Structure Viewer

The NGL panel uses curated experimental structures for the currently annotated
human MT-ATP6 and MT-ND4 variants: RCSB PDB 8H9S chain N and PDB 9I4I chain r.
It first loads the bundled cache under `web/public/structures`, then tries
RCSB's BinaryCIF and mmCIF services. To verify and deliberately refresh the
bundled models:

```bash
bash scripts/fetch_protein_structures.sh
```

Model loading happens once per selected variant; representation, color,
context, chain/full-complex scope, and spin changes reuse the loaded
coordinates. The UI labels experimental B-factor coloring explicitly. The
current viewer accepts only curated four-character PDB identifiers and does not
query a prediction database.

The orange residue is the curated protein-coordinate mapping and the optional
licorice neighborhood shows atoms within 5 Å. Follow the source link to inspect
the model record. If no curated structure exists or every local/network source
fails, the viewer reports each failed source and does not synthesize
coordinates.

## Backend Resource Controls

One analysis job runs at a time by default, and that job may use all available
CPU threads. Set explicit limits for shared deployments:

```bash
MITO_MAX_ACTIVE_JOBS=16 \
MITO_MAX_CONCURRENT_JOBS=2 \
MITO_JOB_THREADS=8 \
MITO_MAX_UPLOAD_BYTES=67108864 \
MITO_CORS_ORIGINS=https://mito.example.org \
MITO_API_KEY='replace-with-a-secret-from-your-secret-manager' \
cargo run -p mito-server
```

`MITO_CORS_ORIGINS` is a comma-separated browser-origin allowlist. Its local
development default allows only `http://127.0.0.1:5173` and
`http://localhost:5173`.

The Compose file binds both ports to loopback. For a non-loopback
`MITO_SERVER_ADDR`, the backend refuses to start unless `MITO_API_KEY` is set.
Analysis endpoints then require `Authorization: Bearer <key>`; `/healthz` and
`/readyz` remain available to orchestrator probes. Do not embed a long-lived API
key in a public browser bundle: terminate TLS and user authentication at a
trusted reverse proxy and inject or validate service credentials there. CORS is
a browser policy, not authorization.

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

VCF `DP` is the quality-filtered locus-callable A/C/G/T molecule depth, `AC` is
alternate molecule depth, `HF=AC/DP`, and `HF_CI95` is a Wilson interval. Default
SNP thresholds are MAPQ 20, base quality 10, and excluded SAM flags `0xF00`
(secondary, QC-fail, duplicate, supplementary). Override them explicitly:

```bash
cargo run -p mito-cli -- analyze -i sample.bam --json \
  --min-mapq 30 --min-base-quality 20 --excluded-snp-flags 3840
```

For the backend use `MITO_MIN_MAPQ`, `MITO_MIN_BASE_QUALITY`, and
`MITO_EXCLUDED_SNP_FLAGS` (decimal or hexadecimal for the flag mask). Analysis
metadata records the effective values.

## Native Benchmark

Build and run the C++ benchmark target:

```bash
cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DMITO_BUILD_BENCHMARKS=ON
cmake --build build-bench -j2
MITO_BENCH_READS=10000 ./build-bench/core/mito_core_benchmark
```
