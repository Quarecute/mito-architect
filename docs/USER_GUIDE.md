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
sudo apt-get install -y build-essential cmake libhts-dev samtools minimap2 bcftools tabix \
  sra-toolkit nodejs npm rustc cargo pkg-config curl gzip jq git ripgrep
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
yay -S --needed htslib samtools minimap2 bcftools sra-tools
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
- bcftools, bgzip, and tabix for VCF interoperability, compression, and index checks.
- SRA Toolkit for bounded public-data acquisition.
- `curl`, `gzip`, `jq`, `git`, `ripgrep`, and GNU coreutils for resource
  refresh, deterministic verification, and immutable evidence indexing.

Dependency roles are intentionally distinct:

| Dependency | Required for | Bundled into application container |
| --- | --- | --- |
| C++20, CMake, pkg-config | native source build | build stage |
| Rust 1.80+, Cargo | CLI/server source build | build stage |
| Node 18+, npm | web source build | build stage |
| htslib | BAM/CRAM engine capability | yes |
| samtools | BAM/CRAM fixture construction and checks | verification image/host only |
| minimap2 | explicit pre-alignment today; first-class competitive pipeline in RC4 | no |
| bcftools | VCF interoperability gate | verification image/host only |
| bgzip + tabix | requested compressed/indexed VCF export and RC3 export gate | no, unless that export is requested |
| SRA Toolkit | optional bounded public validation fixture | no |
| Docker | alternate build/deployment path | host capability, opt-in |

“No” in the last column does not mean the tool is silently downloaded at
runtime. The current server consumes aligned SAM/BAM/CRAM; RC4 will package and
version the minimap2 competitive-mapping stage explicitly.

Without htslib, FASTQ and SAM still work. BAM and CRAM return a clear rebuild
error because those formats need htslib.

FASTQ and text SAM input is strict: malformed separators, invalid numeric SAM
fields, inconsistent sequence/quality lengths, and CIGAR/query-length mismatches
are rejected instead of being silently skipped.

Supplementary-alignment evidence is also strict. `SA` must be a `Z` tag whose
entries use `RNAME,POS,STRAND,CIGAR,MAPQ,NM;`. A malformed `SA` entry aborts the
analysis with `MITO-E1103` before NUMT filtering, SV calling, or parallel feature
extraction. Repair or regenerate the alignment; the engine deliberately does
not fall back to a partial primary-alignment interpretation.

Raw FASTQ is not aligned by the application. Use aligned SAM/BAM/CRAM for
alignment-derived SNP, structural-variant, MAPQ, and coverage results. FASTQ is
currently useful for parser/read-quality checks. Synthetic `snp=` and `sv=`
header controls are disabled by default and require the explicit
`--allow-development-tags` CLI option; never enable it for biological input.

The RC2 evidence graph is available as an explicit draft contract:

```bash
cargo run -p mito-cli -- analyze \
  --input sample.bam \
  --json \
  --evidence-graph \
  --max-evidence-observations 5000000 \
  --max-phase-links 1000000 \
  --evidence-page-size 4096 \
  --evidence-pages-dir analysis-evidence-pages
```

The default for non-negotiating callers remains schema 0.5. Schema 0.6 retains
source alignment fragments, records molecule-assembly decisions, emits
per-molecule callability, normalized SNV/small-indel/SV/complex events,
columnar observation pages, and callable-aware pairwise phase links. Missing molecule/event
pairs are `NOT_CALLABLE`; a support-conditioned SV link is not a whole-sample
denominator. The independent observation and phase-link limits fail with
`MITO-E1601` instead of returning unbounded or partial JSON. See
[`SCHEMA_0.6_DRAFT.md`](SCHEMA_0.6_DRAFT.md) before integrating it.

For protocols with an explicit physical-molecule SAM tag, configure it rather
than assuming QNAME independence:

```bash
cargo run -p mito-cli -- analyze --input sample.bam --json \
  --evidence-graph --molecule-id-tag MI --umi-tag RX --duplex-tag DX \
  --evidence-page-size 4096 --evidence-pages-dir analysis-evidence-pages
```

Only `--molecule-id-tag` changes grouping. UMI and duplex tags are preserved as
audited protocol metadata. Missing/conflicting mandatory IDs and conflicting
metadata remain visible but are excluded from scientific evidence. Do not use
`RX` alone as a physical-molecule key unless the wet-lab protocol explicitly
defines it that way.

The web client sends multipart `evidence_graph=true` and therefore requests
schema 0.6 explicitly; multipart option order is irrelevant. API clients may
also send `evidence_page_size`, `molecule_id_tag`, `umi_tag`, and `duplex_tag`.
Authenticated completed jobs expose:

- `GET /result/<job_id>/summary`: compact UI projection without embedded pages;
- `GET /result/<job_id>/evidence?molecule_id=<exact>&event_id=<exact>&state=<stored>&cursor=0&limit=100`:
  deterministic bounded search across every immutable observation page;
- `GET /result/<job_id>/evidence/<page_index>`: one immutable observation page;
- `GET /result/<job_id>`: the authoritative complete JSON download.

The browser can inspect a physical page or use exact molecule/event/state
filters against the global server index. Search results include stable global,
page, and row coordinates and are limited to 500 rows per request. Empty search
results and a molecule absent from a page are never interpreted as REF or event
absence. Other API
clients receive schema 0.5 unless they negotiate the evidence graph. Operators
can bound jobs with `MITO_MAX_EVIDENCE_OBSERVATIONS`,
`MITO_EVIDENCE_PAGE_SIZE`, and `MITO_MAX_PHASE_LINKS`.

## mtDNA-Server 2 differential

Create an authoritative schema-0.6 result, then compare its SNV/small-indel
callset with an upstream `variants.annotated.txt` file:

```bash
cargo run -p mito-cli -- analyze --input sample.bam \
  --evidence-graph --json > sample.mito-architect.json

cargo run -p mito-cli -- compare-mt-dna-server2 \
  --result sample.mito-architect.json \
  --comparator variants.annotated.txt \
  --sample sample-id \
  --comparator-version 2.1.16 \
  --filter-policy pass-only \
  --output sample.mtdna-server-2.differential.json
```

The adapter is intentionally pinned to mtDNA-Server 2 v2.1.16 and rCRS
`NC_012920.1`. It rejects unrecognized versions, malformed headers or rows,
invalid HF/coverage accounting, duplicate comparison alleles and ambiguous
multi-sample selection. The report records both input SHA-256 values, matched
and tool-only alleles, callset concordance and HF deltas. Optional
`--min-call-concordance` and `--max-mean-hf-delta` development gates write the
report before failing.

Comparator concordance is not a truth benchmark. It must not be reported as
sensitivity, specificity, LoD, clinical validity, or superiority. Those claims
require independently frozen truth and protocols in RC5.

## Validation Data

The committed SAM fixtures provide deterministic SNP/SV/MAPQ/aux-tag truth.
`truth_complex.sam` additionally verifies that forward and reverse observations
of the same three-segment molecule produce one ordered complex-event path.
`truth_snp_edges.sam` and `truth_molecule_edges.sam` pin circular,
multi-allelic, strand, quality, flag-exclusion, and molecule-conflict evidence.
The negative manifest validates strict parser failures:

```bash
cargo run -p mito-cli -- validate-evidence-fixture \
  --input fixtures/truth_molecule_edges.sam \
  --expected-json fixtures/truth_molecule_edges.expected.evidence.json
cargo run -p mito-cli -- validate-error-manifest \
  --manifest fixtures/negative/error_manifest.json
cargo run -p mito-cli -- validate-haplogroup-manifest \
  --manifest fixtures/haplogroup/haplogrep3-3.3.2.expected.json
```

The haplogroup manifest pins the HaploGrep executable, PhyloTree 17.3 tree,
weights, alignment rules, input HSD, tested ranges, and expected results by
SHA-256. Twenty-eight profiles across 23 declared lineage groups require the
same ordered top three as HaploGrep 3.3.2; H2c1, H1b1g, and V7b explicitly
verify three inherited-marker removals. One aligned-SAM case contains two
repeat-equivalent CIGAR representations and requires both to normalize to the
same `523d`, `524d`, and `573.1CCCC` markers and callable coordinates. A second
case requires motif rotation from `8270.1ACCCCCTCT` to
`8289.1CCCCCTCTA`. Three development FASTQ controls pin compound official
alignment rules around 309-310, 188-196, and 56-66. This is representative
differential evidence, not a population accuracy validation.

For aligned input, `haplogroup_markers` records the substitution and small-indel
evidence derived from the read. `haplogroup_callable_ranges` records where that
evidence was callable. Expected tree markers outside those intervals are not
counted as missing. The reported haplogroup score is a weighted similarity used
for ranking, not a calibrated probability or clinical confidence.

## Clinical Assertions and Variant QC

Clinical annotation schema 1.0 keeps every source record separate. In the
Clinical panel, expand a variant to inspect original and normalized
significance, disease, source record and allele IDs, review status, assertion
date, source URL, references, snapshot version, and retrieval date. A red
`conflicting` state means incompatible classified groups are present; never
replace this source review with the summary label alone.

The bundled rows are a small development subset. `update-clinical
--clinvar-live` writes a validated atomic cache from ClinVar's aggregate
`variant_summary` export, but it is not equivalent to all submitter-level SCV
records. The CLI now implements validate-before-publish staging, checksum
manifests, atomic activation, freshness checks, corruption detection, and
verified rollback; source licensing, content approval, and snapshot ownership
remain laboratory responsibilities. The exact resource/result contract is
[`CLINICAL_ANNOTATION_SCHEMA_1.0.md`](CLINICAL_ANNOTATION_SCHEMA_1.0.md).

```bash
cargo run -p mito-cli -- clinical-snapshot stage \
  --store /srv/mito/clinical-snapshots \
  --source /approved/clinical_annotations.tsv \
  --snapshot-id 2026-07-approved \
  --license-id LAB-CLINICAL-RESOURCE-2026 \
  --source-policy laboratory-approved-v1 --activate

cargo run -p mito-cli -- clinical-snapshot status \
  --store /srv/mito/clinical-snapshots --max-age-days 30

cargo run -p mito-cli -- clinical-snapshot rollback \
  --store /srv/mito/clinical-snapshots --max-age-days 30
```

Each aggregate SNP also reports:

- alternate/reference/other/callable molecule support;
- forward/reverse support and an absolute strand-fraction delta;
- count, mean, minimum, and maximum Phred quality per allele group;
- normalized read center proximity and alternate/reference bias delta.

These values help inspect rare calls and molecule distribution. They are not
validated filters and do not establish LoB, LoD, LoQ, or pathogenicity.

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

## Interpreting Complex Rearrangements

The Rearrangements module lists a multi-junction molecule path when one read has
two or more split-alignment junctions. Select the path to inspect its component
canonical SV IDs, traversal order, strand transitions, supporting molecules,
and retained primary/supplementary alignment fragments including `SA`
provenance. The forward path and
exact reverse-complement traversal share one identifier, so read strand does
not split support.

Treat this as observed molecule evidence, not a completed genome assembly. The
engine does not rotate a partially observed path around the circular reference,
infer missing junctions, or claim that the read closes a circle. Confirm complex
architectures with independently characterized truth and inspect the source
alignments before biological interpretation.

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
MITO_MAX_EVIDENCE_OBSERVATIONS=5000000 \
MITO_MAX_PHASE_LINKS=1000000 \
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

## Exporting unified variants and provenance

Schema 0.6 derives SNV and small-indel JSON, TSV, and VCF records from the same
molecule/event evidence. Write all interoperable sidecars with:

```bash
cargo run -p mito-cli -- analyze -i sample.bam --evidence-graph --json \
  --vcf sample.vcf \
  --tsv sample.variants.tsv \
  --bgzip-vcf sample.vcf.gz \
  --provenance-manifest sample.provenance.json
```

`--bgzip-vcf` runs bgzip and tabix and publishes both `.vcf.gz` and `.tbi`
atomically as a pair. It requires the native `bgzip` and `tabix` executables.
`--tsv`, `--bgzip-vcf`, and evidence-page export automatically negotiate schema
0.6. Plain `--vcf` without `--evidence-graph` retains schema-0.5 SNP
compatibility.

In schema 0.6, VCF `DP` is event-callable molecule depth, `AD` is reference and
alternate molecule support, `AC`/`MOLECULE_SUPPORT` is alternate support,
`HF=AC/DP`, and `HF_CI95` is a Wilson interval. `ODC`, `LOWQ`, `CONFLICT`,
`STRAND_SUPPORT`, `MQ`, `BQ`, `NUMT_ASSESSABLE`, and `QC_FLAGS` preserve
observable QC facts. `FILTER=.` means calibrated hard filtering was not applied;
it must not be read as `PASS`.

The provenance manifest records software/schema versions, Git commit when
available, input/reference/resource identities and SHA-256 values, calling
parameters, deterministic policy, command line, authoritative JSON digest, and
all export digests. Timestamps are intentionally omitted so the same invocation
can reproduce the same manifest.

A small indel whose allele representation crosses the 16,569/1 boundary cannot
be represented losslessly as one linear VCF allele. The VCF export therefore
fails explicitly and directs the user to authoritative JSON/TSV; it never drops
the event silently.

Default
SNP thresholds are MAPQ 20, base quality 10, and excluded SAM flags `0xF00`
(secondary, QC-fail, duplicate, supplementary). Override them explicitly:

```bash
cargo run -p mito-cli -- analyze -i sample.bam --json \
  --min-mapq 30 --min-base-quality 20 --excluded-snp-flags 3840
```

For the backend use `MITO_MIN_MAPQ`, `MITO_MIN_BASE_QUALITY`, and
`MITO_EXCLUDED_SNP_FLAGS` (decimal or hexadecimal for the flag mask). Analysis
metadata records the effective values.

Clinical VCF summaries use `CLNSIG`, `CLNSRC`, and `CLNCONFLICT`. Use JSON for
the complete assertion records and snapshot provenance.

## Native Benchmark

Build and run the C++ benchmark target:

```bash
cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DMITO_BUILD_BENCHMARKS=ON
cmake --build build-bench -j2
MITO_BENCH_READS=10000 ./build-bench/core/mito_core_benchmark
```

The benchmark reports total time, result bytes, and monotonic microsecond
timings for reference loading, ingest, resource loading, feature extraction,
event merging, clustering, evidence aggregation, haplogroup assignment, and
serialization. Operational timings are returned separately by the C++
`analyze_profiled()` API and never enter the deterministic scientific JSON.

Run the same-host three-process median regression budget with:

```bash
bash scripts/verify_benchmark_budget.sh
```

The development limits are 1,000 ms for 200 reads, 2,800 ms for 1,000 reads,
and 21,000 ms for 10,000 reads. Override them with
`MITO_BENCH_MAX_MS_200`, `MITO_BENCH_MAX_MS_1000`, and
`MITO_BENCH_MAX_MS_10000` when establishing a documented host-specific
baseline. These smoke budgets are not the RC5 laboratory workload envelope.
