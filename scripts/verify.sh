#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

cmake -S . -B build -DMITO_BUILD_TESTS=ON
cmake --build build -j2
ctest --test-dir build --output-on-failure

cargo check --workspace --offline
cargo test --workspace --offline
cargo clippy --workspace --all-targets --offline -- -D warnings
cargo run -p mito-cli --offline -- doctor

npm run typecheck
npm run build
node web/src/lib/proteinStructureSources.test.mjs
bash scripts/verify_resources.sh

if [[ -f web/public/structures/manifest.sha256 ]]; then
  (
    cd web/public/structures
    sha256sum --check manifest.sha256
  )
else
  echo "Skipping optional offline protein-structure cache validation."
fi

cargo run -p mito-cli --offline -- validate-fixture \
  --input fixtures/truth_snp.sam \
  --expected-vcf fixtures/truth_snp.expected.vcf \
  --expected-passed 3 \
  --expected-numt 1 \
  --expected-snp 3243:A:G

cargo run -p mito-cli --offline -- validate-fixture \
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

cargo run -p mito-cli --offline -- validate-sv-fixture \
  --input fixtures/truth_split.sam \
  --expected-json fixtures/truth_split.expected.svs.json

bash scripts/verify_determinism.sh

if command -v samtools >/dev/null 2>&1 && pkg-config --exists htslib 2>/dev/null; then
  bash scripts/build_bam_fixture.sh
  cargo run -p mito-cli --offline -- validate-fixture \
    --input fixtures/generated/truth_mixed.bam \
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
  cargo run -p mito-cli --offline -- validate-fixture \
    --input fixtures/generated/truth_mixed.cram \
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
else
  echo "error: samtools and htslib are required for the complete verification gate" >&2
  exit 2
fi

if command -v bcftools >/dev/null 2>&1; then
  bcftools view --no-version -Ov fixtures/truth_snp.expected.vcf >/dev/null
  bcftools view --no-version -Ov fixtures/truth_mixed.expected.vcf >/dev/null
else
  echo "error: bcftools is required for VCF interoperability verification" >&2
  exit 2
fi

cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DMITO_BUILD_TESTS=ON -DMITO_BUILD_BENCHMARKS=ON
cmake --build build-bench -j2
MITO_BENCH_READS="${MITO_BENCH_READS:-200}" ./build-bench/core/mito_core_benchmark

if [[ "${MITO_VERIFY_SANITIZERS:-0}" == "1" ]]; then
  bash scripts/verify_sanitizers.sh
fi

HDBSCAN_ROOT="${MITO_HDBSCAN_CPP_ROOT:-$ROOT_DIR/.deps/hdbscan-cpp-master/HDBSCAN-CPP}"
if [[ -f "$HDBSCAN_ROOT/Hdbscan/hdbscan.hpp" ]]; then
  cmake -S . -B build-hdbscan -DCMAKE_BUILD_TYPE=Release -DMITO_BUILD_TESTS=ON \
    -DMITO_HDBSCAN_CPP_ROOT="$HDBSCAN_ROOT"
  cmake --build build-hdbscan -j2
  ctest --test-dir build-hdbscan --output-on-failure
else
  echo "Skipping HDBSCAN-C++ adapter validation (source tree not supplied)."
fi
MITO_BENCH_READS="${MITO_BENCH_READS:-1000}" ./build-bench/core/mito_core_benchmark
MITO_BENCH_READS="${MITO_BENCH_READS:-10000}" ./build-bench/core/mito_core_benchmark

if [[ "${MITO_VERIFY_DOCKER:-0}" == "1" ]]; then
  docker build -t mito-architect .
fi
