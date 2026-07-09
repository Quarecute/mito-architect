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

npm run typecheck
npm run build

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
  --expected-sv insertion:132+15 \
  --expected-sv soft_clip_right:142+13 \
  --expected-mapq 42 \
  --expected-mapq 60 \
  --expected-aux NM=1 \
  --expected-aux MD=2T7

cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DMITO_BUILD_TESTS=ON -DMITO_BUILD_BENCHMARKS=ON
cmake --build build-bench -j2
MITO_BENCH_READS="${MITO_BENCH_READS:-200}" ./build-bench/core/mito_core_benchmark

if [[ "${MITO_VERIFY_DOCKER:-0}" == "1" ]]; then
  docker build -t mito-architect .
fi
