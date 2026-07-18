#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REPETITIONS="${MITO_BENCH_REPETITIONS:-3}"
if [[ ! "$REPETITIONS" =~ ^[0-9]+$ ]] || ((REPETITIONS < 3 || REPETITIONS % 2 == 0)); then
  echo "error: MITO_BENCH_REPETITIONS must be an odd integer >= 3" >&2
  exit 2
fi

cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DMITO_BUILD_TESTS=ON \
  -DMITO_BUILD_BENCHMARKS=ON
cmake --build build-bench -j"${MITO_BUILD_JOBS:-2}"

benchmark_budget_ms() {
  case "$1" in
  200) echo "${MITO_BENCH_MAX_MS_200:-1000}" ;;
  1000) echo "${MITO_BENCH_MAX_MS_1000:-2800}" ;;
  10000) echo "${MITO_BENCH_MAX_MS_10000:-21000}" ;;
  *)
    echo "error: no default budget for $1 reads; set a supported count" >&2
    return 2
    ;;
  esac
}

for reads in 200 1000 10000; do
  budget_ms="$(benchmark_budget_ms "$reads")"
  if [[ ! "$budget_ms" =~ ^[0-9]+$ ]] || ((budget_ms == 0)); then
    echo "error: benchmark budget for $reads reads must be a positive integer" >&2
    exit 2
  fi

  elapsed_values=()
  for ((run = 1; run <= REPETITIONS; ++run)); do
    output="$(MITO_BENCH_READS="$reads" ./build-bench/core/mito_core_benchmark)"
    printf '%s\n' "$output"
    elapsed_ms="$(awk -F= '$1 == "elapsed_ms" { print $2 }' <<<"$output")"
    if [[ ! "$elapsed_ms" =~ ^[0-9]+$ ]]; then
      echo "error: benchmark did not emit one integer elapsed_ms value" >&2
      exit 2
    fi
    elapsed_values+=("$elapsed_ms")
  done

  mapfile -t sorted_values < <(printf '%s\n' "${elapsed_values[@]}" | sort -n)
  median_ms="${sorted_values[REPETITIONS / 2]}"
  if ((median_ms > budget_ms)); then
    echo "error: $reads-read median ${median_ms}ms exceeds ${budget_ms}ms budget" >&2
    exit 1
  fi
  echo "benchmark budget passed: reads=$reads median_ms=$median_ms budget_ms=$budget_ms"
done
