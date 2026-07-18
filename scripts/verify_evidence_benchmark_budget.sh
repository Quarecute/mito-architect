#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REPETITIONS="${MITO_EVIDENCE_BENCH_REPETITIONS:-3}"
READS="${MITO_EVIDENCE_BENCH_READS:-200}"
BUDGET_MS="${MITO_EVIDENCE_BENCH_MAX_MS:-8000}"
MAX_JSON_BYTES="${MITO_EVIDENCE_BENCH_MAX_JSON_BYTES:-60000000}"
MAX_RSS_KIB="${MITO_EVIDENCE_BENCH_MAX_RSS_KIB:-800000}"

for value_name in REPETITIONS READS BUDGET_MS MAX_JSON_BYTES MAX_RSS_KIB; do
  value="${!value_name}"
  if [[ ! "$value" =~ ^[0-9]+$ ]] || ((value == 0)); then
    echo "error: $value_name must be a positive integer" >&2
    exit 2
  fi
done
if ((REPETITIONS < 3 || REPETITIONS % 2 == 0)); then
  echo "error: MITO_EVIDENCE_BENCH_REPETITIONS must be an odd integer >= 3" >&2
  exit 2
fi
cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DMITO_BUILD_TESTS=ON \
  -DMITO_BUILD_BENCHMARKS=ON
cmake --build build-bench -j"${MITO_BUILD_JOBS:-2}"

elapsed_values=()
json_values=()
rss_values=()
for ((run = 1; run <= REPETITIONS; ++run)); do
  output="$(MITO_BENCH_SCHEMA=0.6 MITO_BENCH_READS="$READS" \
    ./build-bench/core/mito_core_benchmark)"
  printf '%s\n' "$output"
  schema="$(awk -F= '$1 == "schema" { print $2 }' <<<"$output")"
  elapsed_ms="$(awk -F= '$1 == "elapsed_ms" { print $2 }' <<<"$output")"
  json_bytes="$(awk -F= '$1 == "json_bytes" { print $2 }' <<<"$output")"
  rss_kib="$(awk -F= '$1 == "peak_rss_kib" { print $2 }' <<<"$output")"
  if [[ "$schema" != "0.6" || ! "$elapsed_ms" =~ ^[0-9]+$ ||
        ! "$json_bytes" =~ ^[0-9]+$ || ! "$rss_kib" =~ ^[0-9]+$ ]]; then
    echo "error: evidence benchmark emitted an invalid result envelope" >&2
    exit 2
  fi
  if ((rss_kib > MAX_RSS_KIB)); then
    echo "error: evidence peak RSS ${rss_kib} KiB exceeds ${MAX_RSS_KIB} KiB budget" >&2
    exit 1
  fi
  if ((json_bytes > MAX_JSON_BYTES)); then
    echo "error: evidence JSON ${json_bytes} bytes exceeds ${MAX_JSON_BYTES} byte budget" >&2
    exit 1
  fi
  elapsed_values+=("$elapsed_ms")
  json_values+=("$json_bytes")
  rss_values+=("$rss_kib")
done

mapfile -t sorted_elapsed < <(printf '%s\n' "${elapsed_values[@]}" | sort -n)
mapfile -t sorted_json < <(printf '%s\n' "${json_values[@]}" | sort -n)
median_ms="${sorted_elapsed[REPETITIONS / 2]}"
median_json="${sorted_json[REPETITIONS / 2]}"
max_rss="$(printf '%s\n' "${rss_values[@]}" | sort -n | tail -n 1)"
if ((median_ms > BUDGET_MS)); then
  echo "error: schema 0.6 median ${median_ms}ms exceeds ${BUDGET_MS}ms budget" >&2
  exit 1
fi
echo "evidence benchmark budget passed: reads=$READS median_ms=$median_ms budget_ms=$BUDGET_MS median_json_bytes=$median_json max_json_bytes=$MAX_JSON_BYTES max_rss_kib=$max_rss rss_budget_kib=$MAX_RSS_KIB"
