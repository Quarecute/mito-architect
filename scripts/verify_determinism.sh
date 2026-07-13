#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT="${MITO_DETERMINISM_INPUT:-$ROOT_DIR/fixtures/truth_mixed.sam}"
THREADS="${MITO_DETERMINISM_THREADS:-1 2 4 8 16}"

for tool in jq sha256sum cargo; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "error: $tool is required for deterministic-output verification" >&2
    exit 2
  fi
done

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

baseline=""
for threads in $THREADS; do
  raw="$TMP_DIR/result-$threads.json"
  canonical="$TMP_DIR/result-$threads.canonical.json"
  cargo run --quiet --manifest-path "$ROOT_DIR/Cargo.toml" -p mito-cli --offline -- \
    analyze --input "$INPUT" --json --threads "$threads" > "$raw"
  jq --sort-keys \
    'del(.metadata.threads, .metadata.requested_threads)' \
    "$raw" > "$canonical"
  digest="$(sha256sum "$canonical" | awk '{print $1}')"
  echo "threads=$threads sha256=$digest"
  if [[ -z "$baseline" ]]; then
    baseline="$digest"
  elif [[ "$digest" != "$baseline" ]]; then
    echo "determinism failure: thread count $threads differs from baseline" >&2
    exit 1
  fi
done

echo "Deterministic result payload verified across threads: $THREADS"
