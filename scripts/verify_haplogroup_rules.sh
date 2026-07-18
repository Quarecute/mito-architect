#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STDERR_FILE="$(mktemp)"
trap 'rm -f "$STDERR_FILE"' EXIT

if MITO_PHYLOTREE_ALIGNMENT_RULES="$ROOT_DIR/fixtures/haplogroup/invalid_alignment_rules.csv" \
  cargo run --quiet --manifest-path "$ROOT_DIR/Cargo.toml" -p mito-cli --offline -- \
  analyze --input "$ROOT_DIR/fixtures/truth_snp.sam" --json \
  >/dev/null 2>"$STDERR_FILE"; then
  echo "error: truncated PhyloTree alignment rules unexpectedly succeeded" >&2
  exit 1
fi

if ! rg -q 'MITO-E1302' "$STDERR_FILE"; then
  echo "error: invalid alignment rules did not preserve MITO-E1302" >&2
  sed -n '1,20p' "$STDERR_FILE" >&2
  exit 1
fi

echo "Truncated PhyloTree alignment-rule resource failed closed with MITO-E1302"
