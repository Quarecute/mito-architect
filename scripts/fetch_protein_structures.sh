#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT_DIR="${MITO_STRUCTURE_DIR:-$ROOT_DIR/web/public/structures}"
mkdir -p "$OUTPUT_DIR"

fetch_model() {
  local pdb_id="$1"
  local lower_id="${pdb_id,,}"
  local destination="$OUTPUT_DIR/$lower_id.bcif"
  local temporary="$destination.part"

  curl --fail --location --retry 3 --connect-timeout 10 --max-time 180 \
    --output "$temporary" "https://models.rcsb.org/$lower_id.bcif"
  test -s "$temporary"
  mv "$temporary" "$destination"
}

fetch_model 8H9S
fetch_model 9I4I

(
  cd "$OUTPUT_DIR"
  sha256sum 8h9s.bcif 9i4i.bcif > manifest.sha256
  sha256sum --check manifest.sha256
)

echo "Cached curated protein structures in $OUTPUT_DIR"
