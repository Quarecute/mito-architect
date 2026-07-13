#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NUCLEAR_FASTA="${MITO_NUCLEAR_FASTA:-}"
OUTPUT_DIR="${MITO_REFERENCE_DIR:-$ROOT_DIR/.data/reference}"

if [[ -z "$NUCLEAR_FASTA" || ! -f "$NUCLEAR_FASTA" ]]; then
  echo "Set MITO_NUCLEAR_FASTA to a versioned human nuclear FASTA." >&2
  exit 2
fi
if ! command -v minimap2 >/dev/null 2>&1; then
  echo "minimap2 is required to build the competitive reference index." >&2
  exit 2
fi

mkdir -p "$OUTPUT_DIR"
COMBINED_FASTA="$OUTPUT_DIR/nuclear-plus-rcrs.fasta"
TEMP_FASTA="$COMBINED_FASTA.part"

# Remove any pre-existing mitochondrial contig before appending the bundled
# NC_012920.1 sequence, preventing ambiguous duplicate reference names.
awk '
  /^>/ {
    name = substr($1, 2)
    keep = !(name == "chrM" || name == "MT" || name == "M" || name == "NC_012920.1")
  }
  keep { print }
' "$NUCLEAR_FASTA" > "$TEMP_FASTA"
awk '1' "$ROOT_DIR/core/data/rcrs.fasta" >> "$TEMP_FASTA"
mv "$TEMP_FASTA" "$COMBINED_FASTA"

minimap2 -d "$OUTPUT_DIR/nuclear-plus-rcrs.mmi" "$COMBINED_FASTA"
{
  minimap2 --version
  sha256sum "$NUCLEAR_FASTA" "$ROOT_DIR/core/data/rcrs.fasta" "$COMBINED_FASTA" \
    "$OUTPUT_DIR/nuclear-plus-rcrs.mmi"
} > "$OUTPUT_DIR/reference-manifest.txt"

echo "Competitive reference: $OUTPUT_DIR/nuclear-plus-rcrs.mmi"
echo "Provenance manifest: $OUTPUT_DIR/reference-manifest.txt"
