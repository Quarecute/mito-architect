#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ACCESSION="${MITO_PUBLIC_ACCESSION:-SRR18110025}"
MAX_SPOTS="${MITO_PUBLIC_MAX_SPOTS:-5000}"
THREADS="${MITO_PUBLIC_THREADS:-4}"
OUT_DIR="${MITO_PUBLIC_OUT_DIR:-$ROOT_DIR/.data/public/$ACCESSION}"
ALIGNMENT_REFERENCE="${MITO_ALIGNMENT_REFERENCE:-$ROOT_DIR/core/data/rcrs.fasta}"

for tool in prefetch fasterq-dump minimap2 samtools; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "error: $tool is required (install sra-tools, minimap2, and samtools)" >&2
    exit 2
  fi
done

case "$MAX_SPOTS" in
  ''|*[!0-9]*) echo "error: MITO_PUBLIC_MAX_SPOTS must be a positive integer" >&2; exit 2 ;;
esac
if [[ "$MAX_SPOTS" -lt 1 ]]; then
  echo "error: MITO_PUBLIC_MAX_SPOTS must be at least 1" >&2
  exit 2
fi

mkdir -p "$OUT_DIR"
if [[ ! -s "$ALIGNMENT_REFERENCE" ]]; then
  echo "error: alignment reference does not exist or is empty: $ALIGNMENT_REFERENCE" >&2
  exit 2
fi

echo "Fetching $ACCESSION (at most $MAX_SPOTS spots)."
echo "This is public human sequence data; follow your institution's data-handling policy."
prefetch --max-size 20G "$ACCESSION"
fasterq-dump --threads "$THREADS" --progress -X "$MAX_SPOTS" -O "$OUT_DIR" "$ACCESSION"

FASTQ="$OUT_DIR/$ACCESSION.fastq"
if [[ ! -s "$FASTQ" ]]; then
  echo "error: fasterq-dump did not create $FASTQ" >&2
  exit 3
fi

minimap2 -t "$THREADS" -ax map-ont "$ALIGNMENT_REFERENCE" "$FASTQ" \
  | samtools sort -@ "$THREADS" -o "$OUT_DIR/$ACCESSION.rcrs.bam" -
samtools index -@ "$THREADS" "$OUT_DIR/$ACCESSION.rcrs.bam"
samtools view -h -o "$OUT_DIR/$ACCESSION.rcrs.sam" "$OUT_DIR/$ACCESSION.rcrs.bam"
samtools quickcheck -v "$OUT_DIR/$ACCESSION.rcrs.bam"

{
  echo "accession=$ACCESSION"
  echo "max_spots=$MAX_SPOTS"
  echo "alignment_reference=$ALIGNMENT_REFERENCE"
  echo "reference_sha256=$(sha256sum "$ALIGNMENT_REFERENCE" | awk '{print $1}')"
  echo "fastq_sha256=$(sha256sum "$FASTQ" | awk '{print $1}')"
  echo "bam_sha256=$(sha256sum "$OUT_DIR/$ACCESSION.rcrs.bam" | awk '{print $1}')"
  minimap2 --version 2>&1 | sed 's/^/minimap2=/'
  samtools --version | sed -n '1s/^/samtools=/'
} > "$OUT_DIR/MANIFEST.txt"

echo "Public validation bundle is ready in $OUT_DIR"
echo "Analyze with: cargo run -p mito-cli -- analyze -i $OUT_DIR/$ACCESSION.rcrs.bam --json"
