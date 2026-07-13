#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${1:-$ROOT_DIR/fixtures/generated}"

if ! command -v samtools >/dev/null 2>&1; then
  echo "error: samtools is required to build the BAM fixture" >&2
  exit 2
fi

mkdir -p "$OUT_DIR"
cp "$ROOT_DIR/core/data/rcrs.fasta" "$OUT_DIR/rcrs.fasta"
samtools faidx "$OUT_DIR/rcrs.fasta"
samtools view -b "$ROOT_DIR/fixtures/truth_mixed.sam" \
  | samtools sort -o "$OUT_DIR/truth_mixed.bam" -
samtools index "$OUT_DIR/truth_mixed.bam"
samtools quickcheck -v "$OUT_DIR/truth_mixed.bam"

samtools view -C -T "$OUT_DIR/rcrs.fasta" \
  -o "$OUT_DIR/truth_mixed.cram" "$ROOT_DIR/fixtures/truth_mixed.sam"
samtools index "$OUT_DIR/truth_mixed.cram"
samtools quickcheck -v "$OUT_DIR/truth_mixed.cram"

echo "created $OUT_DIR/truth_mixed.bam"
echo "created $OUT_DIR/truth_mixed.bam.bai"
echo "created $OUT_DIR/truth_mixed.cram"
echo "created $OUT_DIR/truth_mixed.cram.crai"
