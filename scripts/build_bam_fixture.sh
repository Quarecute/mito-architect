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

samtools view -b "$ROOT_DIR/fixtures/truth_phase_evidence.sam" \
  | samtools sort -o "$OUT_DIR/truth_phase_evidence.bam" -
samtools index "$OUT_DIR/truth_phase_evidence.bam"
samtools quickcheck -v "$OUT_DIR/truth_phase_evidence.bam"
samtools view -C -T "$OUT_DIR/rcrs.fasta" \
  -o "$OUT_DIR/truth_phase_evidence.cram" \
  "$OUT_DIR/truth_phase_evidence.bam"
samtools index "$OUT_DIR/truth_phase_evidence.cram"
samtools quickcheck -v "$OUT_DIR/truth_phase_evidence.cram"

samtools view -b "$ROOT_DIR/fixtures/truth_phase_matrix.sam" \
  | samtools sort -o "$OUT_DIR/truth_phase_matrix.bam" -
samtools index "$OUT_DIR/truth_phase_matrix.bam"
samtools quickcheck -v "$OUT_DIR/truth_phase_matrix.bam"
samtools view -C -T "$OUT_DIR/rcrs.fasta" \
  -o "$OUT_DIR/truth_phase_matrix.cram" \
  "$OUT_DIR/truth_phase_matrix.bam"
samtools index "$OUT_DIR/truth_phase_matrix.cram"
samtools quickcheck -v "$OUT_DIR/truth_phase_matrix.cram"

samtools view -b "$ROOT_DIR/fixtures/truth_protocol_identity.sam" \
  | samtools sort -o "$OUT_DIR/truth_protocol_identity.bam" -
samtools index "$OUT_DIR/truth_protocol_identity.bam"
samtools quickcheck -v "$OUT_DIR/truth_protocol_identity.bam"
samtools view -C -T "$OUT_DIR/rcrs.fasta" \
  -o "$OUT_DIR/truth_protocol_identity.cram" \
  "$OUT_DIR/truth_protocol_identity.bam"
samtools index "$OUT_DIR/truth_protocol_identity.cram"
samtools quickcheck -v "$OUT_DIR/truth_protocol_identity.cram"

samtools view -b "$ROOT_DIR/fixtures/truth_circular_indel.sam" \
  | samtools sort -o "$OUT_DIR/truth_circular_indel.bam" -
samtools index "$OUT_DIR/truth_circular_indel.bam"
samtools quickcheck -v "$OUT_DIR/truth_circular_indel.bam"

samtools view -b "$ROOT_DIR/fixtures/truth_repeat_rotation.sam" \
  | samtools sort -o "$OUT_DIR/truth_repeat_rotation.bam" -
samtools index "$OUT_DIR/truth_repeat_rotation.bam"
samtools quickcheck -v "$OUT_DIR/truth_repeat_rotation.bam"
samtools view -C -T "$OUT_DIR/rcrs.fasta" \
  -o "$OUT_DIR/truth_repeat_rotation.cram" \
  "$OUT_DIR/truth_repeat_rotation.bam"
samtools index "$OUT_DIR/truth_repeat_rotation.cram"
samtools quickcheck -v "$OUT_DIR/truth_repeat_rotation.cram"

echo "created $OUT_DIR/truth_mixed.bam"
echo "created $OUT_DIR/truth_mixed.bam.bai"
echo "created $OUT_DIR/truth_mixed.cram"
echo "created $OUT_DIR/truth_mixed.cram.crai"
echo "created $OUT_DIR/truth_phase_evidence.bam"
echo "created $OUT_DIR/truth_phase_evidence.bam.bai"
echo "created $OUT_DIR/truth_phase_evidence.cram"
echo "created $OUT_DIR/truth_phase_evidence.cram.crai"
echo "created $OUT_DIR/truth_phase_matrix.bam"
echo "created $OUT_DIR/truth_phase_matrix.bam.bai"
echo "created $OUT_DIR/truth_phase_matrix.cram"
echo "created $OUT_DIR/truth_phase_matrix.cram.crai"
echo "created $OUT_DIR/truth_protocol_identity.bam"
echo "created $OUT_DIR/truth_protocol_identity.bam.bai"
echo "created $OUT_DIR/truth_protocol_identity.cram"
echo "created $OUT_DIR/truth_protocol_identity.cram.crai"
echo "created $OUT_DIR/truth_circular_indel.bam"
echo "created $OUT_DIR/truth_circular_indel.bam.bai"
echo "created $OUT_DIR/truth_repeat_rotation.bam"
echo "created $OUT_DIR/truth_repeat_rotation.bam.bai"
echo "created $OUT_DIR/truth_repeat_rotation.cram"
echo "created $OUT_DIR/truth_repeat_rotation.cram.crai"
echo "note: circular-indel CRAM is intentionally not emitted; a linear CRAM"
echo "      slice cannot encode a CIGAR extending past reference position 16569"
