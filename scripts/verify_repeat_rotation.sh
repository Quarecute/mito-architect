#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

project_repeat_rotation() {
  local input="$1"
  local output="$2"
  local threads="${3:-1}"
  cargo run --quiet --manifest-path "$ROOT_DIR/Cargo.toml" -p mito-cli --offline -- \
    analyze --input "$input" --evidence-graph --threads "$threads" --json \
    | jq '{
        schema_version: .metadata.schema_version,
        events: [.events[]
          | select(.type == "SMALL_INSERTION" or .type == "SMALL_DELETION")
          | {id,type,start,end,length,ref,alt,normalization,supporting_molecule_ids,evidence_counts}],
        variants: [.variants[]
          | {event_id,alt_depth,ref_depth,event_absent_depth,callable_depth,supporting_molecule_ids}]
      }' >"$output"
}

assert_same_json() {
  local expected="$1"
  local actual="$2"
  local label="$3"
  jq -S . "$expected" >"$TMP_DIR/expected.sorted.json"
  jq -S . "$actual" >"$TMP_DIR/actual.sorted.json"
  if ! cmp -s "$TMP_DIR/expected.sorted.json" "$TMP_DIR/actual.sorted.json"; then
    echo "error: $label differs from the exact repeat/rotation projection" >&2
    diff -u "$TMP_DIR/expected.sorted.json" "$TMP_DIR/actual.sorted.json" >&2 || true
    exit 1
  fi
}

EXPECTED="$ROOT_DIR/fixtures/truth_repeat_rotation.expected.json"
project_repeat_rotation "$ROOT_DIR/fixtures/truth_repeat_rotation.sam" \
  "$TMP_DIR/canonical.json" 1
assert_same_json "$EXPECTED" "$TMP_DIR/canonical.json" "canonical SAM"

project_repeat_rotation "$ROOT_DIR/fixtures/truth_repeat_rotation_permuted.sam" \
  "$TMP_DIR/permuted.json" 4
assert_same_json "$EXPECTED" "$TMP_DIR/permuted.json" "permuted four-thread SAM"

if [[ -f "$ROOT_DIR/fixtures/generated/truth_repeat_rotation.bam" &&
      -f "$ROOT_DIR/fixtures/generated/truth_repeat_rotation.cram" ]]; then
  project_repeat_rotation "$ROOT_DIR/fixtures/generated/truth_repeat_rotation.bam" \
    "$TMP_DIR/bam.json" 2
  assert_same_json "$EXPECTED" "$TMP_DIR/bam.json" "coordinate-sorted BAM"
  project_repeat_rotation "$ROOT_DIR/fixtures/generated/truth_repeat_rotation.cram" \
    "$TMP_DIR/cram.json" 2
  assert_same_json "$EXPECTED" "$TMP_DIR/cram.json" "coordinate-sorted CRAM"
fi

cargo run --quiet --manifest-path "$ROOT_DIR/Cargo.toml" -p mito-cli --offline -- \
  analyze --input "$ROOT_DIR/fixtures/truth_circular_indel.sam" --evidence-graph --json \
  | jq '{
      events: [.events[] | {id,type,start,end,length,ref,alt,normalization,supporting_molecule_ids,evidence_counts}],
      variants: [.variants[] | {event_id,start,end,ref,alt,alt_depth,ref_depth,event_absent_depth,callable_depth,supporting_molecule_ids}]
    }' >"$TMP_DIR/circular.json"
jq '{
    events: [.events[] | {id,type,start,end,length,ref,alt,normalization,supporting_molecule_ids,evidence_counts}],
    variants: [.variants[] | {event_id,start,end,ref,alt,alt_depth,ref_depth,event_absent_depth,callable_depth,supporting_molecule_ids}]
  }' "$ROOT_DIR/fixtures/truth_circular_indel.expected.json" \
  >"$TMP_DIR/circular.expected.json"
assert_same_json "$TMP_DIR/circular.expected.json" "$TMP_DIR/circular.json" \
  "origin-crossing circular indel"

echo "Exact circular/repeat/motif-rotation projections verified across order, threads, and available alignment formats"
