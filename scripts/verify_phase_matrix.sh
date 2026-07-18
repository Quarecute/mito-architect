#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

OUT_DIR="${MITO_PHASE_MATRIX_VERIFY_DIR:-$ROOT_DIR/fixtures/generated/phase-matrix}"
mkdir -p "$OUT_DIR"
RESULT="$OUT_DIR/result.json"

cargo run --quiet -p mito-cli --offline -- validate-evidence-graph-fixture \
  --input fixtures/truth_phase_matrix.sam \
  --expected-json fixtures/truth_phase_matrix.expected.json \
  --evidence-page-size 5

cargo run --quiet -p mito-cli --offline -- analyze \
  --input fixtures/truth_phase_matrix.sam --evidence-graph --json >"$RESULT"

jq -e '
  ([.events[] | {key: .id, value: .type}] | from_entries) as $types
  | (.molecules | map(.id)) as $molecule_ids
  | any(.phase_links[];
      $types[.event_a_id] == "SNV" and $types[.event_b_id] == "SV_DELETION")
  and any(.phase_links[];
      (($types[.event_a_id] == "SMALL_INSERTION" and $types[.event_b_id] == "SV_DELETION")
       or ($types[.event_b_id] == "SMALL_INSERTION" and $types[.event_a_id] == "SV_DELETION")))
  and any(.phase_links[];
      (($types[.event_a_id] == "SV_DELETION" and $types[.event_b_id] == "SV_INSERTION")
       or ($types[.event_b_id] == "SV_DELETION" and $types[.event_a_id] == "SV_INSERTION")))
  and all(.phase_links[];
      (.supporting_molecule_indices | length) == .both_alternate
      and (.uncertain_molecule_indices | length) == .jointly_uncertain
      and (.qc_flags | type) == "array")
  and any(.phase_links[];
      .event_a_id == "snv:3:T:A"
      and .event_b_id == "sv:deletion:11-20"
      and .jointly_uncertain == 1
      and ([.uncertain_molecule_indices[] as $index | $molecule_ids[$index]] == ["lowq-snv-deletion"])
      and ([.supporting_molecule_indices[] as $index | $molecule_ids[$index]] == ["all-events"]))
' "$RESULT" >/dev/null

echo "Callable-aware SNV/indel/SV phase matrix verified"
