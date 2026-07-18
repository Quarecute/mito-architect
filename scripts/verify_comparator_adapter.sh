#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

OUT_DIR="${MITO_COMPARATOR_VERIFY_DIR:-$ROOT_DIR/fixtures/generated/comparator-adapter}"
mkdir -p "$OUT_DIR"

RESULT="$OUT_DIR/truth-phase.result.json"
REPORT="$OUT_DIR/truth-phase.differential.json"
REPORT_REPEAT="$OUT_DIR/truth-phase.repeat.differential.json"
REPORT_ALL="$OUT_DIR/truth-phase.all-filters.differential.json"
REPORT_EMPTY="$OUT_DIR/truth-phase.empty-comparator.differential.json"
COMPARATOR="fixtures/comparator/mtdna-server-2-v2.1.16.truth-phase.tsv"

cargo run --quiet -p mito-cli --offline -- analyze \
  --input fixtures/truth_phase_evidence.sam \
  --evidence-graph --json >"$RESULT"

cargo run --quiet -p mito-cli --offline -- compare-mt-dna-server2 \
  --result "$RESULT" --comparator "$COMPARATOR" \
  --sample truth_phase_evidence --output "$REPORT" \
  --min-call-concordance 1 --max-mean-hf-delta 0.000000001

cargo run --quiet -p mito-cli --offline -- compare-mt-dna-server2 \
  --result "$RESULT" --comparator "$COMPARATOR" \
  --sample truth_phase_evidence --output "$REPORT_REPEAT"

cmp "$REPORT" "$REPORT_REPEAT"
jq -e '
  .schema_version == "1.0"
  and .report_type == "variant_callset_differential"
  and .mtdna_server_2.version == "2.1.16"
  and .mtdna_server_2.format_contract == "mtdna-server-2-v2.1.16-variants-annotated"
  and .mtdna_server_2.filter_policy == "pass_only"
  and .metrics.matched == 5
  and .metrics.mito_architect_only == 0
  and .metrics.mtdna_server_2_only == 0
  and .metrics.call_concordance == 1
  and .metrics.mean_absolute_hf_delta <= 0.000000001
  and (.matched | length) == 5
  and (.provenance.mito_architect_result.sha256 | test("^[0-9a-f]{64}$"))
  and (.provenance.comparator_result.sha256 | test("^[0-9a-f]{64}$"))
' "$REPORT" >/dev/null

cargo run --quiet -p mito-cli --offline -- compare-mt-dna-server2 \
  --result "$RESULT" --comparator "$COMPARATOR" \
  --sample truth_phase_evidence --filter-policy all --output "$REPORT_ALL"
jq -e '
  .mtdna_server_2.filter_policy == "all"
  and .metrics.matched == 5
  and .metrics.mtdna_server_2_only == 1
  and .metrics.union == 6
  and .metrics.call_concordance > 0.833333333
  and .metrics.call_concordance < 0.833333334
' "$REPORT_ALL" >/dev/null

cargo run --quiet -p mito-cli --offline -- compare-mt-dna-server2 \
  --result "$RESULT" \
  --comparator fixtures/comparator/mtdna-server-2-v2.1.16.empty.tsv \
  --sample truth_phase_evidence --output "$REPORT_EMPTY"
jq -e '
  .metrics.matched == 0
  and .metrics.mito_architect_only == 5
  and .metrics.mtdna_server_2_only == 0
  and .metrics.union == 5
  and .metrics.call_concordance == 0
  and .metrics.mean_absolute_hf_delta == null
' "$REPORT_EMPTY" >/dev/null

if cargo run --quiet -p mito-cli --offline -- compare-mt-dna-server2 \
  --result "$RESULT" \
  --comparator fixtures/comparator/mtdna-server-2-v2.1.16.empty.tsv \
  --output "$OUT_DIR/empty-without-sample.json" >/dev/null 2>&1; then
  echo "error: header-only comparator input was accepted without --sample" >&2
  exit 1
fi

if cargo run --quiet -p mito-cli --offline -- compare-mt-dna-server2 \
  --result "$RESULT" --comparator "$COMPARATOR" \
  --sample truth_phase_evidence --filter-policy all \
  --min-call-concordance 1 --output "$OUT_DIR/expected-gate-failure.json" \
  >/dev/null 2>&1; then
  echo "error: discordant comparator callset unexpectedly passed the concordance gate" >&2
  exit 1
fi
test -s "$OUT_DIR/expected-gate-failure.json"

if cargo run --quiet -p mito-cli --offline -- compare-mt-dna-server2 \
  --result "$RESULT" \
  --comparator fixtures/comparator/mtdna-server-2-v2.1.16.invalid-coverage.tsv \
  --output "$OUT_DIR/invalid-coverage.json" >/dev/null 2>&1; then
  echo "error: invalid comparator coverage accounting was accepted" >&2
  exit 1
fi

if cargo run --quiet -p mito-cli --offline -- compare-mt-dna-server2 \
  --result "$RESULT" \
  --comparator fixtures/comparator/mtdna-server-2-v2.1.16.invalid-level.tsv \
  --output "$OUT_DIR/invalid-level.json" >/dev/null 2>&1; then
  echo "error: invalid comparator allele fraction was accepted" >&2
  exit 1
fi

if cargo run --quiet -p mito-cli --offline -- compare-mt-dna-server2 \
  --result "$RESULT" \
  --comparator fixtures/comparator/mtdna-server-2-v2.1.16.multi-sample.tsv \
  --output "$OUT_DIR/ambiguous-sample.json" >/dev/null 2>&1; then
  echo "error: multi-sample comparator input was accepted without --sample" >&2
  exit 1
fi

if cargo run --quiet -p mito-cli --offline -- compare-mt-dna-server2 \
  --result "$RESULT" --comparator "$COMPARATOR" \
  --comparator-version 2.1.15 --sample truth_phase_evidence \
  --output "$OUT_DIR/unpinned-version.json" >/dev/null 2>&1; then
  echo "error: unpinned comparator format version was accepted" >&2
  exit 1
fi

echo "Pinned mtDNA-Server 2 differential adapter verified"
