#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STORE="$(mktemp -d)"
trap 'rm -rf "$STORE"' EXIT

CLI=(cargo run --quiet --manifest-path "$ROOT_DIR/Cargo.toml" -p mito-cli --offline -- clinical-snapshot)
SOURCE="$ROOT_DIR/core/data/clinical_annotations.tsv"

"${CLI[@]}" stage --store "$STORE" --source "$SOURCE" \
  --snapshot-id governed-a --license-id development-subset \
  --source-policy source-specific-terms --activate \
  | jq -e '.operation == "activate" and .state.active_snapshot == "governed-a" and .verification.verified'

"${CLI[@]}" stage --store "$STORE" --source "$SOURCE" \
  --snapshot-id governed-b --license-id development-subset \
  --source-policy source-specific-terms --activate \
  | jq -e '.state.active_snapshot == "governed-b" and .state.previous_snapshot == "governed-a"'

"${CLI[@]}" rollback --store "$STORE" \
  | jq -e '.operation == "rollback" and .state.active_snapshot == "governed-a"'

"${CLI[@]}" status --store "$STORE" --max-age-days 1 \
  | jq -e '.verified and .active_snapshot == "governed-a" and .verification.manifest.summary.assertion_count > 0'

"${CLI[@]}" verify --store "$STORE" --snapshot-id governed-b \
  | jq -e '.verified and .snapshot_id == "governed-b"'

truncate -s -1 "$STORE/snapshots/governed-b/clinical_annotations.tsv"
if "${CLI[@]}" verify --store "$STORE" --snapshot-id governed-b >/dev/null 2>&1; then
  echo "error: corrupted clinical snapshot unexpectedly verified" >&2
  exit 1
fi

if "${CLI[@]}" stage --store "$STORE" --source "$SOURCE" \
  --snapshot-id ../escape --license-id development-subset \
  --source-policy source-specific-terms >/dev/null 2>&1; then
  echo "error: unsafe clinical snapshot ID unexpectedly staged" >&2
  exit 1
fi

echo "Clinical snapshot staging, activation, rollback, freshness, and corruption checks passed"
