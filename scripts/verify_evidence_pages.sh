#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PAGE_DIR="${MITO_EVIDENCE_PAGE_DIR:-$ROOT_DIR/fixtures/generated/evidence-pages}"
MAX_STORAGE_BYTES="${MITO_EVIDENCE_PAGE_MAX_STORAGE_BYTES:-5000000}"
if [[ ! "$MAX_STORAGE_BYTES" =~ ^[0-9]+$ ]] || ((MAX_STORAGE_BYTES == 0)); then
  echo "error: MITO_EVIDENCE_PAGE_MAX_STORAGE_BYTES must be a positive integer" >&2
  exit 2
fi
mkdir -p "$PAGE_DIR"

cargo run --quiet -p mito-cli --offline -- analyze \
  --input fixtures/truth_phase_evidence.sam --evidence-graph \
  --evidence-page-size 3 --evidence-pages-dir "$PAGE_DIR" --json >/dev/null

jq -e '
  .schema_version == "1.0"
  and .result_schema_version == "0.6"
  and .evidence_encoding.layout == "paged_columnar_molecule_event"
  and .evidence_encoding.observation_page_size == 3
  and (.pages | length) == .evidence_encoding.observation_page_count
  and ([.pages[].count] | add) == .evidence_encoding.observation_count
' "$PAGE_DIR/manifest.json" >/dev/null

while IFS=$'\t' read -r relative_path expected_sha256; do
  actual_sha256="$(sha256sum "$PAGE_DIR/$relative_path" | awk '{print $1}')"
  if [[ "$actual_sha256" != "$expected_sha256" ]]; then
    echo "error: evidence page checksum mismatch for $relative_path" >&2
    exit 1
  fi
done < <(jq -r '.pages[] | [.path, .sha256] | @tsv' "$PAGE_DIR/manifest.json")

manifest_storage_bytes="$(jq '[.pages[].bytes] | add // 0' "$PAGE_DIR/manifest.json")"
actual_storage_bytes="$(jq -r '.pages[].path' "$PAGE_DIR/manifest.json" \
  | while IFS= read -r relative_path; do stat -c '%s' "$PAGE_DIR/$relative_path"; done \
  | awk '{ total += $1 } END { print total + 0 }')"
if [[ ! "$manifest_storage_bytes" =~ ^[0-9]+$ ]] ||
   [[ "$manifest_storage_bytes" != "$actual_storage_bytes" ]]; then
  echo "error: evidence page storage accounting mismatch" >&2
  exit 1
fi
if ((actual_storage_bytes > MAX_STORAGE_BYTES)); then
  echo "error: evidence pages use ${actual_storage_bytes} bytes, above ${MAX_STORAGE_BYTES}" >&2
  exit 1
fi

echo "Deterministic paged evidence sidecars verified: storage_bytes=$actual_storage_bytes budget_bytes=$MAX_STORAGE_BYTES"
