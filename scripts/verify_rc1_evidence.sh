#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
INDEX="$TMP_DIR/evidence.json"
CANONICAL="$TMP_DIR/files.canonical.tsv"

bash "$ROOT_DIR/scripts/generate_rc1_evidence.sh" \
  --output "$INDEX" --allow-dirty >/dev/null

jq -e '
  .schema_version == "1.0"
  and .evidence_kind == "rc1-candidate"
  and (.repository.commit | test("^[0-9a-f]{40}$"))
  and (.repository.tree | test("^[0-9a-f]{40}$"))
  and (.file_set_sha256 | test("^[0-9a-f]{64}$"))
  and (.files | length > 0)
  and (.tools | length == 8)
' "$INDEX" >/dev/null

cd "$ROOT_DIR"
while IFS=$'\t' read -r path expected_sha expected_size expected_mode; do
  if [[ ! -f "$path" ]]; then
    echo "error: evidence index references missing file: $path" >&2
    exit 1
  fi
  actual_sha="$(sha256sum "$path" | awk '{print $1}')"
  actual_size="$(stat -c '%s' "$path")"
  actual_mode="$(stat -c '%a' "$path")"
  if [[ "$actual_sha" != "$expected_sha" || "$actual_size" != "$expected_size" || "$actual_mode" != "$expected_mode" ]]; then
    echo "error: evidence mismatch for $path" >&2
    exit 1
  fi
  printf '%s\t%s\t%s\t%s\n' "$path" "$actual_sha" "$actual_size" "$actual_mode" >> "$CANONICAL"
done < <(jq -r '.files[] | [.path,.sha256,.size_bytes,.mode] | @tsv' "$INDEX")

actual_set_sha="$(sha256sum "$CANONICAL" | awk '{print $1}')"
expected_set_sha="$(jq -r '.file_set_sha256' "$INDEX")"
if [[ "$actual_set_sha" != "$expected_set_sha" ]]; then
  echo "error: RC1 evidence file-set checksum mismatch" >&2
  exit 1
fi

echo "RC1 evidence index structure and every indexed file verified"
