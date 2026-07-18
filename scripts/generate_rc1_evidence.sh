#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT=""
VERIFICATION_LOG=""
ALLOW_DIRTY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --output)
      OUTPUT="${2:-}"
      shift 2
      ;;
    --verification-log)
      VERIFICATION_LOG="${2:-}"
      shift 2
      ;;
    --allow-dirty)
      ALLOW_DIRTY=1
      shift
      ;;
    *)
      echo "usage: $0 --output PATH [--verification-log PATH] [--allow-dirty]" >&2
      exit 2
      ;;
  esac
done

if [[ -z "$OUTPUT" ]]; then
  echo "error: --output is required" >&2
  exit 2
fi

for tool in git jq sha256sum stat sort; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "error: $tool is required to generate an RC1 evidence index" >&2
    exit 2
  fi
done

cd "$ROOT_DIR"
DIRTY=0
if [[ -n "$(git status --porcelain=v1 --untracked-files=all)" ]]; then
  DIRTY=1
fi
if [[ "$DIRTY" -eq 1 && "$ALLOW_DIRTY" -ne 1 ]]; then
  echo "error: refusing to create immutable RC1 evidence from a dirty working tree" >&2
  exit 1
fi
if [[ "$DIRTY" -eq 0 && -z "$VERIFICATION_LOG" ]]; then
  echo "error: an immutable RC1 evidence index requires --verification-log" >&2
  exit 2
fi
if [[ -n "$VERIFICATION_LOG" && ! -f "$VERIFICATION_LOG" ]]; then
  echo "error: verification log does not exist: $VERIFICATION_LOG" >&2
  exit 2
fi

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
FILES_NDJSON="$TMP_DIR/files.ndjson"
TOOLS_NDJSON="$TMP_DIR/tools.ndjson"
CANONICAL="$TMP_DIR/files.canonical.tsv"

mapfile -d '' FILES < <(
  git ls-files --cached --others --exclude-standard -z | sort -z
)
if [[ "${#FILES[@]}" -eq 0 ]]; then
  echo "error: repository evidence file set is empty" >&2
  exit 1
fi

for path in "${FILES[@]}"; do
  if [[ "$path" == *$'\n'* || "$path" == *$'\r'* || "$path" == *$'\t'* ]]; then
    echo "error: evidence paths must not contain control delimiters: $path" >&2
    exit 1
  fi
  if [[ ! -f "$path" ]]; then
    echo "error: evidence path is not a regular file: $path" >&2
    exit 1
  fi
  digest="$(sha256sum "$path" | awk '{print $1}')"
  size="$(stat -c '%s' "$path")"
  mode="$(stat -c '%a' "$path")"
  jq -cn \
    --arg path "$path" \
    --arg sha256 "$digest" \
    --argjson size "$size" \
    --arg mode "$mode" \
    '{path:$path,sha256:$sha256,size_bytes:$size,mode:$mode}' >> "$FILES_NDJSON"
  printf '%s\t%s\t%s\t%s\n' "$path" "$digest" "$size" "$mode" >> "$CANONICAL"
done

record_tool() {
  local name="$1"
  shift
  local version="missing"
  if command -v "$1" >/dev/null 2>&1; then
    version="$("$@" 2>&1 | sed -n '1p')"
  fi
  jq -cn --arg name "$name" --arg version "$version" \
    '{name:$name,version:$version}' >> "$TOOLS_NDJSON"
}

record_tool cmake cmake --version
record_tool cargo cargo --version
record_tool rustc rustc --version
record_tool node node --version
record_tool npm npm --version
record_tool samtools samtools --version
record_tool bcftools bcftools --version
record_tool minimap2 minimap2 --version

FILES_JSON="$(jq -s '.' "$FILES_NDJSON")"
TOOLS_JSON="$(jq -s '.' "$TOOLS_NDJSON")"
FILE_SET_SHA256="$(sha256sum "$CANONICAL" | awk '{print $1}')"
COMMIT="$(git rev-parse HEAD)"
TREE="$(git rev-parse HEAD^{tree})"
SOURCE_DATE_EPOCH="$(git show -s --format=%ct HEAD)"
GENERATED_AT_UNIX="$(date +%s)"
VERIFICATION_SHA256=""
if [[ -n "$VERIFICATION_LOG" ]]; then
  VERIFICATION_SHA256="$(sha256sum "$VERIFICATION_LOG" | awk '{print $1}')"
fi

mkdir -p "$(dirname "$OUTPUT")"
TMP_OUTPUT="${OUTPUT}.tmp-$$"
jq -n \
  --arg schema_version "1.0" \
  --arg evidence_kind "rc1-candidate" \
  --arg commit "$COMMIT" \
  --arg tree "$TREE" \
  --argjson dirty "$([[ "$DIRTY" -eq 1 ]] && echo true || echo false)" \
  --argjson source_date_epoch "$SOURCE_DATE_EPOCH" \
  --argjson generated_at_unix "$GENERATED_AT_UNIX" \
  --arg file_set_sha256 "$FILE_SET_SHA256" \
  --arg verification_log "${VERIFICATION_LOG:-}" \
  --arg verification_sha256 "$VERIFICATION_SHA256" \
  --argjson files "$FILES_JSON" \
  --argjson tools "$TOOLS_JSON" \
  '{
    schema_version:$schema_version,
    evidence_kind:$evidence_kind,
    repository:{
      commit:$commit,
      tree:$tree,
      dirty:$dirty,
      source_date_epoch:$source_date_epoch
    },
    generated_at_unix:$generated_at_unix,
    file_set_sha256:$file_set_sha256,
    verification:{
      log:(if $verification_log == "" then null else $verification_log end),
      sha256:(if $verification_sha256 == "" then null else $verification_sha256 end)
    },
    tools:$tools,
    files:$files
  }' > "$TMP_OUTPUT"
mv "$TMP_OUTPUT" "$OUTPUT"

echo "RC1 evidence index: $OUTPUT"
echo "files=${#FILES[@]} file_set_sha256=$FILE_SET_SHA256 dirty=$DIRTY"
