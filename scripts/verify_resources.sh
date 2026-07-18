#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="$ROOT_DIR/core/data/resource_manifest.tsv"

test -f "$MANIFEST"

checked=0
while IFS=$'\t' read -r name version path expected_sha source license retrieved; do
  if [[ "$name" == "name" || -z "$name" ]]; then
    continue
  fi

  if [[ "$path" == web/* ]]; then
    file="$ROOT_DIR/$path"
  else
    file="$ROOT_DIR/core/data/$path"
  fi

  if [[ ! -f "$file" ]]; then
    echo "Missing required analysis resource: $file" >&2
    exit 1
  fi

  actual_sha="$(sha256sum "$file" | awk '{print $1}')"
  if [[ "$actual_sha" != "$expected_sha" ]]; then
    echo "Checksum mismatch for $name $version: expected $expected_sha, got $actual_sha" >&2
    exit 1
  fi
  checked=$((checked + 1))
done < "$MANIFEST"

if [[ "$checked" -lt 7 ]]; then
  echo "Resource verification checked only $checked files; expected at least 7" >&2
  exit 1
fi

echo "Verified $checked versioned resources"
