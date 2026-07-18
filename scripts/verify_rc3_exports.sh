#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

for program in jq bgzip tabix bcftools sha256sum; do
  if ! command -v "$program" >/dev/null 2>&1; then
    echo "error: $program is required for the RC3 export gate" >&2
    exit 2
  fi
done

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/mito-rc3-exports.XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT

run_export() {
  cargo run --quiet -p mito-cli --offline -- analyze \
    --input fixtures/truth_phase_evidence.sam \
    --evidence-graph \
    --vcf "$WORK_DIR/variants.vcf" \
    --tsv "$WORK_DIR/variants.tsv" \
    --bgzip-vcf "$WORK_DIR/variants.vcf.gz" \
    --provenance-manifest "$WORK_DIR/provenance.json" \
    --json >"$WORK_DIR/result.json"
}

run_export

jq -e '
  .metadata.schema_version == "0.6"
  and ([.events[] | select(.type == "SNV" or .type == "SMALL_INSERTION" or .type == "SMALL_DELETION")] | length) == (.variants | length)
  and ([.variants[].type] | index("SNV") != null)
  and ([.variants[].type] | index("SMALL_INSERTION") != null)
  and ([.variants[].type] | index("SMALL_DELETION") != null)
  and all(.variants[];
    .event_id != null
    and .callable_depth == (.alt_depth + .ref_depth + .other_depth)
    and .filter_status == "NOT_CALIBRATED"
    and (.supporting_molecule_ids | type) == "array")
' "$WORK_DIR/result.json" >/dev/null

expected_lines="$(( $(jq '.variants | length' "$WORK_DIR/result.json") + 1 ))"
actual_lines="$(wc -l <"$WORK_DIR/variants.tsv")"
if [[ "$actual_lines" != "$expected_lines" ]]; then
  echo "error: TSV row count mismatch: expected $expected_lines, got $actual_lines" >&2
  exit 1
fi

bcftools view --no-version -Ov "$WORK_DIR/variants.vcf" >/dev/null
bcftools view --no-version -Ov "$WORK_DIR/variants.vcf.gz" >/dev/null
tabix -l "$WORK_DIR/variants.vcf.gz" | grep -Fx 'NC_012920.1' >/dev/null
gzip -cd "$WORK_DIR/variants.vcf.gz" >"$WORK_DIR/variants.roundtrip.vcf"
cmp "$WORK_DIR/variants.vcf" "$WORK_DIR/variants.roundtrip.vcf"

grep -F '##INFO=<ID=AD,Number=R' "$WORK_DIR/variants.vcf" >/dev/null
grep -F '##INFO=<ID=MOLECULE_SUPPORT,Number=1' "$WORK_DIR/variants.vcf" >/dev/null
grep -F '##INFO=<ID=NUMT_ASSESSABLE,Number=1' "$WORK_DIR/variants.vcf" >/dev/null

jq -e '
  .schema_version == "1.0"
  and .software.result_schema_version == "0.6"
  and .determinism.deterministic_algorithms == true
  and .determinism.random_seed == null
  and (.input.sha256 | length) == 64
  and (.reference.sha256 | length) == 64
  and (.authoritative_result.sha256 | length) == 64
  and ([.exports[].kind] | index("vcf") != null)
  and ([.exports[].kind] | index("variant_tsv") != null)
  and ([.exports[].kind] | index("vcf_bgzip") != null)
  and ([.exports[].kind] | index("vcf_tabix") != null)
' "$WORK_DIR/provenance.json" >/dev/null

while IFS=$'\t' read -r path expected; do
  actual="$(sha256sum "$path" | awk '{print $1}')"
  if [[ "$actual" != "$expected" ]]; then
    echo "error: provenance checksum mismatch for $path" >&2
    exit 1
  fi
done < <(jq -r '.exports[] | [.path, .sha256] | @tsv' "$WORK_DIR/provenance.json")

sha256sum "$WORK_DIR/variants.vcf" "$WORK_DIR/variants.tsv" \
  "$WORK_DIR/variants.vcf.gz" "$WORK_DIR/variants.vcf.gz.tbi" \
  "$WORK_DIR/provenance.json" >"$WORK_DIR/first.sha256"
run_export
sha256sum -c "$WORK_DIR/first.sha256" >/dev/null

cargo run --quiet -p mito-cli --offline -- analyze \
  --input fixtures/truth_circular_indel.sam --evidence-graph --json \
  >"$WORK_DIR/circular.json"
jq -e '
  any(.variants[];
    .type == "SMALL_DELETION"
    and .normalization == "rcrs_circular_small_indel_v1"
    and .vcf_representable == false)
' "$WORK_DIR/circular.json" >/dev/null
if cargo run --quiet -p mito-cli --offline -- analyze \
  --input fixtures/truth_circular_indel.sam --evidence-graph \
  --vcf "$WORK_DIR/circular.vcf" --json \
  >"$WORK_DIR/circular-vcf.stdout" 2>"$WORK_DIR/circular-vcf.stderr"; then
  echo "error: circular-origin indel VCF export must fail without a lossless representation" >&2
  exit 1
fi
grep -F 'has no lossless linear VCF representation' "$WORK_DIR/circular-vcf.stderr" >/dev/null

echo "RC3 unified JSON/TSV/VCF/bgzip/tabix/provenance exports verified"
