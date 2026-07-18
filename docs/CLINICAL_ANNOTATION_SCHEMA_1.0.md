# Clinical Annotation Schema 1.0

Status: frozen software contract for the 0.5 RC1 development line. It defines
lossless source handling and deterministic conflict reporting; it does not make
the bundled development subset complete or clinically validated.

## 1. Invariants

- One TSV row is one source record for exactly one 1-based rCRS SNV.
- Original source significance is never overwritten by another row.
- Empty identifiers, dates, URLs, or review states remain empty. The engine
  does not guess them.
- Variant-level structure fields must agree across rows for the same allele;
  contradictory non-empty mappings fail with `MITO-E1302`.
- Malformed headers, field counts, alleles, positions, residues, and URLs fail
  the entire analysis with `MITO-E1302`; no partially annotated result returns.
- Assertions, sources, diseases, and references are emitted in explicit sorted
  order, independent of TSV row order and thread count.

## 2. TSV contract

The exact tab-separated header is:

```text
position ref alt gene consequence protein residue structure_id structure_chain structure_residue structure_complex source assertion_id allele_id disease clinical_significance review_status assertion_date source_url references resource_version retrieved_at
```

The displayed spaces above represent tabs in the real file. There are exactly
22 fields per row:

| Field | Meaning |
| --- | --- |
| `position`, `ref`, `alt` | 1-based rCRS SNV key; A/C/G/T only |
| `gene` .. `structure_complex` | Optional curated functional and structure mapping; not a source assertion |
| `source` | Required source identity, for example `ClinVar` or `MITOMAP` |
| `assertion_id` | Source record identifier; empty when unavailable |
| `allele_id` | Source allele identifier, distinct from the assertion record |
| `disease` | Disease/phenotype exactly as represented by the source snapshot |
| `clinical_significance` | Original source classification, preserved verbatim |
| `review_status` | Source review state; it is not converted into confidence |
| `assertion_date` | Source-provided evaluation/assertion date |
| `source_url` | Optional HTTP(S) record URL |
| `references` | Semicolon-separated source references |
| `resource_version` | Snapshot or resource version for this row |
| `retrieved_at` | Snapshot retrieval date/time supplied by the producer |

## 3. Normalization and conflicts

The engine adds `normalized_significance` while retaining the original text.
Recognized normalized values are:

`pathogenic`, `likely_pathogenic`,
`pathogenic_or_likely_pathogenic`, `uncertain_significance`, `likely_benign`,
`benign`, `benign_or_likely_benign`, `risk_factor`, `association`,
`drug_response`, `protective`, `affects`, `conflicting`, `not_provided`, and
`other`.

For conflict detection these values are grouped into pathogenic, uncertain,
benign, other-assertion, conflicting, and unclassified groups.

- One record gives `conflict_status = single_assertion`.
- Multiple compatible records give `conflict_status = consistent`.
- An explicit conflicting record, or records spanning incompatible classified
  groups, gives `conflict_status = conflicting` and
  `consensus_significance = conflicting`.
- `not_provided` and unknown `other` values are preserved but do not by
  themselves create a conflict.

`consensus_significance` is an indexing/display summary, not a replacement for
`assertions[]`. Scientific review must inspect source records, review status,
dates, resource versions, and conflicts.

## 4. Result JSON

`metadata.clinical_annotation_schema_version` is `"1.0"`. An annotated entry in
aggregate `variants[]` has:

```json
{
  "annotation": {
    "schema_version": "1.0",
    "conflict_status": "consistent",
    "consensus_significance": "pathogenic",
    "assertions": [
      {
        "source": "ClinVar",
        "assertion_id": "",
        "allele_id": "",
        "disease": "...",
        "clinical_significance": "pathogenic",
        "normalized_significance": "pathogenic",
        "review_status": "not_provided",
        "assertion_date": "",
        "source_url": "",
        "references": [],
        "resource_version": "...",
        "retrieved_at": "..."
      }
    ]
  }
}
```

`pathogenicity`, `phenotype`, `references`, `sources`, `clinvar_allele_id`, and
`mitomap_url` remain compatibility summaries. New consumers use
`consensus_significance`, `conflict_status`, and `assertions[]`.

`reads[].snps` intentionally remains a compact `position/ref/alt` molecule key
and does not repeat derived gene, consequence, structure, or clinical fields.
Resolve that key against `variants[]`; this keeps large molecule results bounded
without losing source information.

VCF exports `CLNSIG`, `CLNSRC`, and the `CLNCONFLICT` flag. VCF cannot carry the
full assertion model; JSON is normative for clinical provenance.

## 5. Resource update behavior

`mito-cli update-clinical` validates schema 1.0 before atomically replacing a
cache. `--clinvar-live` projects GRCh38 MT/M single-nucleotide rows from
ClinVar `variant_summary` and writes record, allele, disease, significance,
review, date, URL, reference, version, and retrieval fields separately.

ClinVar `variant_summary` is an aggregate tabular resource, not a lossless set
of individual submitter SCV assertions. A governed production bundle must use
the licensed/approved source representation selected by the laboratory and
validate its own normalization match rate and conflict retention.

An explicit `MITO_CLINICAL_ANNOTATIONS` override is reported as
`external-unpinned` / `not-verified` in `metadata.resources` unless a future
signed resource-bundle mechanism verifies it. This prevents a custom cache
from being mislabeled as the bundled checksum-pinned subset.

### Governed snapshot lifecycle

`mito-cli clinical-snapshot` provides the software portion of resource
governance; it does not grant a data licence or approve clinical content.

- `stage` requires an immutable operator-selected snapshot ID, a non-placeholder
  licence ID, and a source-policy ID. It validates all 22 fields, records row/
  variant/source counts, copies exact bytes into a private temporary directory,
  revalidates them, records SHA-256, fsyncs data/manifest/directory, and publishes
  the directory with one atomic rename. Existing IDs are never overwritten.
- `activate` verifies schema, governance labels, summary, checksum, timestamp,
  and optional maximum age before atomically replacing `state.json`. The prior
  active ID is retained as the rollback target.
- `verify` and `status` recompute the data checksum and summary. A future-dated,
  stale, malformed, missing, or corrupted snapshot fails closed.
- `rollback` validates the previous snapshot first, then atomically swaps active
  and previous IDs. Failed validation leaves the active state unchanged.

The active path is returned as JSON and must be supplied through
`MITO_CLINICAL_ANNOTATIONS`. The result currently labels any external path as
unpinned because signed bundle-to-engine attestation remains RC5 resource-governance work; the
snapshot manifest is nevertheless the local governance/audit evidence.

```bash
cargo run -p mito-cli -- clinical-snapshot stage \
  --store /srv/mito/clinical-snapshots \
  --source /approved/clinical_annotations.tsv \
  --snapshot-id 2026-07-approved \
  --license-id LAB-CLINICAL-RESOURCE-2026 \
  --source-policy laboratory-approved-v1 --activate

cargo run -p mito-cli -- clinical-snapshot status \
  --store /srv/mito/clinical-snapshots --max-age-days 30
```

## 6. Golden evidence

`fixtures/clinical/conflicting_assertions.tsv` contains two intentionally
incompatible synthetic sources. Its expected JSON proves source, record ID,
allele ID, disease, significance, review status, date, URL, references,
snapshot version, retrieval date, ordering, and conflict status. The canonical
verification gate runs it through the native core via the CLI.
`scripts/verify_clinical_snapshot.sh` additionally proves staging, activation,
freshness, two-snapshot rollback, unsafe-ID rejection, and checksum corruption
failure without changing the user's cache.
