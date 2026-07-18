# Result Schema 0.5 Migration

Status: development contract for engine `0.5.0-dev`; not a released RC1 or a
laboratory-validated format.

## Why the schema changed

Result schema `0.5` adds a required molecule-level representation for observed
multi-junction rearrangement paths. Keeping schema `0.4` would make strict
consumers unable to distinguish the new contract from the 0.4.1 payload.

SV edge schema `1.0` is unchanged. Existing entries under `svs` retain their
canonical IDs and derivation. Complex-event schema `1.0` is an additional
layer over those edges.

## Required additions

- `metadata.complex_sv_event_schema_version`: currently `"1.0"`.
- `metadata.clinical_annotation_schema_version`: currently `"1.0"`.
- `complex_events`: always an array; empty when no molecule carries at least two
  emitted split-alignment junctions.
- `reads[].complex_event_ids`: always an array in 0.5 engine output.
- `clusters[].complex_event_signature`: always an array in 0.5 engine output.
- `reads[].haplogroup_markers`: canonical molecule-level PhyloTree evidence.
- `reads[].haplogroup_range_known`: whether an explicit callable-range contract
  was derived for the molecule.
- `reads[].haplogroup_callable_ranges`: sorted, disjoint inclusive ranges.
- `clusters[].haplogroup_assignment.observed_markers`: strict-majority marker
  consensus used by the scorer.
- `clusters[].haplogroup_assignment.callable_ranges`: strict-majority inclusive
  ranges used to decide which expected markers can be matched or missing.
- `variants[].molecule_support`, `strand_support`, `strand_bias_delta`,
  `allele_quality`, and `read_position`: observational molecule/allele QC with
  explicit `null` when a bias comparison is not estimable.
- Annotated aggregate `variants[]` expose
  `annotation.{schema_version,conflict_status,consensus_significance,assertions}`;
  every assertion retains source and snapshot fields. `reads[].snps` retain
  only the normalized `position/ref/alt` molecule key and do not duplicate
  gene, consequence, structure, or clinical payloads. Resolve the key against
  `variants[]`. See
  [`CLINICAL_ANNOTATION_SCHEMA_1.0.md`](CLINICAL_ANNOTATION_SCHEMA_1.0.md).

Each `complex_events[]` item contains:

- an exact canonical `id`;
- `junction_count` and `segment_count`;
- `canonicalization = "strand_invariant_path"`;
- ordered `junction_ids` and index-aligned `junction_orientations`;
- sorted unique `supporting_reads` molecule IDs.

## Consumer migration

1. Reject unknown major/minor result schemas at strict scientific boundaries.
2. Add `complex_events` to deserialization and preserve unknown future fields.
3. Resolve each `junction_ids[]` value against `svs[].id`; do not reinterpret
   breakpoints independently.
4. Treat a complex event as an observed junction path, not proof of a closed
   circular genome or a fully assembled molecule.
5. Use `reads[].complex_event_ids` for molecule support and
   `clusters[].complex_event_signature` for cluster summaries.
6. Do not interpret a haplogroup score as a probability. Preserve observed
   markers and callable ranges with every exported assignment; markers outside
   the callable range are an invalid producer state.
7. Treat `strand_bias_delta` and `read_position.bias_delta` as observational
   metrics, not pass/fail fields. No RC1 threshold is validated.
8. Use JSON rather than VCF when full clinical assertion provenance is needed;
   VCF contains only summary tags and a conflict flag.

The TypeScript interfaces keep these additions optional so the current UI can
still open archived schema-0.4 results. The native 0.5 producer, CLI export
validator, and server completion validator require the new top-level array.

## Compatibility evidence

- `fixtures/truth_split.expected.svs.json` pins the unchanged simple SV edge
  contract and an empty complex-event array.
- `fixtures/truth_complex.expected.svs.json` pins orientation-aware
  forward/reverse-complement coalescing.
- `fixtures/haplogroup/haplogrep3-3.3.2.expected.json` pins ordered top-three
  classifications for 28 profiles/23 lineage groups, three explicit inherited-
  marker removals, and five projection cases covering 3-prime equivalent
  placement, motif rotation, and compound official-rule rewrites.
- `bash scripts/verify.sh` validates both goldens and deterministic output at
  1, 2, 4, 8, and 16 requested threads.
