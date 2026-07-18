# Result schema 0.6 draft — unified molecule/evidence graph

Status: **RC2 evidence model plus RC3 projection engineering draft, opt-in, not frozen**
CLI: `mito-cli analyze --evidence-graph --json ...`
HTTP multipart: `evidence_graph=true` in any field order
Default for non-negotiating callers: schema 0.5

Schema 0.6 preserves schema-0.5 compatibility fields where meaningful, replaces
the SNP-only `variants[]` table with a unified SNV/small-indel projection, and adds deterministic,
bidirectional traceability between source alignment fragments, assembled
molecules, normalized events, observations, base callability, and pairwise
physical linkage. It is research-use software and does not establish analytical
or clinical validity.

## 1. Top-level contract

Schema 0.6 adds:

```text
evidence_encoding
alignments[]
molecules[]
callability[]
events[]
observation_pages[]
phase_links[]
architectures[]
```

`variants[]` is now a first-class projection over every `SNV`,
`SMALL_INSERTION`, and `SMALL_DELETION` event. It is not an independently
executed caller.

The current evidence scope is `snv_indel_sv_complex_evidence_rc2`:

- passing aligned SNVs;
- CIGAR-derived small insertions and deletions shorter than `sv_min_length`;
- normalized simple SVs from the existing `svs[]` projection;
- strand-invariant complex paths from `complex_events[]`;
- callable-aware event-pair links.

`architectures[]` remains empty because candidate architecture inference,
assignment uncertainty, subsampling stability, and abundance confidence are not
implemented yet. Empty means *projection unavailable*, not *no biological
architectures exist*.

## 2. Stable IDs and deterministic order

- alignment IDs: `alignment:<source-record-index>`;
- molecule IDs: source-policy identifier described below;
- SNVs: `snv:<position>:<REF>:<ALT>`;
- small insertions: `indel:insertion:<anchor>:<inserted-sequence>`;
- small deletions: `indel:deletion:<start>-<end>:<deleted-sequence>`;
- simple structural events: `sv:<canonical-legacy-sv-id>`;
- complex events: the collision-free canonical ordered path ID;
- observations: `observation:<global-page-offset + row-index>`;
- phase links: `phase:<event-a-id>|<event-b-id>`.

Array order is derived only from source order, ordered maps, and normalized IDs.
It must not depend on pointer values, hash iteration, worker completion order, or
thread count.

## 3. Fragment-to-molecule policy

Every input record occurs exactly once in `alignments[]`. Raw sequence and
quality strings are not duplicated in JSON.

| Input state | Identity | Compatibility representative | Evidence eligible |
| --- | --- | --- | ---: |
| one FASTQ record | `fastq_record_proxy` | that record | yes |
| repeated FASTQ name | unique `#record:<index>` suffix plus warning | that record | yes |
| SAM group with one primary candidate | `sam_qname` | unique primary | yes |
| no primary candidate | `sam_qname` | first record plus `NO_PRIMARY_ALIGNMENT` | no |
| multiple primary candidates | `sam_qname` | first primary plus `MULTIPLE_PRIMARY_ALIGNMENTS` | no |
| configured explicit ID tag, one primary per source QNAME | `sam_tag:<TAG>` | deterministic highest-MAPQ primary | yes |
| missing/conflicting configured ID tag | retained `UNASSIGNED:`/`CONFLICT:` ID plus reason | deterministic audit representative | no |
| conflicting configured UMI/duplex values | explicit metadata conflict | deterministic audit representative | no |

All records remain inspectable even when the molecule is not evidence eligible.
Evidence from NUMT-filtered, missing-primary, or multiple-primary molecules is
retained at fragment/callability level but is excluded from variant, event, and
phase projections; `evidence_eligible` records that decision explicitly and the
pairwise phase projection uses `evidence_eligible_only`.

`molecule_id_tag` is the only option that changes physical-molecule grouping.
`umi_tag` and `duplex_tag` preserve and validate declared protocol metadata but
do not silently turn an RX/DX value into a molecule ID. Multiple QNAMEs may be
joined only by the explicitly configured identity tag. Conflicts and missing
required IDs are retained with `exclusion_reasons[]`; partial tags inherited
within one declared group generate warnings. `protocol_flags[]` records
multi-QNAME, supplementary, secondary, SAM-duplicate, multi-reference,
origin-spanning, and length-based concatemer candidates. Candidate flags are
auditable observations, not validated biological classifications. QNAME and
FASTQ policies remain read proxies, not proof of independent original mtDNA
molecules.

Each `molecules[]` row also exposes representative query length, base/MAPQ
quality, NUMT score/evidence, cluster ID, exact alternate event IDs, and a
five-state evidence-count summary. Those counts and IDs are checked against the
observation pages; they are navigation projections, not duplicated evidence.

## 4. Base callability

`callability[]` is the normative source for deciding whether absence may be
interpreted. Each molecule contains:

- `known`, `status`, `callable_bases`, and `callable_fraction`;
- inclusive, normalized reference `ranges[]`;
- per-alignment ranges and status;
- reference-position exclusion counts;
- insertion and soft-clip query-base counts kept in separate units.

An alignment contributes passing bases when it is mitochondrial, mapped, not
secondary/QC-fail/duplicate, above MAPQ, CIGAR-bearing, and has canonical query
bases with valid Phred scores above the configured threshold. Supplementary
fragments may contribute; secondary alignments never do. Deletions, reference
skips, non-canonical bases, invalid/missing qualities, and uncovered positions
are not base-callable.

Ranges are merged for molecule summaries, but absence of an insertion or
deletion is proven only by an uninterrupted callable reference path within one
alignment fragment. The engine does not join unrelated fragments to invent a
reference adjacency.

## 5. Events and observations

Every `events[]` record declares coordinates/alleles where applicable,
normalization, source projection, component event IDs, support IDs, evidence
counts, an explicit `negative_evidence_rule`, and one of:

- `REFERENCE_AND_ALTERNATE`: the current implementation can emit conservative
  reference/absence evidence;
- `ALTERNATE_SUPPORT_ONLY`: missing observations cannot establish absence.

The sparse store is keyed by `(molecule_id,event_id)`. A missing pair is
`NOT_CALLABLE`, declared once in `evidence_encoding`; it is never materialized.

| Observation state | Meaning | Callable denominator |
| --- | --- | ---: |
| `REFERENCE` | passing SNV base equals REF | yes |
| `ALTERNATE` | normalized event is supported | yes |
| `EVENT_ABSENT` | event was assessable and another/reference state was proven | yes |
| `LOW_QUALITY` | aligned locus exists but base eligibility failed | no |
| `CONFLICT` | passing fragments/segments disagree within the molecule | no |
| missing pair (`NOT_CALLABLE`) | event cannot be assessed | no |

Observations carry molecule, event and alignment IDs, evidence source, MAPQ,
strand, and allele/base-quality/read-position fields when applicable. SNVs are
merged across eligible primary and supplementary fragments; deterministic
highest-quality evidence is retained when alleles agree, while disagreement is
`CONFLICT`.

Small-indel absence requires callable reference adjacency (insertion) or a
callable deleted span with both flanks (deletion) in one alignment. Simple
alternative indel at the same normalized type/anchor proves `EVENT_ABSENT` for
the queried allele through `cigar_alternative_small_indel`, but is not mislabeled
as reference-path evidence. Simple
insertion/deletion SVs use the same conservative rule. Inversion, duplication,
ambiguous adjacency, origin topology, soft clips, and complex paths are support-
only. Their rule is `support_only_no_negative_inference`, so a missing pair can
never be promoted to an event-absent observation. Implemented SNV/insertion/
deletion rules name their callable-base, same-fragment adjacency, or
same-fragment span-and-flanks requirement directly in each event.

## 6. Callable-aware phase links

`phase_links[]` is generated only for event pairs where at least one side has an
`ALTERNATE` observation on an evidence-eligible molecule. Counts partition jointly callable
molecules into:

```text
both_alternate
a_alternate_b_absent
a_absent_b_alternate
neither_alternate
```

Their sum equals `jointly_callable`. `jointly_uncertain` is separate.
`supporting_molecule_indices` contains exactly the sorted `molecules[].index`
references counted in `both_alternate`; `uncertain_molecule_indices` contains
exactly the sorted references counted in `jointly_uncertain`. The encoding is
declared by `evidence_encoding.phase_molecule_reference="molecules[].index"`.
Consumers resolve original IDs through the molecule table; repeating long IDs
inside every pair is forbidden. `qc_flags` is empty for a complete,
fully callable comparison, contains `SUPPORT_CONDITIONED` when either event has
no negative-evidence model, and contains
`UNCERTAIN_COOCCURRENCE_EXCLUDED` when uncertain molecules were excluded from
the denominator.
`co_alternate_fraction`, Wilson 95% interval,
`expected_co_alternate_fraction`, and `linkage_delta = P(AB) - P(A)P(B)` are
observational statistics, not calibrated biological confidence or a
pathogenicity claim. The declared null model is
`evidence_encoding.phase_null_model =
independent_marginals_within_jointly_callable`; it is stored once rather than
repeated in every link. Changing it is a schema change.

`COMPLETE_FOR_BOTH_EVENTS` means both event types have implemented negative-
evidence rules. `SUPPORT_CONDITIONED` means at least one side is support-only;
the denominator must not be interpreted as the whole sample.

Pair generation is two-pass: identify ALT-driven candidate pairs, then count all
qualifying states for those pairs. This preserves prior
REF/ABSENT observations without materializing the full event-by-event matrix.

## 7. Fail-closed invariants and limits

Before serialization, the native core verifies:

1. every fragment resolves to exactly one molecule;
2. molecule, event, alignment, and component references resolve;
3. IDs are unique and indices contiguous;
4. callable ranges are ordered, non-overlapping, in-reference, and their length
   equals `callable_bases`;
5. `(molecule,event)` occurs at most once;
6. `NOT_CALLABLE` is not stored explicitly;
7. unified SNV/small-indel AD/DP, uncertainty counts, and supporting IDs equal
   the observation projection;
8. every event support list equals its `ALTERNATE` observation set;
9. each molecule evidence summary and alternate event list round-trips to the
   observations;
10. every ALT-driven event pair occurs exactly once in canonical event order;
11. every phase four-state partition, supporting/uncertain molecule list and
    QC flag set round-trips exactly to the observation store;
12. linkage delta, independent-marginals expectation, co-ALT fraction and
    Wilson interval are recomputed and checked at every strict consumer.

`max_evidence_observations`, `evidence_page_size`, and `max_phase_links` are
independent hard bounds.
Exceeding any bound returns `MITO-E1601`; partial JSON is never returned. CLI
flags and server environment limits are recorded in run metadata.

## 8. Paged observation transport

Observations are serialized as deterministic columnar pages rather than a list
of repeated JSON objects. Each page declares `index`, global `offset`, `count`,
and equal-length columns for molecule/event/alignment IDs, state, allele,
qualities, strand, evidence source, and read position. Page indices and offsets
must be contiguous and their total must equal
`evidence_encoding.observation_count`.

The CLI can write checksum-addressed sidecars with
`--evidence-pages-dir <directory>`. The server exposes an authenticated
`GET /result/:job_id/evidence/:page_index` endpoint and a compact
`GET /result/:job_id/summary` transport projection. The summary contains no
embedded observation pages, declares `observation_storage=remote_http_pages`,
and publishes both page and search endpoints. Authenticated
`GET /result/:job_id/evidence` accepts exact optional `molecule_id`, `event_id`,
and stored `state` filters plus a zero-based match `cursor` and a `limit` from 1
to 500. `NOT_CALLABLE` is rejected because it is a sparse missing-pair state,
not a stored observation. Results are ordered by global observation index and
include global ID, physical page index, row index, and the complete observation
fields. The response reports the exact match count and next cursor.

The server serializes each immutable page once and builds an interned string
dictionary plus integer posting lists once at job completion. An exact query
starts from the smallest applicable posting list rather than parsing the full
result or scanning serialized pages. `GET /result/:job_id` and CLI JSON remain the authoritative
portable schema-0.6 result; the compact summary is a UI transport projection,
not a second scientific contract.

## 9. Unified variant compatibility exports

For schema 0.6, `variants[]` contains event ID/type, normalized coordinates and
alleles, ref/alt/other/callable molecule counts, event-absence/low-quality/
conflict counts, HF and Wilson CI, strand/base/MAPQ/read-position summaries,
homopolymer and multi-allelic context, NUMT assessability, supporting molecule
IDs, and reason-coded observational QC facts. `filter_status` remains
`NOT_CALIBRATED`; QC facts are not hidden hard filters.

The CLI derives JSON, TSV, VCF 4.3, bgzip VCF, tabix, HTML, and provenance from
this projection. VCF `FILTER=.` records that thresholds are not calibrated and
uses INFO fields for QC facts. A small indel whose lossless allele representation
crosses 16,569/1 is rejected from linear VCF export with an explicit error; it
remains available in authoritative JSON/TSV instead of being silently omitted
or rewritten into a misleading allele.

Likewise, a single origin-crossing CIGAR can be represented in SAM/BAM but not
losslessly in a linear CRAM reference slice when its reference consumption
extends past position 16569. Mito-Architect does not claim CRAM parity for that
non-linear representation. Producers must rotate or split the alignment before
CRAM encoding; ordinary BAM/CRAM parity remains a release gate on alignments
valid in both formats.

## 10. Migration and review boundary

- Consumers must negotiate `/metadata/schema_version` and reject unknown
  schemas.
- `variants[].event_id`, `svs[].event_id`, and
  `complex_events[].event_id` link legacy projections to `events[]`.
- `supporting_reads` retains its schema 0.5 compatibility name and currently
  contains assembled source molecule IDs.
- A missing observation must never be converted to REF or event absence.
- Columnar paging reduces repeated-key storage and supports bounded retrieval,
  but the native engine still creates one bounded result document. Declared
  peak-RSS, storage, cancellation, and supported-maximum evidence envelopes are
  required before very large workloads are accepted as production-scale.

Before freeze, two domain reviewers must approve molecule identity, fragment
eligibility, circular indel normalization, event-specific absence rules, phase
statistics, and protocol metadata. Independent truth mixtures must then validate
SNV/indel/SV calling and phase accuracy; code invariants alone do not close the
scientific gate.
