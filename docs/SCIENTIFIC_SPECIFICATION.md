# Scientific Analysis Specification

Status: **RC1 draft for independent review**
Applies to: engine 0.5.0-dev; frozen compatibility schema 0.5 and the explicitly
identified opt-in RC2 schema 0.6 slice
Normative terms: **MUST**, **MUST NOT**, **SHOULD**, and **MAY** have their
usual requirements meaning.

This document freezes the behavior already implemented by the analysis core.
It is not an analytical-validation report and does not establish clinical
performance. Any behavior marked `open` is excluded from the frozen RC1
contract until its implementation, positive/negative golden tests, and review
are complete.

Schema 0.5 is the compatibility baseline, not the final product model. The
roadmap deliberately replaces its implicit feature-token representation with
schema 0.6 alignment fragments, protocol-aware molecules, normalized events,
callable observations, phase links, and candidate architectures. The target
entities and scientific positioning are defined in
[`SCIENTIFIC_POSITIONING.md`](SCIENTIFIC_POSITIONING.md); new semantics become
normative only in a versioned successor to this specification.
The implemented but unfrozen RC2 contract is reviewed separately in
[`SCHEMA_0.6_DRAFT.md`](SCHEMA_0.6_DRAFT.md).

## 1. Coordinate and reference model

- External mtDNA coordinates are one-based and inclusive.
- The bundled reference is rCRS `NC_012920.1`, exactly 16,569 bases. A custom
  reference MUST be a single-record FASTA containing only A, C, G, T, or N.
- Reference-consuming alignment operations wrap modulo the reference length.
  Position 16,570 therefore normalizes to position 1.
- SNP alleles are reported on the reference-forward strand. SAM SEQ is consumed
  according to the SAM alignment contract; the engine MUST NOT apply a second
  reverse-complement transformation to a reverse-strand record.
- A deletion interval is inclusive: `start`, `end`, and
  `length = end - start + 1` for a non-wrapping event.
- An insertion is anchored at the current reference cursor and its length is
  the number of query bases consumed by the insertion.

Open before RC1 freeze: canonical normalization of complex events containing
multiple origin crossings or equivalent split-segment representations.

## 2. Input records and molecule identity

- FASTQ, SAM, BAM, and CRAM are accepted. BAM/CRAM require an htslib-enabled
  build. FASTA is reference-only and MUST fail as read input.
- FASTQ records require a leading `@`, a `+` separator, a non-empty sequence,
  and equal sequence/quality lengths.
- Alignment records require a valid header/record, numeric FLAG/POS/MAPQ, and a
  CIGAR whose query-consuming length equals SEQ when SEQ is present.
- For alignment input, records with the same QNAME are one molecule. The first
  record without secondary (`0x100`) or supplementary (`0x800`) flags is the
  primary; if none exists, the first record is used.
- In schema 0.6, every original record MUST remain a first-class alignment
  fragment. The representative decision MUST expose the identity policy,
  primary-candidate count, ambiguity status, and `NO_PRIMARY_ALIGNMENT` or
  `MULTIPLE_PRIMARY_ALIGNMENTS` warning when applicable. These warnings do not
  prove biological independence or automatically exclude a molecule.
- Supplementary records are folded into the primary molecule as SA evidence.
  An `SA` tag MUST use SAM type `Z` and contain one or more semicolon-terminated
  entries with exactly `RNAME,POS,STRAND,CIGAR,MAPQ,NM`. `RNAME` MUST be
  non-empty, `POS` MUST be positive, `STRAND` MUST be `+` or `-`, `CIGAR` MUST
  contain non-zero operations from the supported `MIDNSHP=X` grammar, `MAPQ`
  MUST be in `[0,255]`, and `NM` MUST be a non-negative integer.
- Every input-record `SA` tag is parsed before molecule filtering and parallel
  feature extraction. Any malformed entry fails the whole analysis with
  `MITO-E1103`; the primary alignment is not used to return a partial
  rearrangement or NUMT conclusion.
  Each molecule contributes at most one A/C/G/T allele to locus depth. If the
  same molecule covers a locus repeatedly with conflicting bases, that
  molecule is non-callable at that locus and MUST NOT expose a read-level SNP,
  contribute a clustering SNP token, or appear in that allele's
  `supporting_reads`.
- Empty inputs and malformed mandatory or supplementary-alignment records MUST
  fail; no partial result is returned. The fail-closed `SA` grammar is pinned by
  `fixtures/negative/error_manifest.json`.

## 3. SNP eligibility and evidence

An aligned molecule is eligible for SNP evidence only when all are true:

1. unmapped flag `0x4` is clear;
2. `(FLAG & excluded_snp_flags) == 0` (default `0xF00`: secondary, QC-fail,
   duplicate, and supplementary);
3. MAPQ is at least `min_mapping_quality` (default 20);
4. POS is non-zero and RNAME is present;
5. RNAME is not recognized as a nuclear chromosome;
6. the base is A/C/G/T with a valid Phred+33 quality at least
   `min_base_quality` (default 10).

For each observed alternate allele at locus `p`:

- `AD_alt` is the number of eligible molecules supporting that alternate;
- `AD_ref` is the number supporting the reference base;
- `AD_other` is the number supporting another A/C/G/T base;
- `DP_callable = AD_alt + AD_ref + AD_other`;
- `HF = AD_alt / DP_callable` when `DP_callable > 0`;
- `CI95` is the two-sided Wilson score interval with
  `z = 1.959963984540054`, clipped to `[0, 1]`.

Each accepted molecule/allele contributes one representative base-quality and
read-position observation. If the same molecule covers a locus repeatedly with
the same base, the highest-Phred observation is selected; equal qualities use
the smallest original query index. A conflicting repeated base makes the
molecule non-callable at that locus.

- `allele_quality.{alternate,reference,other}` reports count, mean, minimum,
  and maximum Phred quality; empty groups use JSON `null`, not zero-valued
  pseudo-observations.
- `strand_support` reports forward/reverse molecule counts from primary SAM
  flag `0x10`, separately for alternate, reference, and other alleles.
- `strand_bias_delta` is the absolute difference between alternate and
  reference forward-strand fractions. It is `null` unless both groups exist.
- Normalized read center proximity is
  `2 * min(q, L - 1 - q) / (L - 1)` for zero-based original query index `q`
  and sequence length `L`; a one-base sequence uses zero.
- `read_position.bias_delta` is the absolute alternate/reference difference of
  mean normalized center proximity and is `null` unless both groups exist.
- `molecule_support` repeats the alternate/reference/other/callable counts in a
  named object so consumers do not reinterpret read counts as molecule counts.

### 3.1 Opt-in RC2 sparse evidence projection

When schema 0.6 is requested, SNV aggregation and the sparse evidence projection
are generated from the same per-molecule observations across eligible primary
and supplementary fragments. Each materialized molecule/event pair is
`REFERENCE`, `ALTERNATE`, `EVENT_ABSENT`, `LOW_QUALITY`, or `CONFLICT`.
An absent sparse pair is normatively `NOT_CALLABLE`, never REF. The native core
MUST verify that alternate/reference/other/callable counts and supporting
molecule IDs exactly reproduce the aggregate variant projection before it
returns JSON. The store MUST stop without a partial result when its configured
observation bound is exceeded. Pairwise phase projection MUST independently
fail when its declared bound is exceeded.

Schema 0.6 now covers SNVs, CIGAR small indels, simple SVs, complex paths,
base-callable intervals, and candidate event-pair phase links. Negative evidence
for an insertion or deletion requires an uninterrupted callable reference path
within one alignment; inversion/duplication/ambiguous/origin/complex events are
support-only until reviewed rules exist. `SUPPORT_CONDITIONED` links MUST NOT be
interpreted as whole-sample phase fractions. `architectures[]` remains an
unimplemented projection, not biological absence.

Every serialized phase link MUST be reproducible from the same sparse
observations: its four jointly callable cells, uncertain count, sorted co-ALT
and uncertain `molecules[].index` references resolvable to exact IDs, Wilson
interval, declared independent-marginals
expectation, linkage delta, assessability and QC flags are exact invariants.
Candidate pairs are ALT-driven and canonical by event order; a missing or extra
pair is a contract failure. These statistics describe observed molecule-level
co-occurrence only. They are not calibrated phase accuracy, false-linkage
probability, clone identity or cell-population evidence.

These are transparent observational QC values. They are not calibrated
probabilities and do not filter or relabel a variant in RC1.

Every alternate at a multi-allelic locus gets a separate variant record with
the same locus-callable depth. Deleted/skipped bases, N, missing/invalid quality,
and conflicting repeated molecule observations do not contribute to callable
depth. The current caller reports observed evidence; it does not yet impose a
validated LoB/LoD/LoQ filter.

## 4. Coverage

- Coverage uses passed, mapped, mitochondrial molecules with a CIGAR.
- Only `M`, `=`, and `X` contribute depth. `D` and `N` advance the reference
  cursor but do not contribute depth. Insertions and clipping do neither.
- Circular spans are decomposed into complete reference cycles plus a bounded
  difference-array remainder.
- `mean_depth` is the arithmetic mean of exact site depths;
  `pct_sites_gt20x` uses strict `depth > 20`; `max_depth` is the exact maximum.
- Display-bin depth is the site mean rounded by adding 0.5 and truncating. It
  MUST NOT replace exact site metrics in scientific decisions.

## 5. NUMT evidence

- Competitive nuclear-plus-mitochondrial alignment evidence is required for an
  assessable NUMT-specificity claim.
- Primary nuclear alignment scores 0.99; supplementary nuclear alignment scores
  at least 0.90. Development heuristics (read-name, extreme GC, excessive
  length) are reported but MUST NOT be presented as validated specificity.
- A molecule is filtered only when filtering is enabled and
  `numt_score > numt_threshold`. Equality does not filter.
- With mtDNA-only input, the result MUST state that NUMT specificity is not
  assessable. Absence of nuclear evidence MUST NOT be reported as absence of a
  NUMT.

Open before RC5: scores and thresholds are engineering defaults, not calibrated
probabilities. Competitive mapping leakage and true-mtDNA loss require blinded
validation.

## 6. Structural rearrangements

- CIGAR `D`/`N`, `I`, and terminal/internal `S` at or above `sv_min_length`
  produce deletion, insertion, and soft-clip evidence.
- Split evidence is ordered by query interval. Opposite-strand adjacent
  segments imply an inversion junction; a positive reference gap implies a
  deletion; a reference overlap implies a duplication; and a positive query
  gap between otherwise adjacent segments implies an insertion.
- Query intervals are normalized to the original molecule orientation. For a
  reverse-strand SAM/SA segment, leading alignment clipping is transformed to
  the corresponding original-query interval before segments are ordered.
- SV event schema `1.0` uses 1-based rCRS coordinates. A deletion interval is
  the inclusive missing reference interval. An insertion anchor is the genomic
  left flank and its length is the unaligned query gap. A duplication interval
  is the inclusive reference overlap. Inversion coordinates are the sorted
  novel-adjacency breakpoints, not the full aligned blocks.
- A same-strand transition against linear reference order is classified using
  the shorter of the circular reference gap and the local reference overlap.
  The former is `circular_origin`; the latter is `duplication`. An exact tie on
  an even-length custom reference is `ambiguous_adjacency`, never guessed.
- Canonical IDs are `deletion:start-end`, `insertion:left_anchor+length`,
  `soft_clip_{left|right}:breakpoint+length`,
  `duplication:start-end`, `inversion:left_breakpoint-right_breakpoint`, and
  `circular_origin:high_breakpoint-low_breakpoint`.
- Equivalent CIGAR, forward-split, and reverse-split events merge under one ID.
  Evidence sources, strand transitions, and supporting molecule IDs are sorted
  and unique. A coordinate/length collision under one ID is an internal error,
  not a best-effort merge.
- A soft clip explained by an emitted split event is not also emitted as an
  independent soft-clip candidate.
- When one molecule has at least two emitted split-alignment junctions, the
  junction IDs form an ordered complex-event path. The lexicographically
  smaller of the oriented forward path and its reverse-complement traversal is
  canonical, so equivalent opposite-strand molecule observations merge without
  discarding junction order or strand transitions.
- Complex-event schema `1.0` uses the exact collision-free ID
  `complex:<junction-id>@<orientation>|<junction-id>@<orientation>[|...]`.
  Supporting molecule IDs are sorted and unique. The component junctions
  remain independently present in `svs`; coalescing never replaces or rewrites
  SV edge schema `1.0`.
- Complex-event paths participate in molecule clustering in addition to their
  component edge tokens. Thus two molecules with the same individual edges but
  a different supported traversal order are not treated as having an identical
  structural signature.
- A complex event is an observed strand-invariant junction path, not proof of a
  closed circular architecture. The current canonicalization does not rotate
  paths, infer unobserved junctions, or combine a CIGAR-only edge with split
  edges.

Open before RC1 freeze: obtain independent review of the circular shortest-path
and complex-path rules. Edge schema 1.0 behavior is pinned by
`truth_split.expected.svs.json`; complex schema 1.0 forward/reverse equivalence
is pinned by `truth_complex.expected.svs.json`. These development fixtures are
not independent biological truth.

## 7. Clustering, haplogroups, and clinical knowledge

- The bounded default clusterer uses unique SNP/SV/complex-path tokens per
  passed molecule, deterministic ordering, and Jaccard distance. A read with no
  variant token receives a length-bucket token.
- Optional HDBSCAN is not part of the frozen production contract until its
  dense-memory bound and label normalization are specified and benchmarked.
- Haplogroup results are ranked evidence against the pinned PhyloTree rCRS 17.3
  bundle. Resource version, source, license, retrieval date, and SHA-256 MUST be
  recorded.
- The supported pinned-tree mutation grammar is: substitution `POSITIONBASE`,
  deletion `POSITIONd`, insertion `POSITION.INDEXBASES`, wildcard insertion
  ordinal `POSITION.XBASES`, and a trailing `!` backmutation. Positions are
  1-based rCRS coordinates, numeric insertion indices are positive, and bases
  are canonical uppercase A/C/G/T. Any unsupported `<poly>` token fails the
  resource load with `MITO-E1302`; it is never silently omitted.
- A backmutation removes the inherited state at the same reference locus or
  insertion ordinal. Substitutions and deletions therefore compete for one
  reference locus, while independent insertion ordinals can coexist.
- `.X` means an unspecified insertion ordinal in the pinned PhyloTree; it is
  not heteroplasmy notation. A wildcard homopolymer marker such as `573.XC`
  matches an observed C-only insertion such as `573.1CCCC`.
- Production alignment evidence contributes substitutions plus CIGAR `D` and
  `I` events up to 50 bp. Inserted bases must all pass the configured base-
  quality threshold. Larger events remain structural-variant evidence and do
  not masquerade as haplogroup-defining small indels.
- Per-molecule callable haplogroup ranges are derived from passing aligned
  A/C/G/T bases and supported small deletions. Cluster markers and callable
  positions require a strict majority (`>50%`); a 50/50 tie produces no
  consensus state. Expected PhyloTree markers outside the consensus callable
  range are neither matched nor reported missing.
- Ranking uses the pinned mutation weights and the weighted Kulczynski form.
  The score is evidence similarity, not a posterior probability, and is not
  clamped to 100 because PhyloTree wildcard weights can make the reference
  implementation exceed 1.0.
- Heteroplasmy is represented by molecule-level clusters and aggregate allele
  evidence, not by changing the PhyloTree token grammar.
- The committed differential corpus pins HaploGrep 3.3.2, PhyloTree rCRS 17.3,
  tree, weight, alignment-rule, and HSD hashes, exact tested ranges, and exact
  top-three classifications for 28 profiles across 23 declared lineage groups,
  including three explicit inherited-marker removals by backmutation. Equal scores
  use PhyloTree source order, followed by name only as a defensive total-order
  fallback.
- Alignment-derived insertions and deletions up to 50 bp are normalized to the
  most 3-prime repeat-equivalent rCRS coordinate without crossing the 16,569/1
  boundary. Inserted motifs are rotated while shifting. Aligned-SAM goldens
  prove that `516d 517d` and an insertion after 567 converge to `523d 524d`
  and an insertion after 573, respectively, and that `8270.1ACCCCCTCT`
  converges to `8289.1CCCCCTCTA`.
- After strict-majority cluster markers and callable ranges are derived, the
  exact checksum-pinned PhyloTree 17.3 alignment rules are applied in source
  order until the marker set converges. Every error-side marker must be
  present; then that set is removed and parseable replacement markers are
  inserted. A numeric replacement is an explicit restoration to rCRS and
  emits no variant marker. Rules containing an ambiguity code that the
  callable A/C/G/T path cannot emit are inapplicable. An empty, truncated,
  malformed, unsupported, or non-convergent rules resource fails with
  `MITO-E1301`/`MITO-E1302`; no unnormalized partial result returns. Exact
  projections pin `309T 310C -> 309d`, the mixed 188-196 rule, and the 56-66
  reference-restoration rule.
- The broader HSD/IUPAC/transition-shorthand input grammar is not claimed
  equivalent. Official-rule projection is a software differential, not a
  population accuracy or mapper-normalization validation.
- Clinical assertion schema 1.0 preserves source, source record ID, allele ID,
  disease, original and normalized significance, review status, assertion
  date, source URL, references, resource version, and retrieval date for every
  row. Assertions are never merged field-wise or overwritten.
- Conflict status is deterministic: explicit source conflicts or incompatible
  pathogenic/uncertain/benign/other-assertion groups emit `conflicting`.
  Unclassified records remain visible but do not alone manufacture a conflict.
  The normative contract is
  [`CLINICAL_ANNOTATION_SCHEMA_1.0.md`](CLINICAL_ANNOTATION_SCHEMA_1.0.md).
- The bundled rows are a checksum-pinned development subset and MUST NOT be
  described as a complete clinical knowledge base. A custom override is
  labeled unpinned/unverified rather than inheriting bundled provenance.

## 8. Determinism and failure behavior

- Normalized input, configuration, reference, and resources MUST yield byte-
  identical canonical JSON at every supported thread count.
- Sorting and tie-breaking MUST be explicit; hash-container iteration MUST NOT
  determine output order.
- Arithmetic overflow, allocation failure, invalid configuration, malformed
  input/reference/resource, missing dependencies, cancellation, and internal
  invariant failure MUST terminate without a partial scientific result.
- External callers MUST branch on the stable codes in
  [`ERROR_CONTRACT.md`](ERROR_CONTRACT.md), not diagnostic text.

## 9. Review and traceability

RC1 closes only after every normative statement has a positive and negative
golden test, every open item is either closed or explicitly excluded, and two
independent scientific reviewers sign a dated review record. Source code and
development fixtures alone cannot satisfy that review gate.
