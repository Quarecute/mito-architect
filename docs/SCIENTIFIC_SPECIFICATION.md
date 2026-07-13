# Scientific Analysis Specification

Status: **RC1 draft for independent review**
Applies to: engine 0.4.1, result schema 0.4
Normative terms: **MUST**, **MUST NOT**, **SHOULD**, and **MAY** have their
usual requirements meaning.

This document freezes the behavior already implemented by the analysis core.
It is not an analytical-validation report and does not establish clinical
performance. Any behavior marked `open` is excluded from the frozen RC1
contract until its implementation, positive/negative golden tests, and review
are complete.

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
- Supplementary records are folded into the primary molecule as SA evidence.
  Each molecule contributes at most one A/C/G/T allele to locus depth. If the
  same molecule covers a locus repeatedly with conflicting bases, that
  molecule is non-callable at that locus and MUST NOT appear in that allele's
  `supporting_reads`.
- Empty inputs and malformed mandatory records MUST fail; no partial result is
  returned.

Open before RC1 freeze: malformed optional SA entries are currently ignored.
RC1 must choose and test either a fail-closed rule or an explicit warning and
per-segment rejection counter.

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

Open before RC2: scores and thresholds are engineering defaults, not calibrated
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

Open before RC1 freeze: coalesce adjacent normalized edges into canonical
multi-junction complex events and obtain independent review of the circular
shortest-path rule. Event schema 1.0 behavior is pinned by
`truth_split.expected.svs.json` but the development fixture is not independent
biological truth.

## 7. Clustering, haplogroups, and clinical knowledge

- The bounded default clusterer uses unique SNP/SV tokens per passed molecule,
  deterministic ordering, and Jaccard distance. A read with no variant token
  receives a length-bucket token.
- Optional HDBSCAN is not part of the frozen production contract until its
  dense-memory bound and label normalization are specified and benchmarked.
- Haplogroup results are ranked evidence against the pinned PhyloTree rCRS 17.3
  bundle. Resource version, source, license, retrieval date, and SHA-256 MUST be
  recorded.
- Current PhyloTree scoring supports the documented SNV/backmutation subset.
  Insertions, deletions, heteroplasmy notation, and the differential HaploGrep
  corpus remain open RC1 work.
- Clinical assertions MUST preserve source identity and resource provenance.
  The bundled rows are a development subset and MUST NOT be described as a
  complete clinical knowledge base.

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
