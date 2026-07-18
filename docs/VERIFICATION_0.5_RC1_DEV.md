# 0.5 / RC1 Development Verification

Date: 2026-07-13
Status: development evidence; not an RC1 release or laboratory validation

## Canonical gate

Command:

```bash
bash scripts/verify.sh
```

Result: passed.

Verified in this run:

- C++ configure/build and `mito_core_smoke_test`;
- Cargo workspace check, tests, and clippy with warnings denied;
- native capability audit: htslib/BAM/CRAM, samtools 1.24, minimap2
  2.30-r1287, bcftools 1.24, and SRA Toolkit 3.4.1;
- TypeScript typecheck and production builds for the visualization library and
  web application;
- local-first NGL structure-source tests and checksums for bundled 8H9S/9I4I
  structures;
- six versioned analysis-resource records;
- SNP and mixed-SV SAM truth fixtures;
- exact SV event-schema 1.0 JSON golden covering CIGAR/split equivalence,
  both split strands, deletion, insertion, duplication, inversion junction,
  circular-origin normalization, sorted provenance/support, and competitive
  NUMT exclusion;
- generated BAM and CRAM equivalence and bcftools VCF parsing;
- optional HDBSCAN-C++ adapter build and smoke test.

The opt-in sanitizer gate was also run separately:

```bash
bash scripts/verify_sanitizers.sh
```

ASan/UBSan configure, build, and `mito_core_smoke_test` passed. Docker was not
enabled in this invocation. All optional gates remain required again on the
exact release candidate where specified by the acceptance programme.

After reverse-strand split query coordinates were normalized to the original
molecule orientation, the complete canonical gate was rerun on the current
tree. This included the browser production build, release benchmarks, optional
HDBSCAN configuration, exact external SV golden, and 1/2/4/8/16-thread
determinism. ASan/UBSan was also rerun after the correction as a separate gate.

## Determinism

Canonical result SHA-256 was identical at 1, 2, 4, 8, and 16 worker threads:

```text
96b0b7a4e83cdcab3d7fca61fab9a63437d99a86dd4d8cbf3d7c77b69e71b4f7
```

## Release-build benchmark smoke

| Reads | Elapsed | JSON bytes |
| ---: | ---: | ---: |
| 200 | 1,653 ms | 901,902 |
| 1,000 | 4,866 ms | 4,333,746 |
| 10,000 | 33,403 ms | 43,048,530 |

These are same-host development smoke measurements, not the RC3 performance
envelope. They do not replace p50/p95 repetition, peak RSS, temporary-storage,
cancellation-latency, architecture, or 1 GB/10 GB/supported-maximum evidence.

## New RC1 assertions

- Stable native codes cross C++ -> C -> Rust without diagnostic-text parsing.
- HTTP failures use error envelope schema 1.0; internal details are not returned
  to clients; an exact test prevents server filesystem paths leaking through
  engine open failures.
- Invalid configuration, missing/malformed/empty/unsupported input, invalid or
  missing reference, and cancellation have exact negative tests.
- Reverse-strand observations are counted on reference coordinates.
- Missing-quality and default secondary/QC-fail/duplicate/supplementary records
  are excluded from SNP evidence.
- Conflicting repeated coverage by one molecule is non-callable and no longer
  leaks into an alternate allele's `supporting_reads`.
- SV event IDs use a dedicated schema version and merge equivalent CIGAR,
  forward-split, and reverse-split evidence without order-dependent provenance.
- Insertions use the genomic left-flank anchor; duplications report an interval
  consistent with length; inversions use junction breakpoints rather than
  entire aligned blocks; origin classification uses circular distance instead
  of a fixed coordinate-window heuristic.
- Explained soft clips are not double-counted as independent candidates, and
  native alignment-coordinate arithmetic is overflow-checked.
- Reverse SAM/SA clipping is transformed back to original molecule query
  coordinates before split segments are ordered.
