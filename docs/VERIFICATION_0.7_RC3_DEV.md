# 0.7 / RC3 Development Verification

Verification date: 2026-07-18
Branch: `sipa/TASK-003`
Scope: final engineering slice; research-use development candidate

## Result

The final repository-implementable RC3 slice passed its local development
gate. The result now provides linked variant, molecule, callable-aware phase,
and rearrangement projections from the same schema 0.6 evidence graph.

This is not analytical or clinical validation. Independent truth materials,
two scientific reviewers, laboratory UAT, and schema disposition remain
external closure gates.

## Implemented contract

- The strict mtDNA-Server 2 v2.1.16 differential remains pinned and explicitly
  non-truth-bearing.
- Exact SNV/indel/SV phase support and uncertainty resolve through stable
  molecule indices and are recomputed by native, CLI, and server validators.
- Completed server jobs pre-serialize immutable evidence pages and build an
  interned string dictionary with integer postings for exact molecule, event,
  and stored-state search.
- The authenticated search endpoint is bounded to 500 rows, cursor-paginated,
  starts from the smallest applicable posting list, and returns stable global,
  page, and row coordinates. `NOT_CALLABLE` cannot be queried as a stored row.
- The browser provides bounded global search, current-page virtualization,
  deterministic match pagination, keyboard-operable tabs, and explicit sparse-
  absence wording.
- The Rearrangements inspector traces an ordered complex path through component
  SVs, strand transitions, supporting molecules, and retained alignment
  fragments. It does not infer missing junctions, circular closure, or a
  biological architecture.
- Exact fixtures prove origin-crossing deletion behavior and collapse left/
  right repeat representations and rotated insertion motifs to the same
  normalized events across record order, requested threads, SAM, BAM, and CRAM.

## Executed gates

| Gate | Result |
| --- | --- |
| `bash scripts/verify.sh` | PASS |
| C++ smoke | 1 passed |
| Rust CLI tests | 4 passed |
| Rust FFI tests | 4 passed |
| Rust server tests | 7 passed |
| TypeScript typecheck and production builds | PASS |
| NGL source-policy tests | 2 passed |
| Analysis navigation/filter tests | 2 passed |
| Comparator, phase, exports, clinical, haplogroup, negative and resource gates | PASS |
| SAM/BAM/CRAM repeat/rotation parity | PASS |
| `bash scripts/verify_sanitizers.sh` | PASS; ASan/UBSan smoke in 24.58 s |
| `bash scripts/verify_benchmark_budget.sh` | PASS |
| `bash scripts/verify_evidence_benchmark_budget.sh` | PASS |
| `cargo fmt --check` and `git diff --check` | PASS |

The canonical schema 0.5 output remained byte-identical at 1, 2, 4, 8, and 16
requested threads:

```text
e639394ea4986db5a3e058de851eaf700b2a3eda37e719ce790c1c47b720cab2
```

## Same-host regression measurements

These are local regression gates, not portable throughput claims:

| Profile | Observed | Fixed budget | Result |
| --- | ---: | ---: | --- |
| schema 0.5, 200 reads | 819 ms median | 1,000 ms | PASS |
| schema 0.5, 1,000 reads | 2,312 ms median | 2,800 ms | PASS |
| schema 0.5, 10,000 reads | 16,402 ms median | 21,000 ms | PASS |
| adverse schema 0.6, 200 reads | 6,736 ms median | 8,000 ms | PASS |
| schema 0.6 JSON | 54,368,423 bytes | 60,000,000 bytes | PASS |
| schema 0.6 peak RSS | 147,048 KiB | 800,000 KiB | PASS |
| deterministic evidence sidecars | 8,397 bytes | 5,000,000 bytes | PASS |

## Remaining release boundary

The engineering slice is complete, but RC3 cannot be scientifically frozen
until the two reviewers disposition schema 0.6 and independent comparator,
phase, circular/rearrangement truth, accessibility, and laboratory UAT evidence
are attached. RC4 architecture inference and competitive mapping, followed by
RC5 blinded analytical validation and publication packaging, remain subsequent
milestones.
