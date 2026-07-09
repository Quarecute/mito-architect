# Production Readiness Plan

This plan defines the gates Mito-Architect must satisfy before it is treated as
production-ready. Each gate needs code, tests, docs, and reproducible command
coverage.

## Gate 1: Reproducible Validation Fixtures

Status: complete for the current offline fixture suite.

Acceptance criteria:

- The repository contains small aligned truth fixtures with expected SNP, SV,
  MAPQ, SAM auxiliary tag, filter, JSON, and VCF outcomes.
- CI/local checks can validate those fixtures without network access.
- Fixture failures report concrete missing or mismatched fields.

## Gate 2: Stable Result Contract

Status: complete for required runtime fields.

Acceptance criteria:

- JSON output has explicit schema-version metadata. Done for native output.
- CLI/server reject malformed result documents before exporting or serving.
  CLI VCF export and server job completion validate required contract fields.
- TypeScript and Rust expectations stay aligned with the JSON contract.

## Gate 3: Cancellable Native Jobs

Status: partial. Native parsing, extraction, fallback clustering, and coverage
are cancellable; an in-progress third-party HDBSCAN call is not interruptible.

Acceptance criteria:

- Server cancellation reaches the native C++ engine through the FFI.
- Cancellation is tested at the Rust boundary.
- Long file reads and third-party calls document their cancellation boundary.

## Gate 4: VCF And Heteroplasmy Export

Status: partial. SNP support export is deterministic, but DP/HF currently use
passed-read count rather than locus-specific callable depth.

Acceptance criteria:

- CLI VCF export is covered by unit and fixture tests.
- VCF support/depth/heteroplasmy semantics are documented.
- Multi-allelic or repeated-position behavior is deterministic.

## Gate 5: Performance Envelope

Status: partial. A Release benchmark target exists and reports runtime/output
size, but representative dataset envelopes and CI regression thresholds are not
yet defined.

Acceptance criteria:

- Native benchmark target runs in Release mode.
- Benchmark output includes read count, elapsed time, and output size.
- The project tracks expected ranges for fixture and synthetic benchmark runs.

## Gate 6: Operational Packaging

Status: code-complete; local Docker execution is blocked until the current user
can access `/var/run/docker.sock`.

Acceptance criteria:

- Docker build is part of the verification checklist. It is available through
  `MITO_VERIFY_DOCKER=1 bash scripts/verify.sh`.
- Native dependency setup is documented and failure modes are explicit.
- Runtime temp files, reports, and caches are isolated from source files.
- Current environment note: `docker.service` starts, but `docker build` requires
  the user session to be in the `docker` group or to run Docker through `sudo`.

## Gate 7: Production Boundary Statement

Status: complete for 0.2.0.

Acceptance criteria:

- README and user guide describe exactly which outputs are validated.
- Known limitations are framed as internal roadmap items, not comparisons.
- Release notes summarize behavior changes and migration risks.
