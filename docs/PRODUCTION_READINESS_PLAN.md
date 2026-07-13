# Production Readiness Plan

This plan defines the gates Mito-Architect must satisfy before it is treated as
production-ready. Each gate needs code, tests, docs, and reproducible command
coverage.

## Gate 1: Reproducible Validation Fixtures

Status: partial. Offline development fixtures exist; representative calibrated
materials and the matrix in `PRODUCTION_ACCEPTANCE.md` are still required.

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

Status: implementation complete, analytical validation pending. SNP aggregates
and VCF now use quality-filtered locus-callable A/C/G/T molecule depth and emit
Wilson 95% heteroplasmy intervals.

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

Status: partial. Container and bounded local execution exist, but durable job
state, authentication, audit logging, observability, and recovery testing do not.

Acceptance criteria:

- Docker build is part of the verification checklist. It is available through
  `MITO_VERIFY_DOCKER=1 bash scripts/verify.sh`.
- Native dependency setup is documented and failure modes are explicit.
- Runtime temp files, reports, and caches are isolated from source files.
- Current environment note: `docker.service` starts, but `docker build` requires
  the user session to be in the `docker` group or to run Docker through `sudo`.

## Gate 7: Production Boundary Statement

Status: complete for the 0.4.1 research-use boundary statement.

Acceptance criteria:

- README and user guide describe exactly which outputs are validated.
- Known limitations are framed as internal roadmap items, not comparisons.
- Release notes summarize behavior changes and migration risks.

## Gate 8: Structure Viewer Provenance

Status: complete for currently mapped MT-ATP6 and MT-ND4 variants. Both use
human experimental cryo-EM entries with verified chain ranges, bundled
checksum-verified coordinates, and two RCSB network formats. Only curated PDB
identifiers are accepted.

Acceptance criteria:

- Coordinate models are resolved from a named public structure source.
- Viewer interactions do not trigger repeated coordinate downloads.
- The UI exposes model provenance, residue selection, local atomic context, and
  confidence coloring.
- Offline or failed model retrieval is explicit; synthetic coordinates are not
  presented as biological structure.

## Gate 9: Biological Benchmark Cohort

Status: partial. A reproducible bounded public ONT acquisition path and offline
truth fixtures exist; calibrated low-frequency and rearrangement reference
materials are still required.

Acceptance criteria:

- Public accessions, reference sequence, checksums, and tool versions are
  captured for every downloaded validation bundle.
- CI remains independent of network and large human datasets.
- Sensitivity, precision, heteroplasmy error, breakpoint tolerance, and NUMT
  leakage are measured against molecule-level truth across relevant depths.
