# Release Notes 0.2.0

Mito-Architect 0.2.0 is a production-readiness hardening release. It keeps the
project focused on its own integrated mtDNA workbench contract: native analysis,
job control, reproducible fixtures, VCF export, web inspection, and operational
packaging.

## Major Changes

- Aligned-read SNP evidence now uses SAM/BAM/CRAM alignment context instead of
  synthetic hash-derived sampling.
- Native analysis jobs can be cancelled cooperatively through the server, Rust
  FFI, C API, and C++ engine.

## Medium Changes

- Result JSON now carries explicit schema and engine version metadata.
- CLI and server validate required result-contract fields before export or job
  completion.
- CLI VCF export reports deterministic SNP support, depth, and heteroplasmy
  fraction.
- Reproducible truth fixtures cover SNP, NUMT filtering, SV signatures, MAPQ,
  SAM auxiliary tags, and VCF output.
- Release-mode native benchmark target is part of the verification gate.
- Docker build can be included in the local verification gate with
  `MITO_VERIFY_DOCKER=1`.

## Smaller Hardening Items

- Upload cancellation state is tracked separately from job metadata.
- Server cleanup removes stale job and cancellation records.
- Frontend polling handles terminal error/cancelled states.
- TypeScript metadata types include the schema version.
- Docker packaging includes native data, web assets, CLI, and server binaries.
- `.gitignore` and `.dockerignore` keep generated artifacts out of source while
  preserving required web assets.
- User setup guidance documents Docker-first and native dependency paths.
- Fixture failures report concrete missing or mismatched expected fields.
- DBSCAN neighbor expansion avoids duplicate queue growth.
- HTML/JSON report paths validate malformed native results earlier.
- Benchmark output includes read count, elapsed time, and JSON size.
- Current scientific and operational boundaries are documented as internal
  roadmap items.

## Migration Notes

- Consumers should treat `metadata.schema_version` as required in result JSON.
- VCF export currently represents SNP support from passed reads; indel and
  multi-sample VCF semantics remain roadmap work.
- Cancellation is cooperative and is observed between major native analysis
  steps.
