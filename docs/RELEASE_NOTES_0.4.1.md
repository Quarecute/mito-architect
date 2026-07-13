# Release Notes 0.4.1

Version 0.4.1 is a deterministic native-workflow and performance hardening
release. The result schema remains 0.4. It is research-use only and does not
replace independent analytical validation.

## Native workflow

- Adds `mito-cli doctor` with engine/schema, compiled htslib capability, and
  samtools/minimap2/bcftools/SRA Toolkit version checks.
- Makes BAM and reference-backed CRAM generated fixtures mandatory in the full
  verification gate when the native toolchain is installed.
- Adds bcftools interoperability checks and a standards-compliant VCF contig
  declaration.
- Exposes engine and htslib capability in the server readiness response and
  refuses to start the full server with a stale no-htslib core.

## Correctness and determinism

- Disables synthetic `snp=` and `sv=` header controls by default; they require
  an explicit development-only flag.
- Exposes and validates the NUMT evidence threshold through C++, versioned C
  ABI, Rust, CLI, server, result metadata, and TypeScript.
- Adds multi-allelic depth, circular-origin SNP, malformed CIGAR, development-
  control, and threshold tests.
- Verifies canonical output identity at 1, 2, 4, 8, and 16 requested threads.
- Adds ASan/UBSan and optional HDBSCAN-C++ adapter build/test gates.

## Performance and UI

- Replaces per-cluster full PhyloTree scans with a mutation inverted index and
  deterministic consensus-assignment cache.
- On the bundled synthetic Release benchmark, 200 reads improved from about
  3.68 s to about 1.58 s on the development host; 1,000 reads completed in
  about 4.50 s. These are development measurements, not a production envelope.
- Removes placeholder multi-sample, single-cell, and timeline tabs. A module is
  now exposed only when a real input schema, computation, tests, and failure
  states exist.

## Verification environment

The completed local gate used samtools/htslib 1.24, minimap2 2.30, bcftools
1.24, and SRA Toolkit 3.4.1. ASan and UBSan passed; LeakSanitizer was disabled
because the execution sandbox uses ptrace, which LSan rejects.
