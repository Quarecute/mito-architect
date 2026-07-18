# Verification Record 0.4.1

Date: 2026-07-12
Host: Arch Linux x86-64, Linux zen 7.1.3, GCC 16.1.1

## Toolchain

| Tool | Verified version |
| --- | --- |
| mito-engine | 0.4.1, result schema 0.4 |
| samtools / htslib | 1.24 / 1.24 |
| minimap2 | 2.30-r1287 |
| bcftools | 1.24 |
| SRA Toolkit | 3.4.1 |
| CMake | 4.4.0 |
| Rust | 1.94.0-nightly |
| Node.js / npm | 26.4.0 / 12.0.0 |

## Passed gates

- C++ build and expanded native smoke suite with system htslib.
- Rust workspace tests and Clippy with warnings denied.
- TypeScript workspace type checking and production builds.
- Local-first RCSB resolution tests and six scientific-resource SHA-256 checks.
- SAM, generated coordinate-sorted BAM, and reference-backed CRAM truth
  comparisons against golden VCF.
- bcftools parsing of both golden VCF files without header warnings.
- Canonical JSON identity at 1, 2, 4, 8, and 16 requested threads.
- Explicit HDBSCAN-C++ adapter Release build and smoke test.
- AddressSanitizer and UndefinedBehaviorSanitizer expanded smoke test.
  LeakSanitizer was disabled because the execution sandbox uses ptrace, which
  causes LSan itself to terminate before program diagnostics.

## Development benchmark

Release build, synthetic SAM, host-local measurement:

| Reads | Elapsed | JSON bytes |
| ---: | ---: | ---: |
| 200 | 1,580 ms | 901,814 |
| 1,000 | 4,497 ms | 4,333,658 |

The earlier 200-read implementation measured about 3,675 ms on the same host.
This result demonstrates the PhyloTree index/cache improvement but is not a
production performance envelope. The 1 GB, 10 GB, and maximum-input matrix in
`PLAN_1.0.md` remains mandatory.

## Not satisfied by this record

- Blinded low-frequency and structural-variant analytical validation.
- Quantified NUMT leakage on a complete versioned competitive reference.
- A complete governed clinical resource and conflict model.
- Durable multi-tenant service state, audit, backup/restore, and recovery drills.
- Independent laboratory, security, and operations approval.
