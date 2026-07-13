# Release Notes 0.3.0

Version 0.3.0 is a research-use production-hardening release. It is not an
analytically or clinically validated release.

## Result contract 0.3

- Adds required aggregate `variants` with locus-callable depth, allele depths,
  heteroplasmy, Wilson 95% intervals, and supporting molecule IDs.
- Records MAPQ, base-quality, and excluded-SAM-flag calling parameters.
- VCF DP/HF now uses aggregate locus evidence; consumers relying on total-read
  denominators must migrate.

## Analysis and provenance

- Adds configurable SNP evidence filters through C++, C ABI, Rust, CLI, and
  backend environment variables.
- Uses primary nuclear and nuclear supplementary alignments as strong NUMT
  evidence when competitive mapping is supplied.
- ClinVar refreshes emit checksum/timestamp metadata and retain bundled curated
  protein-structure mappings through field-wise duplicate merging.
- Cluster labels no longer resemble inferred haplogroups.

## Web and service

- Replaces fabricated NGL fallback coordinates with provenance-bearing models
  and explicit failure states. A hardening patch corrects the human MT-ATP6
  accession mapping and makes experimental RCSB structures local-first.
- Adds molecule/cluster support beside structural residue context.
- Adds API-key protection, safe public-bind enforcement, and health/readiness
  endpoints.

## Remaining release blockers

PhyloTree assignment, split-alignment event graphs, durable job storage,
representative performance envelopes, and external analytical validation remain
mandatory before a production-ready designation. See `PRODUCTION_ACCEPTANCE.md`.
