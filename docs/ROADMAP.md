# Development Roadmap

The detailed release-candidate plan, measurable gates, and work-item definition
of done are maintained in [`PLAN_1.0.md`](PLAN_1.0.md). This file is the compact
milestone index.

Execution status and linked evidence are maintained in
[`RELEASE_STATUS.md`](RELEASE_STATUS.md). RC1 is active; it is not yet closed.

The documentation hierarchy and evidence policy are indexed in
[`README.md`](README.md). Repository-level delivered changes are summarized in
[`../CHANGELOG.md`](../CHANGELOG.md).

This roadmap orders the remaining work by scientific risk. Mito-Architect stays
research-use only until the production acceptance programme has signed evidence
for every mandatory gate.

## 0.3.x: Correctness and viewer reliability — implemented

- Keep curated protein, chain, and residue mappings versioned and testable.
- Expand the optional offline structure cache and expose source/checksum data in
  exported reports.
- Add browser integration tests for WebGL initialization, local-cache success,
  RCSB fallback, residue selection, and total source failure.
- Finish edge-case truth tests for circular loci, reverse alignments,
  multi-allelic sites, and repeated molecule coverage.

Exit condition: the offline fixture suite and viewer tests are deterministic;
no coordinate model is accepted without organism, protein, chain, residue
range, source, and checksum provenance.

## 0.4: Resource-backed analysis contracts — implemented baseline

- Detect and preserve competitive nuclear-plus-mitochondrial alignment evidence;
  explicitly mark mtDNA-only analysis as non-assessable for NUMT specificity.
- Reconstruct split and supplementary alignments as a molecule event graph for
  deletions, insertions, duplications, inversions, and circular-origin events.
- Add PhyloTree-backed haplogroup assignment with ranked alternatives, missing
  diagnostic sites, contamination signals, and resource provenance.
- Snapshot ClinVar and MITOMAP-compatible annotations with normalization,
  conflict preservation, checksums, licenses, and retrieval timestamps.

Exit condition: every reported biological conclusion is resource-backed and
reproducible from the recorded input/reference/resource versions.

The implementation baseline is shipped in 0.4.1. First-class minimap2
orchestration, quantitative NUMT leakage validation, PhyloTree complex indel
nomenclature, and a complete governed clinical resource are carried into the
1.0 programme.

## 0.5 / RC1: Frozen scientific specification

- **Active.** The review draft and error schema are implemented; open semantic
  rules, negative goldens, differential haplogroup evidence, and two reviewer
  approvals remain.
- Canonical SV edge schema 1.0 now merges equivalent CIGAR, forward-split, and
  reverse-split events with exact external JSON evidence; multi-junction event
  coalescing and independent biological truth remain.
- Freeze coordinate, depth, HF/CI, no-call, multi-allelic, circular-origin,
  molecule identity, rearrangement, haplogroup, and clinical conflict semantics.
- Add stable error codes and positive/negative golden output for every rule.
- Review the normative scientific specification independently.

## 0.6 / RC2: Calibrated analytical performance

- Pre-register thresholds, then run blinded heteroplasmy mixtures and
  independently characterized structural-rearrangement/NUMT truth.
- Quantify LoB/LoD/LoQ, precision/recall, HF bias, breakpoint error,
  repeatability/reproducibility, leakage, and haplogroup accuracy.
- Publish and sign supported and unsupported conditions.

## 0.7 / RC3: Scale and performance envelope

- Integrate versioned minimap2 competitive mapping as a first-class stage.
- Bound streaming, clustering, memory, output, temporary storage, and
  cancellation behavior.
- Publish runtime/resource envelopes at 1 GB, 10 GB, and the supported maximum;
  pass fuzzing, sanitizer, C ABI stress, and determinism gates.

## 0.8 / RC4: Durable secure service

- Introduce PostgreSQL metadata, object storage, isolated workers, OIDC/RBAC,
  immutable audit, observability, backup/restore, and failure recovery.
- Pin and sign the supply chain; prove tenant isolation and rollback.

## 0.9 / RC5: Laboratory release candidate

- Complete batch/project/review/approval/resource-governance workflows and
  immutable interoperable exports.
- Complete accessibility, browser/WebGL, operator documentation, and laboratory
  user-acceptance evidence.

## 1.0: Operable production service

- Freeze the exact RC5 artifacts and rerun every scientific, performance,
  security, recovery, workflow, compatibility, and supply-chain gate.
- Publish the evidence index, intended use, supported envelope, limitations,
  migration/rollback instructions, and signed artifacts.
- Complete scientific, security, and operational release review against
  `PRODUCTION_ACCEPTANCE.md`; promote only that reviewed candidate.

Exit condition: every mandatory acceptance gate has linked evidence and named
scientific, security, and operational owners have approved the release.

The release-candidate sequence, deliverables, and evidence gates are specified
in [`PLAN_1.0.md`](PLAN_1.0.md).
