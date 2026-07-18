# Production Acceptance Programme

Mito-Architect remains **research-use only** until every mandatory gate below
has evidence attached to a tagged release. Passing software tests is necessary
but is not equivalent to analytical or clinical validation.

## Phase A: Calling correctness

Implemented in the current working tree:

- SNP evidence is restricted by configurable MAPQ, base quality, and SAM flags.
- Heteroplasmy uses locus-callable A/C/G/T molecule depth, not total reads.
- Aggregate variants include alternate, reference, other, and callable depths;
  Wilson 95% intervals; and supporting molecule IDs.
- VCF DP/HF uses the aggregate contract and records the confidence interval.
- Primary nuclear and nuclear `SA` alignments contribute to NUMT filtering.

Release evidence required:

- Unit truth for circular wraparound, reverse-strand records, deletions spanning
  a queried locus, repeated coverage of a locus by one molecule, missing base
  qualities, multi-allelic loci, and every excluded SAM flag.
- Differential comparison with an established caller on the same alignments.
- Deterministic results at 1, 2, 4, 8, and 16 worker threads.

## Phase B: Resource-backed analysis

Mandatory implementation:

- Competitive mapping to a versioned nuclear plus mitochondrial reference;
  mtDNA-only mapping is not acceptable evidence for NUMT specificity.
- Split/supplementary alignment graph reconstruction for deletions, insertions,
  duplications, inversions, and circular-origin events.
- PhyloTree-backed haplogroup assignment with ranked alternatives, missing
  diagnostic sites, contamination warnings, and resource-version provenance.
- Versioned ClinVar/MITOMAP resource snapshots with retrieval time, checksum,
  license/source metadata, normalization rules, and conflict preservation.

## Phase C: Analytical validation

Use independently characterized reference materials and blinded mixtures.
Evaluate at minimum 0%, 0.5%, 1%, 2%, 5%, 10%, 25%, 50%, 75%, and 100%
heteroplasmy across multiple depths, read-quality strata, tissues, operators,
library preparations, and flow-cell/basecaller versions.

Pre-register release thresholds for:

- SNP precision, recall, false positives per mtDNA genome, bias, repeatability,
  reproducibility, and the 95% limit of blank/detection/quantification.
- SV event precision/recall, breakpoint error, length error, molecule support,
  and complex-event reconstruction accuracy.
- NUMT leakage and true-mtDNA loss under competitive mapping.
- Haplogroup top-1/top-k accuracy and contamination detection.

Thresholds must be selected and signed by the responsible laboratory; they are
not inferred from the included development fixtures.

## Phase D: Performance and operations

- Establish peak RSS, CPU time, output expansion, and cancellation latency for
  representative 1 GB, 10 GB, and expected maximum inputs.
- Replace in-memory job state/results with durable metadata and object storage.
- Add authenticated authorization, tenant/sample isolation, audit logs, TLS at
  the ingress, secret rotation, backup/restore tests, and retention controls.
- Add structured logs, metrics, traces, health/readiness endpoints, alerting,
  and tested recovery from worker crash, disk exhaustion, corrupt upload,
  dependency outage, and interrupted shutdown.
- Produce SBOMs, pinned images, vulnerability scans, signed artifacts, migration
  tests, rollback instructions, and disaster-recovery evidence.

## Phase E: Release decision

A release may be labelled production-ready only when:

1. all mandatory gates above have linked evidence;
2. no known correctness blocker is open;
3. the validation report identifies intended use and limitations;
4. scientific, security, and operational owners sign the release; and
5. the UI and exported reports show resource versions and analysis parameters.

Clinical use additionally requires the laboratory's applicable accreditation,
quality-management, regulatory, and reporting processes. The repository cannot
grant that status by itself.
