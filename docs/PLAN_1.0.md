# Mito-Architect 1.0 Roadmap

## 1. Product boundary

Version 1.0 is a reproducible research instrument for molecule-level long-read
mtDNA analysis. It must call and phase heteroplasmy and structural
rearrangements, assess NUMT evidence after competitive mapping, assign
haplogroups, attach versioned clinical knowledge, and expose every conclusion
with molecule and resource provenance.

The repository can establish software quality and reproducibility. A laboratory
must supply independently characterized material, pre-register analytical
thresholds, and approve the validation report. Clinical diagnostic use also
requires the laboratory's quality system, accreditation, and jurisdictional
review.

## 2. Non-negotiable engineering rules

1. The same normalized input, configuration, reference bundle, and resource
   bundle produces the same canonical result at every supported thread count.
2. No development control, fabricated coordinate, guessed clinical assertion,
   or mtDNA-only NUMT claim is presented as biological evidence.
3. Every output records engine/schema version, parameters, input/reference
   identity, resource version and SHA-256, and limitations that affect
   interpretation.
4. Malformed data fails closed with a stable machine-readable error; partial
   scientific conclusions are not silently returned.
5. Fast paths must preserve the same golden output as scalar/reference paths.
6. A feature is exposed in the production UI only after its input contract,
   result schema, tests, and failure states exist.

## 3. Current baseline: 0.4.1 hardening

| Capability | Current state | Evidence |
| --- | --- | --- |
| FASTQ/SAM/BAM/CRAM | Implemented | parser tests; generated BAM and CRAM truth validation |
| Determinism | Implemented locally | canonical equality at 1, 2, 4, 8, and 16 threads |
| SNP/HF evidence | Implemented baseline | locus-callable DP, allele depths, HF, Wilson 95% CI |
| NUMT assessment | Implemented baseline | competitive header/evidence contract; non-assessable warning otherwise |
| Rearrangements | Implemented edge baseline | canonical event schema 1.0 and exact CIGAR/forward/reverse split JSON golden; complex coalescing remains |
| Haplogroups | Implemented baseline | PhyloTree rCRS 17.3 SNV/backmutation scorer with ranked provenance |
| Clinical annotation | Development subset | versioned local rows; complete licensed/conflict-preserving bundle required |
| 3D protein view | Implemented mapped subset | bundled RCSB models, checksum verification, curated residue mapping |
| Service operations | Development service | bounded jobs/API key; durable state, tenants, audit, restore not implemented |
| Analytical validation | Not complete | laboratory truth mixtures and signed thresholds required |

## 4. Release sequence

### 0.5 / RC1 — Frozen scientific specification

Objective: make every caller decision independently reviewable before tuning.

Current execution state: active. The normative draft is in
[`SCIENTIFIC_SPECIFICATION.md`](SCIENTIFIC_SPECIFICATION.md), error schema 1.0
is in [`ERROR_CONTRACT.md`](ERROR_CONTRACT.md), and the evidence ledger is in
[`RELEASE_STATUS.md`](RELEASE_STATUS.md). These artifacts begin RC1; they do not
close its reviewer and golden-corpus gate.

Engineering work:

- Make `scripts/bootstrap.sh --check` a required environment audit and publish
  a tested dependency matrix for every supported native and container host.
- Write normative coordinate, strand, circular-origin, molecule identity,
  callable-depth, allele-depth, HF, CI, multi-allelic, and no-call semantics.
- Replace ad-hoc error text at external boundaries with error codes and a
  versioned error schema while retaining actionable messages.
- Extend PhyloTree parsing to insertions, deletions, heteroplasmy notation, and
  the complete supported mutation grammar; document deviations from HaploGrep.
- Add strand-bias, read-position-bias, molecule-support, and per-allele quality
  metrics without converting them into unvalidated hard filters.
- Complete multi-junction complex-event coalescing on top of canonical SV edge
  schema 1.0; preserve the exact equivalence golden.
- Freeze a clinical normalization/conflict model that preserves source,
  assertion, review status, disease, date, allele ID, and source URL separately.

Required tests:

- Circular SNPs and SVs across position 16569/1, both strands, overlapping
  CIGAR operations, repeated molecule coverage, every excluded SAM flag,
  missing qualities, malformed SA tags, and multi-allelic sites.
- Golden JSON/VCF/TSV for every normative rule.
- Differential haplogroup corpus against a pinned HaploGrep release.

Exit gate:

- Scientific specification reviewed by two named reviewers.
- Every rule maps to at least one positive and one negative golden test.
- No biological output field lacks a documented derivation.

### 0.6 / RC2 — Calibrated analytical performance

Objective: quantify where the instrument is accurate and where it must refuse
or qualify interpretation.

Data design:

- Blinded mixtures at 0%, 0.5%, 1%, 2%, 5%, 10%, 25%, 50%, 75%, and 100% HF.
- Multiple depths, molecule lengths, quality strata, tissues, operators,
  extraction/library protocols, flow cells, and supported basecaller versions.
- Independently characterized deletions, insertions, duplications, inversions,
  circular-origin events, and complex mixtures.
- Competitive nuclear-plus-mitochondrial truth containing known NUMTs and true
  mtDNA molecules.

Metrics to pre-register:

| Domain | Mandatory metrics |
| --- | --- |
| SNP | precision, recall, FP/genome, HF bias, MAE, repeatability, reproducibility |
| Detection | limit of blank, 95% LoD, LoQ, reportable range |
| SV | event precision/recall, breakpoint/length error, support calibration |
| NUMT | leakage, true-mtDNA loss, low-HF false positives |
| Haplogroup | top-1/top-3 accuracy, quality calibration, contamination detection |
| Clinical | normalized match rate, conflict retention, stale-resource behavior |

Exit gate:

- Thresholds were fixed before unblinding.
- Results meet those thresholds on material independent of development fixtures.
- A signed validation report identifies supported and unsupported conditions.

### 0.7 / RC3 — Scale and performance envelope

Objective: make resource use predictable on laboratory workloads.

Engineering work:

- Integrate minimap2 competitive mapping as a first-class, versioned pipeline
  stage rather than relying on an undocumented external command.
- Stream parser-to-feature work in bounded batches; page or column-store
  molecule details instead of retaining one monolithic JSON document.
- Keep sparse token clustering as the bounded default. Permit HDBSCAN only
  behind measured memory limits; never materialize an unbounded dense matrix.
- Maintain mutation-inverted haplogroup scoring and cache identical consensus
  signatures.
- Add input/output byte counters, peak RSS, phase timings, cancellation latency,
  and bounded temporary-storage accounting.
- Profile before optimization; preserve canonical golden output after each
  optimization.

Benchmark matrix:

- 1 GB, 10 GB, and the declared maximum supported input.
- 1, 2, 4, 8, and 16 threads on each supported CPU architecture.
- Sparse and high-diversity molecule populations, high SV burden, and maximum
  auxiliary-tag payload.
- Cold/warm resource cache and local/network-isolated viewer operation.

Exit gate:

- Published p50/p95 wall time, peak RSS, temp bytes, result bytes, and cancel
  latency with regression budgets.
- Canonical output hashes match across supported thread counts.
- ASan/UBSan, parser fuzzing, and C++/C ABI stress tests have no unresolved
  correctness finding.

### 0.8 / RC4 — Durable secure service

Objective: survive process, host, and operator failures without losing
provenance or crossing sample boundaries.

Architecture:

- PostgreSQL for job/project/sample metadata and immutable audit events.
- Checksummed S3-compatible object storage for inputs, results, reports, and
  resource bundles.
- Separate API, scheduler, and least-privilege workers with idempotent state
  transitions and explicit retry/cancellation policy.
- OIDC authentication, role-based authorization, tenant/sample isolation, TLS
  ingress, secret rotation, quotas, retention, and deletion workflow.
- Structured logs, OpenTelemetry traces, Prometheus metrics, dashboards,
  dependency/readiness probes, and actionable alerts.

Failure drills:

- Worker crash and restart, API restart during upload, scheduler restart,
  corrupt object, database failover, disk exhaustion, object-store outage,
  dependency timeout, interrupted shutdown, backup restore, migration rollback,
  and full release rollback.

Supply chain:

- Pinned base images and dependencies, SBOM, license report, vulnerability scan,
  signed artifacts, provenance attestation, migration test, and rollback image.

Exit gate:

- Threat model and security review approved.
- Restore point and recovery time objectives demonstrated in a recorded drill.
- No high-severity unresolved vulnerability or tenant-isolation finding.

### 0.9 / RC5 — Laboratory workflow and release candidate

Objective: make the validated engine usable without developer intervention.

Workflow:

- Project/sample manifests, batch import, explicit reference/resource bundle,
  QC rules, review states, comments, approval, amendment, and immutable export.
- Multi-sample comparison only after a true multi-sample schema and statistics
  exist; single-cell and longitudinal views follow the same rule.
- Resource-update workflow with staging, diff, validation, approval, activation,
  rollback, and old-result reproducibility.
- JSON, VCF, TSV, and standalone HTML compatibility tests and schema migration
  tooling.
- Accessibility, keyboard operation, supported-browser/WebGL matrix, large-data
  UI profiling, and explicit offline/failure behavior.

Documentation:

- Operator guide, scientific method, validation report, intended use,
  limitations, data governance, security model, backup/restore, incident
  response, troubleshooting, and release checklist.

Exit gate:

- Representative laboratory users complete defined workflows without developer
  intervention.
- Every mandatory acceptance item links to an artifact, test, report, or drill.

## 5. 1.0 release decision

Release 1.0 only when all statements are true:

1. All RC exit gates pass on the exact release candidate.
2. No open correctness blocker or silent-data-loss defect exists.
3. Blinded validation meets pre-registered thresholds.
4. Determinism and performance envelopes are published for supported hosts.
5. Backup, restore, migration, rollback, and recovery drills pass.
6. UI and exports expose complete input, parameter, reference, resource, and
   limitation provenance.
7. Named scientific, security, and operational owners approve the release.

## 6. Critical path and execution order

Work proceeds in this order; later tracks may prototype in parallel but cannot
close before their prerequisite gate:

| Order | Release | Blocking deliverable | Evidence that closes it |
| --- | --- | --- | --- |
| 1 | 0.5 / RC1 | Frozen semantics and stable errors | normative spec, differential corpus, positive/negative goldens, two reviewers |
| 2 | 0.6 / RC2 | Calibrated scientific claims | pre-registration, blinded truth data, metric tables, signed validation report |
| 3 | 0.7 / RC3 | Bounded compute and first-class mapping | versioned competitive pipeline, benchmark matrix, fuzz/sanitizer/stress evidence |
| 4 | 0.8 / RC4 | Durable and secure operation | threat model, isolation tests, observability, backup/restore and failure drills |
| 5 | 0.9 / RC5 | Laboratory-ready workflow | user acceptance, governed resources, export compatibility, complete operator docs |
| 6 | 1.0 | Release decision | exact-RC evidence index and named scientific/security/operations approvals |

The scientific critical path is RC1 -> blinded-material acquisition -> RC2.
The operational critical path is storage/auth architecture -> RC4. Resource
licensing and validation-material procurement begin during RC1 because code
cannot compress those lead times.

For each RC we cut one immutable candidate, run its entire gate on that exact
artifact, attach hashes and raw results, and either promote it or cut a new
candidate. Evidence from an older binary does not close a newer release.

## 7. What automation can and cannot complete

The repository can automate host dependency checks, locked project dependency
resolution, compilation, deterministic fixtures, interoperability checks,
sanitizers, benchmarks, resource checksums, and release-evidence packaging.

It cannot manufacture independent truth. The following need explicit owners
and external inputs: blinded mixtures, independently characterized structural
variants, a licensed/governed clinical bundle, a complete competitive nuclear
reference approved for the protocol, deployment identity/storage services, and
scientific/security/operations sign-off. These are planned deliverables, not
exceptions to the definition of done.

## 8. Work tracking format

Every roadmap item must have:

- an owner and reviewer;
- a linked specification or threat/validation requirement;
- positive, negative, malformed-input, and cancellation tests where applicable;
- a benchmark or stated reason why performance is irrelevant;
- documentation and migration impact;
- an evidence artifact attached to the release candidate.

“Code merged” is not an exit criterion. The unit of completion is verified,
reviewable evidence on the release candidate.
