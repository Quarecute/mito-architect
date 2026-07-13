# Documentation Index

Mito-Architect documentation is divided into user, scientific, engineering,
and release-evidence tracks. Status is gate-based: a roadmap milestone is not
complete until its documented exit conditions have evidence for the exact
candidate artifact.

## Start here

| Need | Document |
| --- | --- |
| Install and run the instrument | [`USER_GUIDE.md`](USER_GUIDE.md) |
| See what is implemented today | [`RELEASE_STATUS.md`](RELEASE_STATUS.md) |
| Understand the short milestone sequence | [`ROADMAP.md`](ROADMAP.md) |
| Execute the detailed path to 1.0 | [`PLAN_1.0.md`](PLAN_1.0.md) |
| Review mandatory production gates | [`PRODUCTION_ACCEPTANCE.md`](PRODUCTION_ACCEPTANCE.md) |
| Review scientific definitions | [`SCIENTIFIC_SPECIFICATION.md`](SCIENTIFIC_SPECIFICATION.md) |
| Integrate machine-readable failures | [`ERROR_CONTRACT.md`](ERROR_CONTRACT.md) |

## Evidence and history

- [`../CHANGELOG.md`](../CHANGELOG.md): repository-level change summary.
- [`VERIFICATION_0.5_RC1_DEV.md`](VERIFICATION_0.5_RC1_DEV.md): latest local
  development evidence and exact verification boundary.
- [`VERIFICATION_0.4.1.md`](VERIFICATION_0.4.1.md): 0.4.1 verification record.
- `RELEASE_NOTES_*.md`: changes and limitations for each software baseline.
- [`../fixtures/README.md`](../fixtures/README.md): deterministic fixtures,
  generated BAM/CRAM workflow, and bounded public-data workflow.

## Planning hierarchy

The documents intentionally have different roles:

1. `ROADMAP.md` is the compact RC1-to-1.0 milestone index.
2. `PLAN_1.0.md` defines engineering work, test matrices, and exit gates.
3. `RELEASE_STATUS.md` records current evidence and open closure conditions.
4. `PRODUCTION_ACCEPTANCE.md` is the release decision checklist.

`PRODUCTION_READINESS_PLAN.md` is retained for the historical 0.2-to-0.4
hardening programme. New work should be planned in `PLAN_1.0.md` and reflected
in `RELEASE_STATUS.md`.

## Evidence policy

- Development smoke tests do not constitute analytical validation.
- Evidence from an older binary does not close a newer release candidate.
- Scientific thresholds must be pre-registered before blinded evaluation.
- Clinical use requires the laboratory's quality system, accreditation, data
  governance, and jurisdiction-specific review in addition to software gates.
