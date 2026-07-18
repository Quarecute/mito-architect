# Scientific Product Positioning

Status: normative product-direction document for the path from `0.5.0-dev` to
research version 1.0. It defines the scientific problem Mito-Architect is built
to solve; it does not claim that the planned methods are already validated.

## 1. Product identity

Mito-Architect is a long-read platform for reconstructing the molecular
architecture of human mitochondrial DNA. It combines standard sample-level
variant analysis with physical molecule-level linkage, structural-event paths,
and uncertainty-aware reconstruction of coexisting mtDNA architectures.

The product must answer three connected questions from one evidence model:

1. **Variant-centric:** which normalized events are present, at what frequency,
   and with what callable depth and quality evidence?
2. **Molecule-centric:** which events are physically observed together on the
   same independently defined molecule, and which positions were actually
   callable?
3. **Architecture-centric:** which reproducible molecule types coexist in the
   sample, what defines them, how abundant are they, and how stable is the
   inference?

The core product statement is:

> Variant-centric analysis describes the composition of changes;
> molecule-centric analysis describes their physical organization;
> architecture-centric analysis describes the reproducible structure of the
> molecular population.

Mito-Architect is not positioned as another generic SNV/indel web report and
must not be described as the first tool to phase or sequence individual mtDNA
molecules.

## 2. Defensible novelty

The intended contribution is the integrated scientific model, not any one
commodity feature:

> Mito-Architect builds linked variant-centric and molecule-centric projections
> from a common callable-aware evidence graph, jointly phases small variants
> and structural events, canonicalizes multi-junction circular mtDNA paths, and
> reconstructs coexisting molecular architectures with explicit uncertainty and
> bidirectional traceability to supporting alignments.

The strongest claims, once independently validated, are:

- callable-aware joint phasing of SNV, small indel, SV, and complex-event paths;
- explicit distinction between reference, alternate, event-absent,
  not-callable, low-quality, and conflicting observations;
- strand-invariant and rotation-aware canonicalization of circular
  multi-junction paths;
- architecture signatures, abundance confidence intervals, ambiguous
  assignments, and stability under resampling;
- bidirectional `variant <-> molecule <-> architecture <-> alignment`
  traceability;
- one deterministic result contract shared by CLI, service, reports, and UI.

These are hypotheses and engineering objectives until the corresponding
independent mixture, SV, NUMT, and resampling gates pass.

## 3. Relationship to existing tools

### mtDNA-Server 2

mtDNA-Server 2 is the primary comparator for standard variant-centric
capability. Its published workflow combines input validation, FastQC/MultiQC,
mutserve2 SNV calling, Mutect2 indel calling, fusion mode, coverage estimation,
HaploGrep 3, haplocheck, annotation, Nextflow execution, and interactive
reporting. Mito-Architect does not claim novelty for reproducing these features.

The intended distinction is the primary analytical unit:

| Question | mtDNA-Server 2 published focus | Mito-Architect 1.0 focus |
| --- | --- | --- |
| Which variants are present? | primary result | required parity projection |
| What is their HF? | primary result | derived from molecule-callable evidence |
| Which events are on the same molecule? | not the primary result model | first-class result |
| Are SNP and SV physically linked? | not the primary result model | first-class phasing edge |
| Which molecule architectures coexist? | not the primary result model | primary scientific output |
| Can an aggregate be traced to an alignment? | aggregate/sample reports | bidirectional evidence graph |

The mtDNA-Server 2 article states that its current output annotates detected
NUMT positions and that NUMT-derived false positives are not removed by the
published low-coverage model; predefined-NUMT realignment is described as
future work. Mito-Architect therefore keeps NUMT specificity non-assessable
unless competitive nuclear-plus-mitochondrial evidence is present. This is a
different contract, not yet a validated superiority claim.

### baldur and full-length native mtDNA sequencing

The published baldur workflow already detects heteroplasmy, physically phases
selected SNVs, and characterizes complex deletions in full-length native mtDNA
reads. It proves that molecule-level mtDNA analysis is scientifically valuable
and also prevents Mito-Architect from claiming that molecule phasing itself is
new.

Mito-Architect must differentiate itself through callable-aware joint SNV/SV
phasing, generalized event traceability, circular multi-junction
canonicalization, and validated population-level architecture inference rather
than through the existence of a per-read matrix alone.

### iMiGseq

iMiGseq provides full-length mtDNA analysis and complete haplotyping at the
individual-molecule level, including genetic linkage and large deletions. It is
another reason to avoid an unsupported "first molecule-centric mtDNA tool"
claim.

### Himito and future comparators

Himito is a 2025 bioRxiv preprint describing a graph-based long-read mtDNA
toolkit with NUMT filtering, haplotype assembly, variant calling, and methylation
analysis. Because it is a preprint and the field is moving, its methods,
version, availability, and overlap with Mito-Architect must be re-evaluated
when the benchmark protocol is frozen. The comparator registry is a versioned
release resource, not a static marketing table.

## 4. Common evidence model

The target analysis graph is:

```text
SAM/BAM/CRAM records
        |
        v
alignment fragments ---- provenance / quality / flags / query interval
        |
        v
molecule assembler ----- protocol-aware identity and exclusion decisions
        |
        v
normalized observations + normalized small/SV/complex events
        |
        v
callable-aware sparse molecule-event evidence store
        |
        +------------------+------------------+
        v                  v                  v
variant projection   molecule projection   phase links
        \                  |                  /
         +-----------------+-----------------+
                           v
                candidate architectures
                           |
                           v
              one sample result and report
```

Variant and molecule results must never be produced by independent callers that
can silently disagree. Every aggregate count must be derivable from the same
normalized observation set and resolvable to molecule and alignment evidence.

## 5. Required scientific entities

The 1.0 contract must define at least:

- `AlignmentFragment`: one parsed primary, secondary, supplementary, duplex, or
  paired evidence fragment with immutable source identity;
- `Molecule`: the protocol-aware unit assembled from one or more fragments;
- `CallableInterval`: an interval and reasoned state for what can be assessed;
- `Event`: a normalized SNV, indel, SV edge, or complex circular path;
- `Observation`: one molecule-event state with evidence and provenance;
- `PhaseLink`: joint-callability and co-occurrence counts for an event pair;
- `Architecture`: a stable, explicitly uncertain group of molecular profiles;
- `Interpretation`: haplogroup, clinical, protein, or other resource-backed
  projection that never rewrites the underlying observation.

Required observation states are semantically distinct:

```text
REFERENCE       locus callable; reference allele supported
ALTERNATE       locus callable; normalized alternate supported
EVENT_ABSENT    structural event assessable and not observed
NOT_CALLABLE    molecule does not make the event assessable
LOW_QUALITY     evidence exists but fails the declared quality contract
CONFLICT        incompatible evidence exists within the assembled molecule
```

`NOT_CALLABLE`, `LOW_QUALITY`, and `CONFLICT` must never be coerced to
`REFERENCE` or `EVENT_ABSENT`.

## 6. Product boundaries

Version 1.0 is a production-grade **research instrument** available as a local
CLI/container and a single-laboratory web application. It is not a clinical
diagnostic device and does not infer cellular clones from bulk sequencing.

The permitted language before validation is:

- molecular profile;
- candidate molecular architecture;
- candidate haplotype;
- molecule group;
- ambiguous or unassigned molecule.

The following claims require cell-level evidence and are not inferred from bulk
long reads:

- cell clone or subclone;
- mitochondrial lineage within a specific cell population;
- tissue-specific causal pathogenicity;
- diagnostic verdict.

A hardened multi-tenant hosted service is an additional deployment profile. It
requires durable storage, OIDC/RBAC, tenant isolation, audit, backup/restore,
and security review, but it is not allowed to displace the scientific critical
path to research 1.0.

## 7. Validation principle

Standard variant performance must be compared with mtDNA-Server 2 and other
appropriate callers on the same frozen inputs. Mito-Architect does not need to
win every standard metric, but must meet pre-registered usability thresholds.

Its additional claims require known molecular mixtures and metrics unavailable
to a purely variant-centric benchmark:

- joint-callable phase precision/recall and false-linkage rate;
- molecule assignment accuracy and phase-switch error;
- architecture precision/recall, adjusted Rand index, and normalized mutual
  information;
- architecture abundance error and confidence-interval coverage;
- rare-architecture detection limit;
- stability under molecule and coverage subsampling;
- circular breakpoint/path equivalence and false path-merging rate;
- NUMT leakage and true-mtDNA loss under competitive mapping.

Thresholds and exclusions must be fixed before blinded truth is inspected.
Development fixtures prove software behavior only.

## 8. Primary literature used for positioning

- Weissensteiner et al. (2024), mtDNA-Server 2:
  <https://doi.org/10.1093/nar/gkae296>
- Keraite et al. (2022), full-length native mtDNA and baldur:
  <https://doi.org/10.1038/s41467-022-33530-3>
- Bi et al. (2023), iMiGseq:
  <https://doi.org/10.1093/nar/gkad208>
- Su et al. (2025), Himito preprint:
  <https://doi.org/10.1101/2025.11.03.686348>

Source claims must be rechecked against the cited version when a paper,
benchmark, or release claim is prepared.
