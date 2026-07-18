# RC1 Scientific Rule Traceability

Status: live development ledger for result schema 0.5. `covered` means the
software rule has executable local evidence; it does not mean analytical or
clinical validation. RC1 still requires two independent reviewers.

| Rule / contract | Positive evidence | Negative / boundary evidence | State |
| --- | --- | --- | --- |
| Strict FASTQ/SAM/CIGAR ingest | native smoke; `truth_*.sam` | malformed FASTQ and CIGAR cases in `fixtures/negative/error_manifest.json` | covered |
| Fail-closed `SA:Z` grammar | split/complex SV goldens | field-count, position, strand, CIGAR, MAPQ, NM, terminator cases | covered |
| QNAME molecule collapse | `truth_split.sam` primary + supplementary records | supplementary/default excluded flags | covered; mapper differential pending |
| Circular rCRS SNP coordinates | `truth_snp_edges.sam` position 16,569/1 span | invalid reference/input coordinates fail | covered |
| Callable depth and multi-allelic AD/HF | `truth_snp_edges.expected.evidence.json` | missing quality, excluded flags, N/non-callable paths | covered |
| Repeated-locus conflict | `truth_molecule_edges.sam` | conflicting repeat contributes neither allele nor support | covered |
| Wilson HF interval | SNP/molecule/QC evidence goldens | zero/non-callable observations excluded | covered |
| Molecule/strand support | all three evidence goldens | absent comparison emits `null`, no inferred strand | covered |
| Per-allele Phred summary | `truth_qc_metrics.expected.evidence.json` | empty allele groups emit explicit nulls | covered |
| Read-position representative selection | repeated same-allele Q20/Q40 case in `truth_qc_metrics.sam` | conflicting repeat is non-callable | covered |
| Observational bias semantics | exact strand/read-position deltas in QC golden | no hard filter or classification mutation exists | covered; RC5 calibration required |
| Competitive NUMT assessment | `truth_split.sam`, native smoke | mt-only/FASTQ explicitly non-assessable | development baseline; blinded leakage truth pending |
| Canonical SV edge schema 1.0 | `truth_split.expected.svs.json` | coordinate collision is internal failure | covered; independent biological truth pending |
| Complex path schema 1.0 | `truth_complex.expected.svs.json` | partial/circular path interpretation remains unsupported | software equivalence covered; truth pending |
| PhyloTree grammar and ordering | 28-profile/23-lineage pinned HaploGrep differential | malformed tree tokens fail `MITO-E1302` | covered representative subset |
| PhyloTree backmutation | H2c1/H1b1g/V7b remove `152C`/`1438G`/`16298C`, exact HaploGrep top three | validator rejects a removed marker reappearing in winner matched/missing | covered three explicit cases |
| Official alignment-rule normalization | two CIGAR 3-prime cases plus three compound marker-set projections | malformed/truncated/unsupported/non-convergent resource fails closed; reference restoration emits no mutation | covered software differential; mapper review pending |
| Clinical assertion preservation | `clinical/conflicting_assertions.expected.json` | no field-wise overwrite | covered schema 1.0 |
| Clinical conflicts/resource validation | synthetic pathogenic-vs-benign golden | four-case clinical negative manifest | covered software contract; governed bundle pending |
| Clinical override provenance | result resource entry becomes external-unpinned/not-verified | override never inherits bundled hash | covered boundary; signed bundles pending |
| Clinical snapshot lifecycle | stage/activate/status/two-way rollback golden | unsafe ID and post-stage corruption fail; rollback validates before state write | covered mechanism; licensed content pending |
| Result/error schema boundary | CLI/server validators and Rust tests | missing schema metadata rejected | covered |
| HTTP route contract | in-process real Router goldens for health/readiness/status/result/download | auth, malformed UUID/multipart, unknown/pending/completed states | covered development baseline |
| Thread determinism | `scripts/verify_determinism.sh`, 1/2/4/8/16 | canonical digest mismatch fails | covered on current host only |
| Immutable candidate evidence | per-file/tool/commit/tree/log index generator | dirty tree or missing release log refused | mechanism covered; exact clean RC1 run pending |

## Open RC1 closure items

1. Add independent mapper-level/biological review for the committed official
   rules and backmutation corpus.
2. Add independently reviewed complex/circular rearrangement truth rather than
   only software-representation equivalence.
3. Assign two scientific reviewers and record dated disposition of every row.
4. Acquire/freeze a governed, licensed clinical snapshot and validate match
   rate and conflict retention through the implemented lifecycle.
5. Rerun all evidence on one clean immutable RC1 artifact and generate the
   implemented evidence index with its exact verification-log hash.
