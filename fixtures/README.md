# Validation Data

The repository separates deterministic truth fixtures from public biological
data. Small truth files are committed so CI remains fast and offline. Public
human reads are downloaded only on demand and stay outside version control.

## Committed fixtures

| Input | Purpose | Expected output |
| --- | --- | --- |
| `truth_snp.sam` | Pathogenic m.3243A>G, MAPQ, aux tags, and NUMT accounting | `truth_snp.expected.vcf` |
| `truth_mixed.sam` | Reference mismatch, deletion, insertion, right soft clip, unmapped read | `truth_mixed.expected.vcf` |
| `truth_split.sam` | Canonical CIGAR/split/reverse-strand SV equivalence, circular origin, and competitive NUMT evidence | `truth_split.expected.svs.json` |
| `truth_complex.sam` | Forward/reverse three-segment molecule equivalence and ordered multi-junction coalescing | `truth_complex.expected.svs.json` |
| `truth_snp_edges.sam` | Circular 16,569/1 SNP and multi-allelic callable-depth semantics | `truth_snp_edges.expected.evidence.json` |
| `truth_molecule_edges.sam` | Reverse-strand calls, missing quality, default flag exclusions, and repeated-locus molecule conflict | `truth_molecule_edges.expected.evidence.json` |
| `truth_qc_metrics.sam` | Exact strand, Phred, read-position, multi-allelic, and repeated-same-allele representative selection | `truth_qc_metrics.expected.evidence.json` |
| `truth_evidence_graph.sam` | Schema 0.6 fragment retention, supplementary membership, missing/multiple-primary ambiguity, sparse SNV states, and aggregate parity | `truth_evidence_graph.expected.json` |
| `truth_phase_evidence.sam` | Callable-aware SNV/SNV and SNV/small-indel phase partitions, alternative insertion alleles, non-joinable split-fragment adjacency, explicit indel absence, and pairwise confidence metrics | `truth_phase_evidence.expected.json` |
| `truth_phase_matrix.sam` | Exact SNV/small-indel/SV-deletion/SV-insertion pair classes, sparse NOT_CALLABLE exclusion, uncertain co-occurrence, supporting/uncertain molecule IDs, and null-statistic round-trip | `truth_phase_matrix.expected.json` |
| `truth_circular_indel.sam` | Origin-adjacent small deletion, callable reference support, circular normalization, and explicit non-representability in linear VCF | `truth_circular_indel.expected.json` |
| `truth_protocol_identity.sam` | Explicit `MI` grouping, UMI/duplex provenance, multi-QNAME assembly, missing identity, conflicting UMI, duplicate-only exclusion, and partial-tag inheritance | `truth_protocol_identity.expected.json` |
| `comparator/mtdna-server-2-v2.1.16.truth-phase.tsv` | Pinned `variants.annotated.txt` callset/HF differential, plus non-PASS policy behavior | comparator adapter JSON assertions |
| `comparator/mtdna-server-2-v2.1.16.invalid-coverage.tsv` | Invalid forward/reverse/total coverage accounting must fail closed | process-level negative gate |
| `comparator/mtdna-server-2-v2.1.16.invalid-level.tsv` | Non-finite/out-of-range allele fraction must fail closed | process-level negative gate |
| `comparator/mtdna-server-2-v2.1.16.multi-sample.tsv` | Ambiguous multi-sample selection must require `--sample` | process-level negative gate |
| `comparator/mtdna-server-2-v2.1.16.empty.tsv` | Header-only no-call comparator requires an explicit sample and produces a valid empty callset | comparator adapter JSON assertions |
| `clinical/conflicting_assertions.tsv` | Synthetic incompatible clinical sources with complete assertion provenance | `clinical/conflicting_assertions.expected.json` |
| `haplogroup/haplogrep3-3.3.2.hsd` | Exact 28-profile input run through official HaploGrep 3.3.2 with PhyloTree rCRS 17.3 | `haplogroup/haplogrep3-3.3.2.expected.json` |
| `haplogroup/truth_haplogroup_indels.sam` | Two repeat-equivalent CIGAR representations normalized to `523d`, `524d`, and `573.1CCCC`, including callable range | same haplogroup manifest |
| `haplogroup/truth_haplogroup_rotated_indel.sam` | Repeat-motif insertion normalized from `8270.1ACCCCCTCT` to `8289.1CCCCCTCTA` | same haplogroup manifest |
| `haplogroup/truth_haplogroup_rule_309.fastq` | Official compound rule `309T 310C -> 309d` | same haplogroup manifest |
| `haplogroup/truth_haplogroup_rule_188.fastq` | Official mixed insertion/SNP/deletion rewrite around 188-196 | same haplogroup manifest |
| `haplogroup/truth_haplogroup_rule_56.fastq` | Official compound rewrite with explicit reference restoration at 65 | same haplogroup manifest |
| `haplogroup/invalid_alignment_rules.csv` | Truncated official-rule resource must fail closed with `MITO-E1302` | process-level negative gate |
| `tiny.fastq` | Strict four-record FASTQ parser/read-quality and documented CLI smoke | valid schema-0.5 JSON |

`tiny.fastq` is not biological truth and its header tags are test controls. Raw
FASTQ has no coordinates, so it cannot exercise alignment-derived SNP, SV,
coverage, or mapping-quality behavior.

## Negative input corpus

`negative/error_manifest.json` is the versioned error golden. Its 11 cases pin
`MITO-E1103` for malformed FASTQ separators, primary CIGAR/query disagreement,
and invalid `SA:Z` field count, position, strand, CIGAR, MAPQ, edit distance,
termination, type, and empty value. Run it with:

```bash
cargo run -p mito-cli -- validate-error-manifest \
  --manifest fixtures/negative/error_manifest.json
```

Diagnostics remain human-facing and may become more specific; the manifest
compares the stable error-schema version and code, never message text.

## Haplogroup differential corpus

`haplogroup/haplogrep3-3.3.2.expected.json` pins the official Linux release
asset checksum, exact PhyloTree 17.3 tree, weights, and alignment-rule checksums,
Kulczynski metric, tested ranges, top three reference hits, and reference
quality. The default validator gates exact marker/range projection and ordered
top-three agreement for 28 profiles covering 23 declared lineage groups and
three explicit inherited-marker removals. Five additional projection cases
prove small-indel extraction, repeat-equivalent 3-prime normalization,
inserted-motif rotation from CIGAR/SEQ/QUAL, and three compound official-rule
rewrites. The compound FASTQ headers are explicit development controls, not
biological evidence.

Run the offline committed gate with:

```bash
cargo run -p mito-cli -- validate-haplogroup-manifest \
  --manifest fixtures/haplogroup/haplogrep3-3.3.2.expected.json
```

The corpus is representative differential evidence, not a population
haplogroup-accuracy validation set. RC1 still requires independent scientific
review and mapper-level biological differential evidence; RC2 requires
independent accuracy and quality calibration.

The SNP evidence goldens also pin named molecule counts, forward/reverse
support, per-allele Phred summaries, and normalized read-position/strand bias
deltas. These are reported QC observations, not validated filters. The
clinical golden intentionally uses `TEST-A`/`TEST-B` synthetic sources so no
fixture assertion can be mistaken for patient-facing medical knowledge.

The schema 0.6 third-slice goldens are exercised through the complete C++ -> C
ABI -> Rust CLI boundary. They cover callable-aware event evidence, columnar
pages, event-specific negative-evidence rules, and explicit protocol identity:

```bash
cargo run -p mito-cli -- validate-evidence-graph-fixture \
  --input fixtures/truth_evidence_graph.sam \
  --expected-json fixtures/truth_evidence_graph.expected.json
cargo run -p mito-cli -- validate-evidence-graph-fixture \
  --input fixtures/truth_protocol_identity.sam \
  --expected-json fixtures/truth_protocol_identity.expected.json \
  --molecule-id-tag MI --umi-tag RX --duplex-tag DX
bash scripts/verify_phase_matrix.sh
bash scripts/verify_comparator_adapter.sh
bash scripts/verify_repeat_rotation.sh
```

These are software contract and comparator-format fixtures. They are not
biological validation of molecule identity, low-frequency detection, phase
accuracy, false-linkage rate, or comparator superiority.

## Reproducible BAM and CRAM fixtures

Do not commit binary copies of the SAM fixture. Build coordinate-sorted BAM and
reference-backed CRAM forms with the installed samtools version:

```bash
bash scripts/build_bam_fixture.sh
cargo run -p mito-cli -- analyze -i fixtures/generated/truth_mixed.bam --json
cargo run -p mito-cli -- analyze -i fixtures/generated/truth_mixed.cram --json
```

The generated BAM/BAI and CRAM/CRAI are ignored by Git. The build script also
creates BAM/CRAM pairs for `truth_phase_evidence`, `truth_phase_matrix`, and
`truth_protocol_identity` and `truth_repeat_rotation`, plus a BAM form of
`truth_circular_indel`. `truth_repeat_rotation` proves that left/right repeat
representations and cyclically rotated insertion motifs collapse to the same
three normalized events across record order, requested threads, SAM, BAM, and
CRAM. The
origin-crossing CIGAR is intentionally not converted to CRAM: a CRAM slice is
linear and cannot losslessly encode reference consumption beyond position
16569. CRAM producers must rotate or split such an alignment before encoding.
The canonical gate requires every representable form to reproduce the same
normalized schema-0.6 golden, including arbitrary `MI`/`RX`/`DX` aux-tag
preservation. Use
`bash scripts/verify_determinism.sh` to compare canonical output at 1, 2, 4, 8,
and 16 requested threads.

The RC3 export gate verifies unified JSON/TSV/VCF counts, bcftools parsing,
bgzip/tabix round-trip, provenance checksums, deterministic reruns, and the
fail-closed circular-origin VCF boundary:

```bash
bash scripts/verify_rc3_exports.sh
```

## Public Oxford Nanopore long-read data

`SRR18110025` is the public Oxford Nanopore GM12878 run from BioProject
`PRJNA809571`. It contains targeted, amplification-free long mtDNA reads used in
the nCATS-mtDNA study. The complete run is about 2.22 Gbp, so the acquisition
script defaults to the first 5,000 spots and records checksums and tool versions.

Install `sra-tools`, `minimap2`, and `samtools`, then run:

```bash
bash scripts/fetch_public_fixture.sh
```

Override the bounded sample size or worker count as needed:

```bash
MITO_PUBLIC_MAX_SPOTS=20000 MITO_PUBLIC_THREADS=8 \
  bash scripts/fetch_public_fixture.sh
```

Outputs are written under `.data/public/SRR18110025/`: raw FASTQ, rCRS-aligned
SAM, coordinate-sorted BAM, BAI, and a provenance manifest. Public human data
still require appropriate institutional handling even when access is open.

This public run is useful for realistic read lengths, mapping, phasing, and
coverage behavior. It is not a calibrated low-frequency or deletion truth set;
analytical sensitivity and specificity must be validated separately with
mixtures or reference materials whose molecule-level truth is known.

For NUMT assessment, provide a prebuilt combined human nuclear plus rCRS
reference instead of mapping only to rCRS:

```bash
MITO_ALIGNMENT_REFERENCE=/resources/hg38-plus-rcrs.mmi \
  bash scripts/fetch_public_fixture.sh
```

Competitive mapping is required to expose primary nuclear alignments and `SA`
supplementary alignments to the core NUMT filter. An mtDNA-only alignment cannot
establish NUMT specificity.

Build a versioned combined index without duplicating the mitochondrial contig:

```bash
MITO_NUCLEAR_FASTA=/resources/GRCh38.primary_assembly.genome.fa \
  bash scripts/build_competitive_reference.sh
```

The script records input, combined-reference, and minimap2-index checksums plus
the minimap2 version under `.data/reference/reference-manifest.txt`.
