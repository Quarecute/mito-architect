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
| `tiny.fastq` | Parser, read-quality, tagged development calls, and error handling | JSON smoke output |

`tiny.fastq` is not biological truth and its header tags are test controls. Raw
FASTQ has no coordinates, so it cannot exercise alignment-derived SNP, SV,
coverage, or mapping-quality behavior.

## Reproducible BAM and CRAM fixtures

Do not commit binary copies of the SAM fixture. Build coordinate-sorted BAM and
reference-backed CRAM forms with the installed samtools version:

```bash
bash scripts/build_bam_fixture.sh
cargo run -p mito-cli -- analyze -i fixtures/generated/truth_mixed.bam --json
cargo run -p mito-cli -- analyze -i fixtures/generated/truth_mixed.cram --json
```

The generated BAM/BAI and CRAM/CRAI are ignored by Git. Use
`bash scripts/verify_determinism.sh` to compare canonical output at 1, 2, 4, 8,
and 16 requested threads.

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
