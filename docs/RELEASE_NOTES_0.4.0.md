# Release Notes 0.4.0

Version 0.4.0 completes the planned software implementation for the 0.3.x and
0.4 milestones. It remains research-use only: software verification is not a
substitute for blinded analytical validation or an operated production service.

## Analysis contract 0.4

- Collapses primary and supplementary SAM/BAM records by molecule name while
  retaining raw-alignment and molecule counts.
- Reconstructs deletion, insertion, duplication, inversion, and circular-origin
  events from primary plus `SA` split-alignment segments and records evidence
  source, orientation, and segment count.
- Distinguishes competitive nuclear-plus-mitochondrial NUMT assessment from
  mtDNA-only or unaligned inputs and records per-molecule evidence.
- Adds resource-backed haplogroup candidates using the complete versioned
  PhyloTree rCRS 17.3 XML and weights, with matched, missing, extra, quality,
  and contamination-warning fields.
- Emits a checksum-bearing resource manifest in every result.

## Clinical and structure provenance

- Preserves ClinVar allele IDs and MITOMAP links in read-level and aggregate
  variant annotations.
- Bundles the curated RCSB BinaryCIF models used by MT-ATP6 and MT-ND4 mappings,
  verifies their checksums, and loads them locally before network fallbacks.
- Restricts the NGL viewer to curated four-character PDB identifiers. It does
  not query a protein-prediction database.
- Adds a mapped-chain default view, optional full-complex view, 5 Å residue
  context, B-factor/secondary/uniform coloring, provenance, and PNG export.

## Reproducibility and tests

- Vendors PhyloTree 17.3 XML, weights, and license and verifies all mandatory
  resource checksums offline.
- Adds deterministic tests for structure source resolution, competitive NUMT
  evidence, molecule collapsing, split-event classes, and haplogroup ranking.
- Adds scripts for reproducible SAM-to-BAM fixtures, a bounded public ONT data
  bundle, and construction of a nuclear-plus-rCRS minimap2 index.

## Known scientific limitations

- The current PhyloTree scorer handles simple SNVs and backmutations; complex
  indel nomenclature and HaploGrep-equivalent quality calibration remain 1.0
  validation work.
- The bundled clinical table is a curated development subset, not a complete
  licensed clinical knowledgebase.
- Low-frequency sensitivity, SV accuracy, NUMT leakage, haplogroup accuracy,
  and cross-run reproducibility have not yet passed the blinded validation
  programme in `PRODUCTION_ACCEPTANCE.md`.
- Job storage, user isolation, audit trails, recovery, and production
  observability remain future work.
