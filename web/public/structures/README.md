# Optional offline protein structures

This directory ships the curated RCSB BinaryCIF coordinate files and checksum
manifest used by the NGL viewer. It loads these local files first, then tries
the same immutable PDB entries at RCSB. Run
`bash scripts/fetch_protein_structures.sh` to verify and deliberately refresh
the bundled copies.

Curated mappings:

- Human MT-ATP6: PDB `8H9S`, chain `N`, residues 1-226.
- Human MT-ND4: PDB `9I4I`, chain `r`, residues 1-459.

These are experimental cryo-EM coordinates. The viewer interprets their
B-factor field as an experimental B-factor and does not contact a prediction
database.
