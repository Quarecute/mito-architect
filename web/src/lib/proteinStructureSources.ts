export type StructureFileExtension = 'pdb' | 'cif' | 'bcif';

export interface StructureSourceCandidate {
  url: string;
  ext: StructureFileExtension;
  label: string;
  entryUrl: string;
  chain: string;
  coordinateDescription: string;
  confidenceSemantics: 'bfactor';
}

export function normalizePdbId(structureId: string): string | undefined {
  const trimmed = structureId.trim();
  return /^[0-9][A-Z0-9]{3}$/i.test(trimmed) ? trimmed.toUpperCase() : undefined;
}

export function rcsbCandidates(pdbId: string, chain: string): StructureSourceCandidate[] {
  const normalized = normalizePdbId(pdbId);
  if (!normalized) {
    return [];
  }
  const lowerId = normalized.toLowerCase();
  const common = {
    ext: 'bcif' as const,
    entryUrl: `https://www.rcsb.org/structure/${normalized}`,
    chain,
    coordinateDescription: 'experimental cryo-EM coordinates',
    confidenceSemantics: 'bfactor' as const
  };
  return [
    {
      ...common,
      url: `/structures/${lowerId}.bcif`,
      label: `RCSB PDB ${normalized} (local cache)`
    },
    {
      ...common,
      url: `https://models.rcsb.org/${lowerId}.bcif`,
      label: `RCSB PDB ${normalized}`
    },
    {
      ...common,
      url: `https://files.rcsb.org/download/${normalized}.cif`,
      ext: 'cif',
      label: `RCSB PDB ${normalized} (mmCIF fallback)`
    }
  ];
}
