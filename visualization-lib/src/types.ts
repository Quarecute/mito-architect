export type LayerName = 'genes' | 'coverage' | 'clusters' | 'svs' | 'snps';

/** Run-level metadata emitted by the native analysis core. */
export interface MitoMetadata {
  schema_version?: string;
  engine_version: string;
  sample: string;
  input_path?: string;
  reference_path?: string;
  reference_accession?: string;
  reference_length: number;
  threads?: number;
  algorithm_notes?: string[];
}

/** NUMT and input-read accounting. */
export interface FilterStats {
  input_reads: number;
  passed_reads: number;
  numt_filtered_reads: number;
  numt_threshold: number;
}

/** rCRS gene interval used by circular tracks and gene search. */
export interface GeneAnnotation {
  name: string;
  start: number;
  end: number;
  strand: '+' | '-' | string;
  biotype: string;
}

/** Binned coverage depth over inclusive rCRS coordinates. */
export interface CoverageBin {
  start: number;
  end: number;
  depth: number;
}

/** Large deletion/insertion call with read-level support. */
export interface StructuralVariant {
  id: string;
  type: string;
  start: number;
  end: number;
  length: number;
  known_event: boolean;
  supporting_reads: string[];
  annotation?: ClinicalAnnotation;
}

/** Optional protein-structure mapping for non-synonymous mtDNA variants. */
export interface ProteinStructureMapping {
  structure_id?: string;
  chain?: string;
  residue_index?: number;
  residue_label?: string;
  complex?: string;
}

/** Single-nucleotide variant call and optional annotation payload. */
export interface SnpCall {
  position: number;
  ref: string;
  alt: string;
  gene?: string;
  consequence?: string;
  protein?: string;
  residue?: string;
  annotation?: ClinicalAnnotation;
  structure?: ProteinStructureMapping;
}

/** Clinical annotation resolved from the local MITOMAP/ClinVar cache. */
export interface ClinicalAnnotation {
  phenotype?: string;
  pathogenicity?: string;
  references?: string[];
  source?: 'MITOMAP' | 'ClinVar' | 'local-cache' | string;
  sources?: string[];
}

/** Summary statistics for the coverage layer. */
export interface CoverageMetrics {
  mean_depth: number;
  pct_sites_gt20x?: number;
  pct_bins_gt20x?: number;
  max_depth: number;
  mapping_quality_histogram: Array<{ range: string; count: number }>;
}

export interface ReadFeature {
  id: string;
  length: number;
  mean_quality: number;
  numt_score: number;
  filtered_numt: boolean;
  cluster_id: number;
  snps: SnpCall[];
  sv_ids: string[];
}

export interface ClusterSummary {
  id: number;
  label: string;
  haplogroup?: string;
  size: number;
  consensus_haplotype: string;
  reads: string[];
  sv_signature: Array<{ sv_id: string; support: number }>;
}

/** Full analysis payload shared by core, CLI, server, report, and UI. */
export interface MitoAnalysisData {
  metadata: MitoMetadata;
  filter_stats: FilterStats;
  coverage_metrics?: CoverageMetrics;
  genes: GeneAnnotation[];
  coverage: CoverageBin[];
  svs: StructuralVariant[];
  clusters: ClusterSummary[];
  reads: ReadFeature[];
}

export interface MitoCircosOptions {
  width?: number;
  height?: number;
  theme?: 'dark' | 'light';
  showControls?: boolean;
  activeLayers?: Partial<Record<LayerName, boolean>>;
  palette?: string[];
}

export interface LayerState {
  genes: boolean;
  coverage: boolean;
  clusters: boolean;
  svs: boolean;
  snps: boolean;
}
