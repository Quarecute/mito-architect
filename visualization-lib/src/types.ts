export type LayerName = 'genes' | 'coverage' | 'clusters' | 'svs' | 'snps';

/** Run-level metadata emitted by the native analysis core. */
export interface MitoMetadata {
  schema_version?: string;
  sv_event_schema_version?: string;
  engine_version: string;
  sample: string;
  input_path?: string;
  reference_path?: string;
  reference_accession?: string;
  reference_length: number;
  threads?: number;
  algorithm_notes?: string[];
  calling_parameters?: {
    min_mapping_quality: number;
    min_base_quality: number;
    excluded_snp_flags: number;
    numt_threshold?: number;
    development_tags_enabled?: boolean;
  };
  resources?: AnalysisResource[];
}

export interface AnalysisResource {
  name: string;
  version: string;
  path: string;
  sha256: string;
  source: string;
  license: string;
  retrieved: string;
}

/** NUMT and input-read accounting. */
export interface FilterStats {
  input_reads: number;
  passed_reads: number;
  numt_filtered_reads: number;
  numt_threshold: number;
  input_alignment_records?: number;
  input_molecules?: number;
  numt_assessment?: {
    mode: 'competitive_alignment' | 'mt_only_or_unknown' | 'unaligned_fastq' | string;
    nuclear_contigs_present: boolean;
    specificity_assessable: boolean;
  };
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
  /** @deprecated Use evidence_sources; retained for result-schema 0.4 readers. */
  evidence_source?: 'cigar' | 'split_alignment' | 'development_tag' | 'combined' | string;
  /** Sorted unique provenance categories merged under the canonical event ID. */
  evidence_sources?: Array<'cigar' | 'split_alignment' | 'development_tag' | string>;
  /** @deprecated Use orientations; retained for result-schema 0.4 readers. */
  orientation?: string;
  /** Sorted unique strand transitions observed for split-alignment evidence. */
  orientations?: string[];
  segment_count?: number;
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

/** Locus-level SNP evidence after MAPQ, base-quality, flag, and NUMT filters. */
export interface AggregateVariant extends SnpCall {
  alt_depth: number;
  ref_depth: number;
  other_depth: number;
  callable_depth: number;
  heteroplasmy: number;
  ci95_low: number;
  ci95_high: number;
  supporting_reads: string[];
}

/** Clinical annotation resolved from the local MITOMAP/ClinVar cache. */
export interface ClinicalAnnotation {
  phenotype?: string;
  pathogenicity?: string;
  references?: string[];
  source?: 'MITOMAP' | 'ClinVar' | 'local-cache' | string;
  sources?: string[];
  clinvar_allele_id?: string;
  mitomap_url?: string;
}

/** Summary statistics for the coverage layer. */
export interface CoverageMetrics {
  mean_depth: number;
  pct_sites_gt20x?: number;
  pct_bins_gt20x?: number;
  max_depth: number;
  mapping_quality_histogram: Array<{ mapq: number; count: number }>;
}

export interface ReadFeature {
  id: string;
  length: number;
  mean_quality: number;
  numt_score: number;
  filtered_numt: boolean;
  numt_evidence?: string[];
  cluster_id: number;
  mapping_quality?: number;
  flags?: number;
  reference_name?: string;
  aux_tags?: Record<string, string>;
  outlier?: boolean;
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
  haplogroup_assignment?: {
    resource: string;
    quality: number;
    contamination_warning: boolean;
    candidates: Array<{
      name: string;
      score: number;
      matched: string[];
      missing: string[];
      extra: string[];
    }>;
  };
}

/** Full analysis payload shared by core, CLI, server, report, and UI. */
export interface MitoAnalysisData {
  metadata: MitoMetadata;
  filter_stats: FilterStats;
  coverage_metrics?: CoverageMetrics;
  genes: GeneAnnotation[];
  coverage: CoverageBin[];
  svs: StructuralVariant[];
  variants?: AggregateVariant[];
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
