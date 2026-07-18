export type LayerName = 'genes' | 'coverage' | 'clusters' | 'svs' | 'snps';

/** Run-level metadata emitted by the native analysis core. */
export interface MitoMetadata {
  schema_version?: string;
  sv_event_schema_version?: string;
  complex_sv_event_schema_version?: string;
  clinical_annotation_schema_version?: string;
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
    max_evidence_observations?: number;
    max_phase_links?: number;
    evidence_page_size?: number;
    molecule_id_tag?: string;
    umi_tag?: string;
    duplex_tag?: string;
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
  evidence_eligible_molecules?: number;
  ambiguous_molecules?: number;
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
  /** Schema 0.6 unified evidence-graph reference. */
  event_id?: string;
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

/** Strand-invariant ordered path of two or more junctions observed on one molecule. */
export interface ComplexStructuralEvent {
  id: string;
  /** Schema 0.6 unified evidence-graph reference. */
  event_id?: string;
  junction_count: number;
  segment_count: number;
  canonicalization: 'strand_invariant_path' | string;
  /** Canonical SV edge IDs in molecule traversal order, normalized against path reversal. */
  junction_ids: string[];
  /** Strand transition paired by index with each canonical junction ID. */
  junction_orientations: string[];
  supporting_reads: string[];
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

/** Unified schema 0.6 SNV/small-indel projection derived from the evidence graph. */
export interface AggregateVariant extends SnpCall {
  /** Schema 0.6 normalized event reference. */
  event_id?: string;
  type?: 'SNV' | 'SMALL_INSERTION' | 'SMALL_DELETION' | string;
  start?: number;
  end?: number;
  length?: number;
  normalization?: string;
  negative_evidence_rule?: EvidenceEvent['negative_evidence_rule'];
  assessability?: EvidenceEvent['assessability'];
  vcf_position?: number;
  vcf_representable?: boolean;
  alt_depth: number;
  ref_depth: number;
  other_depth: number;
  event_absent_depth?: number;
  low_quality_depth?: number;
  conflict_depth?: number;
  callable_depth: number;
  heteroplasmy: number;
  ci95_low: number;
  ci95_high: number;
  supporting_reads: string[];
  supporting_molecule_ids?: string[];
  filter_status?: 'NOT_CALIBRATED' | string;
  qc_flags?: string[];
  numt_assessability?: 'ASSESSABLE' | 'NOT_ASSESSABLE' | string;
  multi_allelic?: boolean;
  homopolymer_context?: {
    reference_base: string | null;
    run_length: number;
  };
  molecule_support?: {
    alternate: number;
    reference: number;
    other: number;
    callable: number;
  };
  strand_support?: {
    alt_forward: number;
    alt_reverse: number;
    ref_forward: number;
    ref_reverse: number;
    other_forward: number;
    other_reverse: number;
  };
  /** Absolute alternate-vs-reference forward-strand fraction difference; observational only. */
  strand_bias_delta?: number | null;
  allele_quality?: {
    alternate: AlleleQualitySummary;
    reference: AlleleQualitySummary;
    other: AlleleQualitySummary;
  };
  mapping_quality?: {
    alternate: NumericQualitySummary;
    reference: NumericQualitySummary;
    other: NumericQualitySummary;
  };
  read_position?: {
    definition: 'normalized_center_proximity' | string;
    alternate_mean: number | null;
    reference_mean: number | null;
    other_mean: number | null;
    bias_delta: number | null;
  };
}

export interface NumericQualitySummary {
  count: number;
  mean: number | null;
  min: number | null;
  max: number | null;
}

export interface AlleleQualitySummary {
  count: number;
  mean_phred: number | null;
  min_phred: number | null;
  max_phred: number | null;
}

/** One source-preserving clinical record; empty optional source fields are not inferred. */
export interface ClinicalAssertion {
  source: string;
  assertion_id: string;
  allele_id: string;
  disease: string;
  clinical_significance: string;
  normalized_significance: string;
  review_status: string;
  assertion_date: string;
  source_url: string;
  references: string[];
  resource_version: string;
  retrieved_at: string;
}

/** Deterministic summary plus lossless source assertions from the local clinical cache. */
export interface ClinicalAnnotation {
  schema_version?: string;
  phenotype?: string;
  /** @deprecated Prefer consensus_significance and inspect assertions/conflict_status. */
  pathogenicity?: string;
  consensus_significance?: string;
  conflict_status?: 'single_assertion' | 'consistent' | 'conflicting' | string;
  assertions?: ClinicalAssertion[];
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
  haplogroup_markers?: string[];
  haplogroup_range_known?: boolean;
  haplogroup_callable_ranges?: Array<{ start: number; end: number }>;
  sv_ids: string[];
  complex_event_ids?: string[];
}

export interface ClusterSummary {
  id: number;
  label: string;
  haplogroup?: string;
  size: number;
  consensus_haplotype: string;
  reads: string[];
  sv_signature: Array<{ sv_id: string; support: number }>;
  complex_event_signature?: Array<{ event_id: string; support: number }>;
  haplogroup_assignment?: {
    resource: string;
    quality: number;
    contamination_warning: boolean;
    observed_markers?: string[];
    callable_ranges?: Array<{ start: number; end: number }>;
    candidates: Array<{
      name: string;
      score: number;
      matched: string[];
      missing: string[];
      extra: string[];
    }>;
  };
}

export type ObservationState =
  | 'REFERENCE'
  | 'ALTERNATE'
  | 'EVENT_ABSENT'
  | 'NOT_CALLABLE'
  | 'LOW_QUALITY'
  | 'CONFLICT';

/** Original SAM/BAM/CRAM record retained before molecule assembly. */
export interface AlignmentFragment {
  id: string;
  source_record_index: number;
  molecule_id: string;
  molecule_index: number;
  role:
    | 'primary_candidate'
    | 'secondary'
    | 'supplementary'
    | 'secondary_and_supplementary'
    | 'unaligned_read'
    | string;
  selected_representative: boolean;
  flags: number;
  strand: '+' | '-' | string;
  mapping_quality: number;
  reference_name: string;
  reference_start: number;
  cigar: string;
  query_length: number;
  base_qualities_available: boolean;
  aux_tags: Record<string, string>;
}

/** Deterministic result of the configured fragment-to-molecule policy. */
export interface AssembledMolecule {
  id: string;
  index: number;
  identity_policy: 'sam_qname' | 'fastq_record_proxy' | string;
  assembly_status:
    | 'unique_primary'
    | 'fallback_without_primary'
    | 'first_of_multiple_primaries'
    | 'single_unaligned_fragment'
    | string;
  primary_candidate_count: number;
  ambiguous: boolean;
  analysis_eligible: boolean;
  evidence_eligible: boolean;
  callability_status: string;
  callable_bases: number;
  callable_fraction: number;
  query_length: number;
  mean_base_quality: number;
  mapping_quality: number;
  numt_score: number;
  numt_evidence: string[];
  cluster_id: number;
  alternate_event_ids: string[];
  evidence_state_counts: {
    alternate: number;
    reference: number;
    event_absent: number;
    low_quality: number;
    conflict: number;
  };
  representative_alignment_id: string;
  source_qnames: string[];
  protocol_metadata: Record<string, string>;
  protocol_flags: string[];
  exclusion_reasons: string[];
  alignment_ids: string[];
  warnings: string[];
}

export interface ReferenceRange {
  start: number;
  end: number;
}

export interface AlignmentCallability {
  alignment_id: string;
  eligible: boolean;
  status: string;
  callable_bases: number;
  inserted_query_bases: number;
  soft_clipped_query_bases: number;
  reference_exclusion_counts: Record<string, number>;
  disrupted_adjacency_anchors: number[];
  ranges: ReferenceRange[];
}

/** Per-molecule base-callability projection; absence is never inferred from a coverage gap. */
export interface MoleculeCallability {
  molecule_id: string;
  status: string;
  known: boolean;
  basis: 'passing_aligned_reference_bases' | string;
  callable_bases: number;
  callable_fraction: number;
  ranges: ReferenceRange[];
  alignments: AlignmentCallability[];
}

export interface EvidenceCounts {
  alternate: number;
  reference: number;
  event_absent: number;
  callable: number;
  low_quality: number;
  conflict: number;
}

/** Normalized schema 0.6 SNV, small-indel, SV, or complex-path event. */
export interface EvidenceEvent {
  id: string;
  index: number;
  type:
    | 'SNV'
    | 'SMALL_INSERTION'
    | 'SMALL_DELETION'
    | 'COMPLEX_SV_PATH'
    | `SV_${string}`
    | string;
  start: number | null;
  end: number | null;
  length: number;
  ref: string | null;
  alt: string | null;
  normalization: string;
  source_projection: 'variants' | 'small_indels' | 'svs' | 'complex_events' | string;
  negative_evidence_rule:
    | 'callable_base_allele'
    | 'same_fragment_callable_reference_adjacency'
    | 'same_fragment_callable_deleted_span_with_flanks'
    | 'support_only_no_negative_inference'
    | string;
  assessability: 'REFERENCE_AND_ALTERNATE' | 'ALTERNATE_SUPPORT_ONLY' | string;
  component_event_ids: string[];
  supporting_molecule_ids: string[];
  evidence_counts: EvidenceCounts;
}

/** Sparse molecule/event evidence; an absent pair is explicitly NOT_CALLABLE. */
export interface EvidenceObservation {
  id: string;
  molecule_id: string;
  event_id: string;
  alignment_id: string;
  state: ObservationState;
  observed_allele: string | null;
  base_quality: number | null;
  mapping_quality: number;
  strand: '+' | '-' | string;
  evidence_source: string;
  read_position: number | null;
}

/** Columnar page used to keep observation keys and browser allocations bounded. */
export interface EvidenceObservationPage {
  index: number;
  offset: number;
  count: number;
  columns: {
    molecule_id: string[];
    event_id: string[];
    alignment_id: string[];
    state: ObservationState[];
    observed_allele: Array<string | null>;
    base_quality: Array<number | null>;
    mapping_quality: number[];
    strand: string[];
    evidence_source: string[];
    read_position: Array<number | null>;
  };
}

export interface EvidenceEncoding {
  layout: 'paged_columnar_molecule_event' | string;
  scope: 'snv_indel_sv_complex_evidence_rc2' | string;
  observation_storage: 'embedded_columnar_pages' | 'remote_http_pages' | string;
  /** Server transport projection; canonical CLI JSON keeps embedded pages. */
  observation_page_endpoint?: string;
  /** Exact bounded server-side filters over every immutable observation page. */
  observation_search_endpoint?: string;
  missing_pair_state: 'NOT_CALLABLE';
  phase_molecule_policy: 'evidence_eligible_only' | string;
  phase_molecule_reference: 'molecules[].index' | string;
  phase_null_model: 'independent_marginals_within_jointly_callable' | string;
  observation_limit: number;
  observation_count: number;
  observation_page_size: number;
  observation_page_count: number;
  phase_link_limit: number;
}

/** Callable-aware pairwise physical linkage projection. */
export interface PhaseLink {
  id: string;
  event_a_id: string;
  event_b_id: string;
  assessability: 'COMPLETE_FOR_BOTH_EVENTS' | 'SUPPORT_CONDITIONED' | string;
  jointly_callable: number;
  jointly_uncertain: number;
  both_alternate: number;
  a_alternate_b_absent: number;
  a_absent_b_alternate: number;
  neither_alternate: number;
  co_alternate_fraction: number;
  co_alternate_ci95_low: number;
  co_alternate_ci95_high: number;
  expected_co_alternate_fraction: number;
  linkage_delta: number;
  /** Stable references into molecules[].index; avoids repeating long IDs per pair. */
  supporting_molecule_indices: number[];
  uncertain_molecule_indices: number[];
  qc_flags: string[];
}

/** Full analysis payload shared by core, CLI, server, report, and UI. */
export interface MitoAnalysisData {
  metadata: MitoMetadata;
  filter_stats: FilterStats;
  coverage_metrics?: CoverageMetrics;
  genes: GeneAnnotation[];
  coverage: CoverageBin[];
  svs: StructuralVariant[];
  complex_events?: ComplexStructuralEvent[];
  variants?: AggregateVariant[];
  clusters: ClusterSummary[];
  reads: ReadFeature[];
  evidence_encoding?: EvidenceEncoding;
  alignments?: AlignmentFragment[];
  molecules?: AssembledMolecule[];
  callability?: MoleculeCallability[];
  events?: EvidenceEvent[];
  observations?: EvidenceObservation[];
  observation_pages?: EvidenceObservationPage[];
  phase_links?: PhaseLink[];
  architectures?: unknown[];
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
