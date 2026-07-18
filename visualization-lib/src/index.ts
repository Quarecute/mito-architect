import { MitoCircos } from './MitoCircos';

export { MitoCircos };
export type {
  AggregateVariant,
  AlignmentCallability,
  AlignmentFragment,
  AlleleQualitySummary,
  AssembledMolecule,
  ClinicalAnnotation,
  ClinicalAssertion,
  ClusterSummary,
  ComplexStructuralEvent,
  CoverageBin,
  CoverageMetrics,
  FilterStats,
  EvidenceEncoding,
  EvidenceCounts,
  EvidenceEvent,
  EvidenceObservation,
  GeneAnnotation,
  LayerName,
  LayerState,
  MitoAnalysisData,
  MitoCircosOptions,
  MitoMetadata,
  MoleculeCallability,
  ObservationState,
  PhaseLink,
  ProteinStructureMapping,
  ReadFeature,
  ReferenceRange,
  SnpCall,
  StructuralVariant
} from './types';

declare global {
  interface Window {
    MitoCircos?: typeof MitoCircos;
  }
}

if (typeof window !== 'undefined') {
  window.MitoCircos = MitoCircos;
}
