import { MitoCircos } from './MitoCircos';

export { MitoCircos };
export type {
  ClinicalAnnotation,
  ClusterSummary,
  CoverageBin,
  CoverageMetrics,
  FilterStats,
  GeneAnnotation,
  LayerName,
  LayerState,
  MitoAnalysisData,
  MitoCircosOptions,
  MitoMetadata,
  ProteinStructureMapping,
  ReadFeature,
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
