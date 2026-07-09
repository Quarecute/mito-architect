import type { LayerName, MitoAnalysisData } from '@mito-architect/visualization-lib';
import { create } from 'zustand';

type LayerMap = Record<LayerName, boolean>;

interface MitoState {
  selectedFile?: File;
  jobId?: string;
  data?: MitoAnalysisData;
  activeLayers: LayerMap;
  selectedCluster?: number;
  selectedSvId?: string;
  geneQuery: string;
  minQuality: number;
  setSelectedFile: (file?: File) => void;
  setJobId: (jobId?: string) => void;
  setData: (data?: MitoAnalysisData) => void;
  setLayer: (layer: LayerName, active: boolean) => void;
  setSelectedCluster: (clusterId?: number) => void;
  setSelectedSvId: (svId?: string) => void;
  setGeneQuery: (query: string) => void;
  setMinQuality: (quality: number) => void;
}

export const useMitoStore = create<MitoState>((set) => ({
  activeLayers: {
    genes: true,
    coverage: true,
    clusters: true,
    svs: true,
    snps: true
  },
  geneQuery: '',
  minQuality: 0,
  setSelectedFile: (selectedFile) => set({ selectedFile }),
  setJobId: (jobId) => set({ jobId }),
  setData: (data) => set({ data }),
  setLayer: (layer, active) =>
    set((state) => ({
      activeLayers: {
        ...state.activeLayers,
        [layer]: active
      }
    })),
  setSelectedCluster: (selectedCluster) => set({ selectedCluster }),
  setSelectedSvId: (selectedSvId) => set({ selectedSvId }),
  setGeneQuery: (geneQuery) => set({ geneQuery }),
  setMinQuality: (minQuality) => set({ minQuality })
}));
