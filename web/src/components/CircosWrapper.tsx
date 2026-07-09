import { MitoCircos } from '@mito-architect/visualization-lib';
import type { LayerName, MitoAnalysisData } from '@mito-architect/visualization-lib';
import { useEffect, useRef } from 'react';

interface CircosWrapperProps {
  data: MitoAnalysisData;
  activeLayers: Record<LayerName, boolean>;
  onClusterSelect: (clusterId: number) => void;
  onSvSelect: (svId: string) => void;
}

export default function CircosWrapper({
  data,
  activeLayers,
  onClusterSelect,
  onSvSelect
}: CircosWrapperProps) {
  const hostRef = useRef<HTMLDivElement>(null);
  const instanceRef = useRef<MitoCircos>();

  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;
    const instance = new MitoCircos(host, data, {
      activeLayers,
      showControls: false
    }).render();
    instanceRef.current = instance;

    const handleCluster = (event: Event) => {
      const detail = (event as CustomEvent<{ id: number }>).detail;
      onClusterSelect(detail.id);
    };
    const handleSv = (event: Event) => {
      const detail = (event as CustomEvent<{ id: string }>).detail;
      onSvSelect(detail.id);
    };
    host.addEventListener('mito:cluster-select', handleCluster);
    host.addEventListener('mito:sv-select', handleSv);
    return () => {
      host.removeEventListener('mito:cluster-select', handleCluster);
      host.removeEventListener('mito:sv-select', handleSv);
      instance.destroy();
    };
  }, [activeLayers, data, onClusterSelect, onSvSelect]);

  return <div ref={hostRef} className="min-h-[520px]" />;
}
