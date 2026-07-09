import type { ClusterSummary, GeneAnnotation, LayerName } from '@mito-architect/visualization-lib';
import { Layers, Search, SlidersHorizontal } from 'lucide-react';
import { useMemo } from 'react';
import { useMitoStore } from '../lib/store';

const LAYERS: Array<[LayerName, string]> = [
  ['genes', 'Genes'],
  ['coverage', 'Coverage'],
  ['clusters', 'Clusters'],
  ['svs', 'SVs'],
  ['snps', 'SNPs']
];

interface SidebarProps {
  genes: GeneAnnotation[];
  clusters: ClusterSummary[];
}

export default function Sidebar({ genes, clusters }: SidebarProps) {
  const activeLayers = useMitoStore((state) => state.activeLayers);
  const setLayer = useMitoStore((state) => state.setLayer);
  const geneQuery = useMitoStore((state) => state.geneQuery);
  const setGeneQuery = useMitoStore((state) => state.setGeneQuery);
  const minQuality = useMitoStore((state) => state.minQuality);
  const setMinQuality = useMitoStore((state) => state.setMinQuality);
  const selectedCluster = useMitoStore((state) => state.selectedCluster);
  const setSelectedCluster = useMitoStore((state) => state.setSelectedCluster);

  const filteredGenes = useMemo(() => {
    const query = geneQuery.trim().toLowerCase();
    return query ? genes.filter((gene) => gene.name.toLowerCase().includes(query)) : genes;
  }, [geneQuery, genes]);

  return (
    <aside className="grid content-start gap-4">
      <section className="glass-panel rounded-lg border border-line p-4 shadow-tool">
        <h2 className="mb-3 flex items-center gap-2 text-sm font-semibold tracking-normal">
          <Layers className="h-4 w-4 text-aqua" aria-hidden />
          Layers
        </h2>
        <div className="grid grid-cols-2 gap-2">
          {LAYERS.map(([layer, label]) => (
            <label
              key={layer}
              className={[
                'flex cursor-pointer items-center justify-between gap-3 rounded-md border px-3 py-2 text-sm',
                activeLayers[layer] ? 'border-aqua/70 bg-aqua/10 text-text' : 'border-line bg-panel2 text-muted'
              ].join(' ')}
            >
              <span>{label}</span>
              <input
                type="checkbox"
                checked={activeLayers[layer]}
                onChange={(event) => setLayer(layer, event.currentTarget.checked)}
                className="h-4 w-4 accent-teal-300"
              />
            </label>
          ))}
        </div>
      </section>

      <section className="glass-panel rounded-lg border border-line p-4 shadow-tool">
        <h2 className="mb-3 flex items-center gap-2 text-sm font-semibold tracking-normal">
          <Search className="h-4 w-4 text-aqua" aria-hidden />
          Gene Search
        </h2>
        <input
          value={geneQuery}
          onChange={(event) => setGeneQuery(event.target.value)}
          placeholder="MT-ND5"
          className="w-full rounded-md border border-line bg-panel2 px-3 py-2 text-sm outline-none focus:border-aqua focus:shadow-focus"
        />
        <div className="mt-3 max-h-40 overflow-auto pr-1 text-sm scrollbar-thin">
          {filteredGenes.map((gene) => (
            <div key={gene.name} className="flex justify-between border-b border-line py-2 last:border-0">
              <span className="font-medium">{gene.name}</span>
              <span className="text-muted">{gene.start}-{gene.end}</span>
            </div>
          ))}
        </div>
      </section>

      <section className="glass-panel rounded-lg border border-line p-4 shadow-tool">
        <h2 className="mb-3 flex items-center gap-2 text-sm font-semibold tracking-normal">
          <SlidersHorizontal className="h-4 w-4 text-aqua" aria-hidden />
          Read Filter
        </h2>
        <label className="text-xs uppercase tracking-normal text-muted">Minimum mean quality</label>
        <input
          type="range"
          min={0}
          max={45}
          step={1}
          value={minQuality}
          onChange={(event) => setMinQuality(Number(event.currentTarget.value))}
          className="mt-2 w-full accent-teal-300"
        />
        <div className="mt-2 text-sm font-semibold">{minQuality}</div>
      </section>

      <section className="glass-panel rounded-lg border border-line p-4 shadow-tool">
        <h2 className="mb-3 text-sm font-semibold tracking-normal">Clusters</h2>
        <div className="grid gap-2">
          <button
            type="button"
            onClick={() => setSelectedCluster(undefined)}
            className={[
              'rounded-md border px-3 py-2 text-left text-sm',
              selectedCluster === undefined ? 'border-aqua bg-aqua/10 text-aqua shadow-focus' : 'border-line bg-panel2 text-muted'
            ].join(' ')}
          >
            All clusters
          </button>
          {clusters.map((cluster) => (
            <button
              key={cluster.id}
              type="button"
              onClick={() => setSelectedCluster(cluster.id)}
              className={[
                'rounded-md border px-3 py-2 text-left text-sm',
                selectedCluster === cluster.id ? 'border-aqua bg-aqua/10 text-aqua shadow-focus' : 'border-line bg-panel2 text-muted hover:border-aqua/50'
              ].join(' ')}
            >
              <span className="flex items-center justify-between gap-3">
                <span>
                  <span className="font-semibold">{cluster.label}</span>
                  <span className="ml-2">{cluster.size} reads</span>
                </span>
                {cluster.haplogroup ? (
                  <span className="text-xs text-muted">{cluster.haplogroup}</span>
                ) : null}
              </span>
            </button>
          ))}
        </div>
      </section>
    </aside>
  );
}
