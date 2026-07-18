import type { MitoAnalysisData, ReadFeature } from '@mito-architect/visualization-lib';
import { ArrowDownUp } from 'lucide-react';
import { useMemo, useState } from 'react';
import { FixedSizeList, ListChildComponentProps } from 'react-window';
import { useMitoStore } from '../lib/store';

interface VariantTableProps {
  data: MitoAnalysisData;
}

type SortKey = 'id' | 'cluster_id' | 'mean_quality' | 'length';
const TABLE_GRID =
  'grid-cols-[minmax(260px,1.7fr)_120px_120px_120px_minmax(240px,1.6fr)]';

export default function VariantTable({ data }: VariantTableProps) {
  const [sortKey, setSortKey] = useState<SortKey>('cluster_id');
  const selectedCluster = useMitoStore((state) => state.selectedCluster);
  const minQuality = useMitoStore((state) => state.minQuality);

  const rows = useMemo(() => {
    return [...data.reads]
      .filter((read) => !read.filtered_numt)
      .filter((read) => selectedCluster === undefined || read.cluster_id === selectedCluster)
      .filter((read) => read.mean_quality >= minQuality)
      .sort((a, b) => compareReads(a, b, sortKey));
  }, [data.reads, minQuality, selectedCluster, sortKey]);

  return (
    <section className="rounded-lg border border-line bg-panel">
      <div className="flex flex-wrap items-center justify-between gap-3 border-b border-line px-4 py-3">
        <h2 className="text-sm font-semibold tracking-normal">Reads and Variants</h2>
        <div className="flex items-center gap-2 text-xs text-muted">
          <ArrowDownUp className="h-4 w-4" aria-hidden />
          {rows.length} rows
        </div>
      </div>
      <div className={`grid ${TABLE_GRID} gap-2 border-b border-line px-4 py-2 text-xs font-semibold uppercase tracking-normal text-muted`}>
        <Header align="left" label="Read" value="id" sortKey={sortKey} setSortKey={setSortKey} />
        <Header align="center" label="Cluster" value="cluster_id" sortKey={sortKey} setSortKey={setSortKey} />
        <Header align="right" label="Q mean" value="mean_quality" sortKey={sortKey} setSortKey={setSortKey} />
        <Header align="right" label="Length" value="length" sortKey={sortKey} setSortKey={setSortKey} />
        <div className="text-left">SNPs / SVs</div>
      </div>
      <FixedSizeList height={320} width="100%" itemCount={rows.length} itemSize={42} itemData={rows}>
        {Row}
      </FixedSizeList>
    </section>
  );
}

function Header({
  label,
  value,
  sortKey,
  setSortKey,
  align
}: {
  label: string;
  value: SortKey;
  sortKey: SortKey;
  setSortKey: (value: SortKey) => void;
  align: 'left' | 'center' | 'right';
}) {
  return (
    <button
      type="button"
      onClick={() => setSortKey(value)}
      className={[
        'w-full',
        align === 'left' ? 'text-left' : align === 'center' ? 'text-center' : 'text-right',
        sortKey === value ? 'text-aqua' : 'text-muted'
      ].join(' ')}
    >
      {label}
    </button>
  );
}

function Row({ index, style, data }: ListChildComponentProps<ReadFeature[]>) {
  const read = data[index];
  return (
    <div
      style={style}
      className={`grid ${TABLE_GRID} items-center gap-2 border-b border-line px-4 text-sm`}
    >
      <div className="truncate font-medium">{read.id}</div>
      <div className="text-center">{read.cluster_id < 0 ? 'Outlier' : `C${read.cluster_id + 1}`}</div>
      <div className="text-right tabular-nums">{read.mean_quality.toFixed(1)}</div>
      <div className="text-right tabular-nums">{read.length}</div>
      <div className="truncate text-muted">
        {read.snps.length} SNPs
        {read.sv_ids.length > 0 ? ` | ${read.sv_ids.join(', ')}` : ''}
      </div>
    </div>
  );
}

function compareReads(a: ReadFeature, b: ReadFeature, sortKey: SortKey): number {
  if (sortKey === 'id') return a.id.localeCompare(b.id);
  return Number(a[sortKey]) - Number(b[sortKey]);
}
