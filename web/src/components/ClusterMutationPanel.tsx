import type {
  ClusterSummary,
  MitoAnalysisData,
  ReadFeature,
  StructuralVariant
} from '@mito-architect/visualization-lib';
import { Crosshair, Dna, LocateFixed } from 'lucide-react';
import { useMemo, useState } from 'react';
import { useMitoStore } from '../lib/store';

interface ClusterMutationPanelProps {
  data: MitoAnalysisData;
}

interface MutationSummary {
  key: string;
  position: number;
  ref: string;
  alt: string;
  support: number;
  frequency: number;
  readIds: string[];
}

interface SvSummary {
  sv: StructuralVariant;
  support: number;
  frequency: number;
  readIds: string[];
}

export default function ClusterMutationPanel({ data }: ClusterMutationPanelProps) {
  const selectedCluster = useMitoStore((state) => state.selectedCluster);
  const selectedSvId = useMitoStore((state) => state.selectedSvId);
  const setSelectedSvId = useMitoStore((state) => state.setSelectedSvId);

  const summary = useMemo(() => {
    const selectedReads = data.reads.filter(
      (read) =>
        !read.filtered_numt &&
        (selectedCluster === undefined || read.cluster_id === selectedCluster)
    );
    const cluster = data.clusters.find((item) => item.id === selectedCluster);
    return {
      cluster,
      reads: selectedReads,
      mutations: summarizeSnps(selectedReads),
      svs: summarizeSvs(selectedReads, data.svs)
    };
  }, [data.reads, data.svs, data.clusters, selectedCluster]);

  const label = selectedCluster === undefined ? 'All clusters' : summary.cluster?.label ?? `H${selectedCluster + 1}`;
  const denominator = Math.max(1, summary.reads.length);
  const selectedSv = data.svs.find((sv) => sv.id === selectedSvId);

  return (
    <section className="rounded-lg border border-line bg-panel">
      <div className="flex flex-wrap items-start justify-between gap-3 border-b border-line px-4 py-3">
        <div>
          <h2 className="flex items-center gap-2 text-sm font-semibold tracking-normal">
            <Crosshair className="h-4 w-4 text-aqua" aria-hidden />
            Cluster Mutation Drill-Down
          </h2>
          <p className="mt-1 text-sm text-muted">
            {label} | {summary.reads.length} molecules | frequencies are within the current selection
          </p>
        </div>
        <div className="rounded-md border border-line bg-panel2 px-3 py-2 text-sm font-semibold text-aqua">
          {selectedCluster === undefined ? 'Global view' : 'Selected from plot/sidebar'}
        </div>
      </div>

      <div className="grid gap-4 p-4 xl:grid-cols-[minmax(0,1.2fr)_minmax(0,0.8fr)]">
        <div className="grid gap-4">
          <GenomeDistribution
            referenceLength={data.metadata.reference_length}
            mutations={summary.mutations}
            svs={summary.svs}
          />

          <div className="overflow-hidden rounded-md border border-line">
            <div className="grid grid-cols-[0.9fr_0.7fr_0.7fr_1.4fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
              <div>Mutation</div>
              <div>Support</div>
              <div>Frequency</div>
              <div>Reads</div>
            </div>
            <div className="max-h-64 overflow-auto scrollbar-thin">
              {summary.mutations.length === 0 ? (
                <EmptyRow text="No SNPs in this selection" />
              ) : (
                summary.mutations.map((mutation) => (
                  <div
                    key={mutation.key}
                    className="grid grid-cols-[0.9fr_0.7fr_0.7fr_1.4fr] gap-2 border-b border-line px-3 py-2 text-sm last:border-0"
                  >
                    <div className="font-semibold">
                      {mutation.position} {mutation.ref}&gt;{mutation.alt}
                    </div>
                    <div>{mutation.support}/{denominator}</div>
                    <div>{formatPercent(mutation.frequency)}</div>
                    <div className="truncate text-muted" title={mutation.readIds.join(', ')}>
                      {mutation.readIds.join(', ')}
                    </div>
                  </div>
                ))
              )}
            </div>
          </div>
        </div>

        <div className="grid content-start gap-4">
          <div className="rounded-md border border-line bg-panel2 p-4">
            <h3 className="flex items-center gap-2 text-sm font-semibold tracking-normal">
              <Dna className="h-4 w-4 text-aqua" aria-hidden />
              Cluster Signature
            </h3>
            <dl className="mt-3 grid grid-cols-2 gap-3 text-sm">
              <Metric label="Molecules" value={summary.reads.length.toString()} />
              <Metric label="SNP sites" value={summary.mutations.length.toString()} />
              <Metric label="SV calls" value={summary.svs.length.toString()} />
              <Metric
                label="Mean Q"
                value={
                  summary.reads.length === 0
                    ? '0.0'
                    : (
                        summary.reads.reduce((sum, read) => sum + read.mean_quality, 0) /
                        summary.reads.length
                      ).toFixed(1)
                }
              />
            </dl>
          </div>

          <div className="overflow-hidden rounded-md border border-line">
            <div className="border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
              Structural Variants
            </div>
            <div className="max-h-56 overflow-auto scrollbar-thin">
              {summary.svs.length === 0 ? (
                <EmptyRow text="No SVs in this selection" />
              ) : (
                summary.svs.map((item) => (
                  <button
                    key={item.sv.id}
                    type="button"
                    onClick={() => setSelectedSvId(item.sv.id)}
                    className={[
                      'w-full border-b px-3 py-2 text-left text-sm last:border-0',
                      selectedSvId === item.sv.id
                        ? 'border-aqua bg-teal-300/10'
                        : 'border-line hover:bg-panel2'
                    ].join(' ')}
                  >
                    <div className="flex items-center justify-between gap-3">
                      <span className="font-semibold">
                        {item.sv.type} {item.sv.start}-{item.sv.end}
                      </span>
                      <span className={item.sv.known_event ? 'text-amber' : 'text-magenta'}>
                        {formatPercent(item.frequency)}
                      </span>
                    </div>
                    <div className="mt-1 text-muted">
                      {item.support}/{denominator} molecules | {item.sv.length} bp
                    </div>
                  </button>
                ))
              )}
            </div>
          </div>

          <div className="rounded-md border border-line bg-panel2 p-4">
            <h3 className="text-sm font-semibold tracking-normal">Selected SV</h3>
            {selectedSv ? (
              <div className="mt-3 text-sm">
                <div className="font-semibold">
                  {selectedSv.type} {selectedSv.start}-{selectedSv.end}
                </div>
                <div className="mt-1 text-muted">
                  {selectedSv.length} bp | {selectedSv.supporting_reads.length} supporting reads
                </div>
                <div className="mt-2 text-muted">
                  {selectedSv.known_event ? 'Known common deletion candidate' : 'Novel or uncatalogued event'}
                </div>
              </div>
            ) : (
              <div className="mt-3 text-sm text-muted">Click an SV chord in the plot or an SV row here.</div>
            )}
          </div>
        </div>
      </div>
    </section>
  );
}

function GenomeDistribution({
  referenceLength,
  mutations,
  svs
}: {
  referenceLength: number;
  mutations: MutationSummary[];
  svs: SvSummary[];
}) {
  const width = 900;
  const height = 94;
  const labelPadding = 34;
  const axisWidth = width - labelPadding * 2;
  const [hoveredSvId, setHoveredSvId] = useState<string>();
  const coordinateSpan = Math.max(1, referenceLength - 1);
  const x = (position: number) =>
    labelPadding +
    ((Math.max(1, Math.min(referenceLength, position)) - 1) / coordinateSpan) * axisWidth;
  const tickPositions = [0, 0.25, 0.5, 0.75, 1].map((fraction) =>
    1 + Math.round(coordinateSpan * fraction)
  );

  return (
    <div className="rounded-md border border-line bg-panel2 p-3">
      <div className="mb-2 flex items-center justify-between gap-3 text-sm">
        <h3 className="flex items-center gap-2 font-semibold tracking-normal">
          <LocateFixed className="h-4 w-4 text-aqua" aria-hidden />
          Genome Distribution
        </h3>
        <span className="text-xs text-muted">rCRS coordinates, 1-{referenceLength}</span>
      </div>
      <svg viewBox={`0 0 ${width} ${height}`} className="h-24 w-full overflow-visible">
        <line
          x1={labelPadding}
          x2={width - labelPadding}
          y1="42"
          y2="42"
          stroke="#334155"
          strokeWidth="8"
          strokeLinecap="round"
        />
        {tickPositions.map((position) => (
          <g key={position} transform={`translate(${x(position)} 0)`}>
            <line y1="28" y2="58" stroke="#64748b" strokeWidth="1" />
            <text y="76" textAnchor="middle" fill="#94a3b8" fontSize="16">
              {position}
            </text>
          </g>
        ))}
        {svs.map((item, index) => (
          <rect
            key={item.sv.id}
            x={x(item.sv.start)}
            y={12 + index * 8}
            width={Math.max(4, x(item.sv.end) - x(item.sv.start))}
            height="10"
            rx="2"
            fill={item.sv.known_event ? '#f97316' : '#e879f9'}
            opacity={hoveredSvId === item.sv.id ? 1 : 0.38 + item.frequency * 0.55}
            stroke={hoveredSvId === item.sv.id ? '#f8fafc' : 'transparent'}
            strokeWidth={hoveredSvId === item.sv.id ? 3 : 0}
            onPointerEnter={() => setHoveredSvId(item.sv.id)}
            onPointerLeave={() => setHoveredSvId(undefined)}
          />
        ))}
        {mutations.map((mutation) => (
          <g key={mutation.key} transform={`translate(${x(mutation.position)} 0)`}>
            <line
              y1="24"
              y2="60"
              stroke={mutation.alt === 'T' ? '#ef4444' : mutation.alt === 'G' ? '#f59e0b' : '#22c55e'}
              strokeWidth={Math.max(2, 2 + mutation.frequency * 8)}
              strokeLinecap="round"
            />
            <circle
              cy="42"
              r={Math.max(3, 3 + mutation.frequency * 6)}
              fill={mutation.alt === 'T' ? '#ef4444' : mutation.alt === 'G' ? '#f59e0b' : '#22c55e'}
              stroke="#08111f"
              strokeWidth="1"
            />
          </g>
        ))}
      </svg>
    </div>
  );
}

function summarizeSnps(reads: ReadFeature[]): MutationSummary[] {
  const byMutation = new Map<string, MutationSummary>();
  for (const read of reads) {
    for (const snp of read.snps) {
      const key = `${snp.position}:${snp.ref}:${snp.alt}`;
      const current =
        byMutation.get(key) ??
        ({
          key,
          position: snp.position,
          ref: snp.ref,
          alt: snp.alt,
          support: 0,
          frequency: 0,
          readIds: []
        } satisfies MutationSummary);
      current.support += 1;
      current.readIds.push(read.id);
      byMutation.set(key, current);
    }
  }
  return [...byMutation.values()]
    .map((mutation) => ({
      ...mutation,
      frequency: reads.length === 0 ? 0 : mutation.support / reads.length
    }))
    .sort((a, b) => a.position - b.position || b.support - a.support);
}

function summarizeSvs(reads: ReadFeature[], svs: StructuralVariant[]): SvSummary[] {
  const byId = new Map(svs.map((sv) => [sv.id, sv]));
  const support = new Map<string, { count: number; readIds: string[] }>();
  for (const read of reads) {
    for (const svId of read.sv_ids) {
      const current = support.get(svId) ?? { count: 0, readIds: [] };
      current.count += 1;
      current.readIds.push(read.id);
      support.set(svId, current);
    }
  }
  return [...support.entries()]
    .map(([svId, item]) => {
      const sv =
        byId.get(svId) ??
        ({
          id: svId,
          type: 'sv',
          start: 1,
          end: 1,
          length: 0,
          known_event: false,
          supporting_reads: item.readIds
        } satisfies StructuralVariant);
      return {
        sv,
        support: item.count,
        frequency: reads.length === 0 ? 0 : item.count / reads.length,
        readIds: item.readIds
      };
    })
    .sort((a, b) => a.sv.start - b.sv.start || b.support - a.support);
}

function Metric({ label, value }: { label: string; value: string }) {
  return (
    <div>
      <dt className="text-xs uppercase tracking-normal text-muted">{label}</dt>
      <dd className="mt-1 text-lg font-semibold">{value}</dd>
    </div>
  );
}

function EmptyRow({ text }: { text: string }) {
  return <div className="px-3 py-4 text-sm text-muted">{text}</div>;
}

function formatPercent(value: number): string {
  return `${Math.round(value * 1000) / 10}%`;
}
