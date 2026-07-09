import type {
  ClinicalAnnotation,
  MitoAnalysisData,
  ProteinStructureMapping,
  ReadFeature
} from '@mito-architect/visualization-lib';
import {
  Activity,
  Atom,
  ChartNoAxesCombined,
  Clock3,
  Database,
  Dna,
  GitBranch,
  ScanSearch,
  Table2
} from 'lucide-react';
import { useEffect, useMemo, useState } from 'react';
import NglProteinViewer, { type ProteinViewerVariant } from './NglProteinViewer';

interface AnalysisModulesPanelProps {
  data: MitoAnalysisData;
}

type ModuleTab =
  | 'haplogroups'
  | 'clinical'
  | 'coverage'
  | 'comparison'
  | 'singleCell'
  | 'protein'
  | 'timeline';

interface VariantSummary {
  key: string;
  label: string;
  type: string;
  position: number;
  support: number;
  frequency: number;
  gene?: string;
  consequence?: string;
  protein?: string;
  residue?: string;
  annotation?: ClinicalAnnotation;
  structure?: ProteinStructureMapping;
}

const MODULES: Array<{ id: ModuleTab; label: string; Icon: typeof Dna }> = [
  { id: 'haplogroups', label: 'Haplogroups', Icon: Dna },
  { id: 'clinical', label: 'Clinical', Icon: Database },
  { id: 'coverage', label: 'Coverage', Icon: Activity },
  { id: 'comparison', label: 'Multi-sample', Icon: Table2 },
  { id: 'singleCell', label: 'Single-cell', Icon: ScanSearch },
  { id: 'protein', label: '3D Protein', Icon: Atom },
  { id: 'timeline', label: 'Timeline', Icon: Clock3 }
];

export default function AnalysisModulesPanel({ data }: AnalysisModulesPanelProps) {
  const [activeTab, setActiveTab] = useState<ModuleTab>('haplogroups');
  const passedReads = useMemo(() => data.reads.filter((read) => !read.filtered_numt), [data.reads]);
  const variants = useMemo(() => summarizeVariants(passedReads, data), [passedReads, data]);
  const coverageMetrics = useMemo(() => summarizeCoverage(data), [data]);
  const active = MODULES.find((module) => module.id === activeTab) ?? MODULES[0];

  return (
    <section className="rounded-lg border border-line bg-panel">
      <div className="flex flex-wrap items-center justify-between gap-3 border-b border-line px-4 py-3">
        <h2 className="flex items-center gap-2 text-sm font-semibold tracking-normal">
          <ChartNoAxesCombined className="h-4 w-4 text-aqua" aria-hidden />
          Analysis Modules
        </h2>
        <div className="flex flex-wrap gap-2">
          {MODULES.map(({ id, label, Icon }) => (
            <button
              key={id}
              type="button"
              onClick={() => setActiveTab(id)}
              className={[
                'inline-flex items-center gap-2 rounded-md border px-3 py-2 text-sm',
                activeTab === id ? 'border-aqua text-aqua' : 'border-line text-muted hover:bg-panel2'
              ].join(' ')}
            >
              <Icon className="h-4 w-4" aria-hidden />
              {label}
            </button>
          ))}
        </div>
      </div>

      <div className="p-4">
        <div className="mb-4 flex items-center gap-2 text-sm font-semibold">
          <active.Icon className="h-4 w-4 text-aqua" aria-hidden />
          {active.label}
        </div>
        {activeTab === 'haplogroups' && <HaplogroupPanel data={data} />}
        {activeTab === 'clinical' && <ClinicalPanel variants={variants} />}
        {activeTab === 'coverage' && <CoveragePanel data={data} metrics={coverageMetrics} />}
        {activeTab === 'comparison' && <ComparisonPanel data={data} variants={variants} />}
        {activeTab === 'singleCell' && <SingleCellPanel reads={passedReads} />}
        {activeTab === 'protein' && <ProteinPanel variants={variants} />}
        {activeTab === 'timeline' && <TimelinePanel data={data} variants={variants} />}
      </div>
    </section>
  );
}

function HaplogroupPanel({ data }: { data: MitoAnalysisData }) {
  return (
    <div className="overflow-hidden rounded-md border border-line">
      <div className="grid grid-cols-[0.7fr_0.9fr_0.6fr_1.4fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
        <div>Cluster</div>
        <div>Haplogroup</div>
        <div>Molecules</div>
        <div>SV signature</div>
      </div>
      {data.clusters.length === 0 ? (
        <EmptyRow text="No clusters" />
      ) : (
        data.clusters.map((cluster) => (
          <div
            key={cluster.id}
            className="grid grid-cols-[0.7fr_0.9fr_0.6fr_1.4fr] gap-2 border-b border-line px-3 py-2 text-sm last:border-0"
          >
            <div className="font-semibold">{cluster.label}</div>
            <div className={cluster.haplogroup === 'unassigned' ? 'text-muted' : 'text-aqua'}>
              {cluster.haplogroup ?? 'unassigned'}
            </div>
            <div>{cluster.size}</div>
            <div className="truncate text-muted">
              {cluster.sv_signature.length === 0
                ? 'none'
                : cluster.sv_signature.map((sv) => `${sv.sv_id} n=${sv.support}`).join(', ')}
            </div>
          </div>
        ))
      )}
    </div>
  );
}

function ClinicalPanel({ variants }: { variants: VariantSummary[] }) {
  return (
    <div className="overflow-hidden rounded-md border border-line">
      <div className="grid grid-cols-[1fr_0.6fr_0.6fr_1fr_1fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
        <div>Variant</div>
        <div>Support</div>
        <div>Frequency</div>
        <div>Pathogenicity</div>
        <div>Phenotype / References</div>
      </div>
      {variants.length === 0 ? (
        <EmptyRow text="No variants" />
      ) : (
        variants.map((variant) => (
          <div
            key={variant.key}
            className="grid grid-cols-[1fr_0.6fr_0.6fr_1fr_1fr] gap-2 border-b border-line px-3 py-2 text-sm last:border-0"
          >
            <div className="font-semibold">{variant.label}</div>
            <div>{variant.support}</div>
            <div>{formatPercent(variant.frequency)}</div>
            <div className={variant.annotation?.pathogenicity ? 'text-amber' : 'text-muted'}>
              {variant.annotation?.pathogenicity ?? 'not annotated'}
            </div>
            <div className="truncate text-muted">
              {variant.annotation?.phenotype ?? 'local MITOMAP/ClinVar cache required'}
              {variant.annotation?.sources?.length ? ` | ${variant.annotation.sources.join('+')}` : ''}
              {variant.annotation?.references?.length ? ` | ${variant.annotation.references.join(', ')}` : ''}
            </div>
          </div>
        ))
      )}
    </div>
  );
}

function CoveragePanel({
  data,
  metrics
}: {
  data: MitoAnalysisData;
  metrics: ReturnType<typeof summarizeCoverage>;
}) {
  const preview = data.coverage.slice(0, 80);
  return (
    <div className="grid gap-4 lg:grid-cols-[320px_minmax(0,1fr)]">
      <div className="grid gap-3 sm:grid-cols-3 lg:grid-cols-1">
        <Metric label="Mean depth" value={metrics.meanDepth.toFixed(2)} />
        <Metric label="Bins >20x" value={formatPercent(metrics.pctGt20 / 100)} />
        <Metric label="Max depth" value={metrics.maxDepth.toString()} />
      </div>
      <div className="rounded-md border border-line bg-panel2 p-3">
        <div className="mb-2 text-xs font-semibold uppercase tracking-normal text-muted">Coverage histogram</div>
        <div className="flex h-28 items-end gap-px overflow-hidden">
          {preview.map((bin) => (
            <div
              key={`${bin.start}-${bin.end}`}
              className="min-w-[3px] flex-1 rounded-t-sm bg-sky"
              title={`${bin.start}-${bin.end}: ${bin.depth}`}
              style={{ height: `${Math.max(4, (bin.depth / Math.max(1, metrics.maxDepth)) * 104)}px` }}
            />
          ))}
        </div>
      </div>
    </div>
  );
}

function ComparisonPanel({ data, variants }: { data: MitoAnalysisData; variants: VariantSummary[] }) {
  const sample = data.metadata.sample || 'sample';
  return (
    <div className="grid gap-4 xl:grid-cols-[minmax(0,1fr)_320px]">
      <div className="overflow-hidden rounded-md border border-line">
        <div className="grid grid-cols-[1.3fr_0.7fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <div>Variant</div>
          <div>{sample}</div>
        </div>
        {variants.slice(0, 14).map((variant) => (
          <div
            key={variant.key}
            className="grid grid-cols-[1.3fr_0.7fr] gap-2 border-b border-line px-3 py-2 text-sm last:border-0"
          >
            <div className="truncate font-semibold">{variant.label}</div>
            <div className={variant.support > 0 ? 'text-aqua' : 'text-muted'}>{variant.support > 0 ? 'present' : 'absent'}</div>
          </div>
        ))}
        {variants.length === 0 && <EmptyRow text="No variant matrix rows" />}
      </div>
      <div className="rounded-md border border-line bg-panel2 p-4 text-sm">
        <Metric label="Samples" value="1" />
        <Metric label="Jaccard self" value="1.000" />
        <Metric label="UPGMA status" value="waiting for multiple samples" />
      </div>
    </div>
  );
}

function SingleCellPanel({ reads }: { reads: ReadFeature[] }) {
  const cells = new Map<string, number>();
  for (const read of reads) {
    const match = read.id.match(/\b(?:CB|CR):Z:([^\s]+)/);
    if (match) {
      cells.set(match[1], (cells.get(match[1]) ?? 0) + 1);
    }
  }
  return (
    <div className="grid gap-3 sm:grid-cols-3">
      <Metric label="Cell barcodes" value={cells.size.toString()} />
      <Metric label="Barcode reads" value={[...cells.values()].reduce((sum, count) => sum + count, 0).toString()} />
      <Metric label="Cluster status" value={cells.size > 0 ? 'ready' : 'no CB/CR tags'} />
    </div>
  );
}

function ProteinPanel({ variants }: { variants: VariantSummary[] }) {
  const proteinVariants = variants.filter((variant) => variant.protein || variant.structure?.structure_id);
  const [selectedKey, setSelectedKey] = useState<string>();
  const selected = proteinVariants.find((variant) => variant.key === selectedKey) ?? proteinVariants[0];

  useEffect(() => {
    if (!selectedKey && proteinVariants[0]) {
      setSelectedKey(proteinVariants[0].key);
    } else if (selectedKey && !proteinVariants.some((variant) => variant.key === selectedKey)) {
      setSelectedKey(proteinVariants[0]?.key);
    }
  }, [proteinVariants, selectedKey]);

  return (
    <div className="grid gap-4 lg:grid-cols-[320px_minmax(0,1fr)]">
      <div className="grid content-start gap-3">
        <Metric label="Mapped residues" value={proteinVariants.length.toString()} />
        <Metric label="Viewer" value={selected?.structure?.structure_id ?? 'NGL-ready'} />
        <Metric label="Example" value={selected ? selected.residue ?? selected.label : 'no protein variant'} />
        <div className="overflow-hidden rounded-md border border-line bg-panel2">
          <div className="border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
            Protein variants
          </div>
          <div className="max-h-52 overflow-auto scrollbar-thin">
            {proteinVariants.length === 0 ? (
              <EmptyRow text="No protein-mapped variants" />
            ) : (
              proteinVariants.map((variant) => (
                <button
                  key={variant.key}
                  type="button"
                  onClick={() => setSelectedKey(variant.key)}
                  className={[
                    'w-full border-b px-3 py-2 text-left text-sm last:border-0',
                    selected?.key === variant.key ? 'border-aqua bg-teal-300/10' : 'border-line hover:bg-panel'
                  ].join(' ')}
                >
                  <div className="font-semibold">{variant.protein ?? variant.gene}</div>
                  <div className="mt-1 truncate text-muted">{variant.residue ?? variant.label}</div>
                </button>
              ))
            )}
          </div>
        </div>
      </div>
      <div className="rounded-md border border-line bg-panel2 p-4">
        {selected ? (
          <div className="grid gap-4">
            <NglProteinViewer variant={toProteinViewerVariant(selected)} />
            <dl className="grid gap-3 sm:grid-cols-2 xl:grid-cols-4">
              <Metric label="Variant" value={selected.label} />
              <Metric label="Residue" value={selected.residue ?? 'mapped'} />
              <Metric label="Structure" value={selected.structure?.structure_id ?? 'local model'} />
              <Metric label="Frequency" value={formatPercent(selected.frequency)} />
            </dl>
          </div>
        ) : (
          <div className="text-sm text-muted">
            Non-synonymous consequence calls and local structure mappings are required before residue highlights are populated.
          </div>
        )}
      </div>
    </div>
  );
}

function toProteinViewerVariant(variant: VariantSummary): ProteinViewerVariant {
  return {
    key: variant.key,
    label: variant.label,
    gene: variant.gene,
    protein: variant.protein,
    residue: variant.residue,
    frequency: variant.frequency,
    structure: variant.structure
  };
}

function TimelinePanel({ data, variants }: { data: MitoAnalysisData; variants: VariantSummary[] }) {
  const metadata = data.metadata as typeof data.metadata & { collection_date?: string };
  return (
    <div className="grid gap-4 lg:grid-cols-[320px_minmax(0,1fr)]">
      <div className="grid gap-3">
        <Metric label="Samples with dates" value={metadata.collection_date ? '1' : '0'} />
        <Metric label="Tree tips" value="1" />
        <Metric label="Variant events" value={variants.length.toString()} />
      </div>
      <div className="rounded-md border border-line bg-panel2 p-4">
        <div className="flex items-center gap-3 text-sm">
          <GitBranch className="h-4 w-4 text-aqua" aria-hidden />
          <span className="font-semibold">{data.metadata.sample}</span>
          <span className="text-muted">{metadata.collection_date ?? 'date not provided'}</span>
        </div>
        <div className="mt-4 h-1 rounded-full bg-line">
          <div className="h-1 w-full rounded-full bg-aqua" />
        </div>
      </div>
    </div>
  );
}

function Metric({ label, value }: { label: string; value: string }) {
  return (
    <div className="min-w-0 rounded-md border border-line bg-panel2 p-3">
      <div className="text-xs font-semibold uppercase tracking-normal text-muted">{label}</div>
      <div className="mt-1 break-words text-lg font-semibold leading-snug tracking-normal [overflow-wrap:anywhere]">
        {value}
      </div>
    </div>
  );
}

function EmptyRow({ text }: { text: string }) {
  return <div className="px-3 py-3 text-sm text-muted">{text}</div>;
}

function summarizeVariants(reads: ReadFeature[], data: MitoAnalysisData): VariantSummary[] {
  const denominator = Math.max(1, reads.length);
  const snps = new Map<string, VariantSummary>();

  for (const read of reads) {
    for (const snp of read.snps) {
      const key = `snp:${snp.position}:${snp.ref}:${snp.alt}`;
      const existing = snps.get(key);
      if (existing) {
        existing.support += 1;
        existing.frequency = existing.support / denominator;
      } else {
        snps.set(key, {
          key,
          label: `${snp.position} ${snp.ref}>${snp.alt}`,
          type: 'SNP',
          position: snp.position,
          support: 1,
          frequency: 1 / denominator,
          gene: snp.gene,
          consequence: snp.consequence,
          protein: snp.protein,
          residue: snp.residue,
          annotation: snp.annotation,
          structure: snp.structure
        });
      }
    }
  }

  const svs = data.svs.map((sv) => ({
    key: `sv:${sv.id}`,
    label: `${sv.type} ${sv.start}-${sv.end}`,
    type: sv.type,
    position: sv.start,
    support: sv.supporting_reads.length,
    frequency: sv.supporting_reads.length / denominator,
    annotation: sv.annotation
  }));

  return [...snps.values(), ...svs].sort((a, b) => a.position - b.position || a.label.localeCompare(b.label));
}

function summarizeCoverage(data: MitoAnalysisData) {
  if (data.coverage_metrics) {
    return {
      meanDepth: data.coverage_metrics.mean_depth,
      pctGt20: data.coverage_metrics.pct_sites_gt20x ?? data.coverage_metrics.pct_bins_gt20x ?? 0,
      maxDepth: data.coverage_metrics.max_depth
    };
  }

  const bins = data.coverage;
  const maxDepth = Math.max(0, ...bins.map((bin) => bin.depth));
  const meanDepth = bins.length === 0 ? 0 : bins.reduce((sum, bin) => sum + bin.depth, 0) / bins.length;
  const pctGt20 = bins.length === 0 ? 0 : (bins.filter((bin) => bin.depth > 20).length * 100) / bins.length;
  return { meanDepth, pctGt20, maxDepth };
}

function formatPercent(value: number): string {
  return `${(value * 100).toFixed(1)}%`;
}
