import { useQuery } from '@tanstack/react-query';
import { Activity, Binary, Dna, Filter, Gauge } from 'lucide-react';
import { useCallback, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import AnalysisModulesPanel from '../components/AnalysisModulesPanel';
import ClusterMutationPanel from '../components/ClusterMutationPanel';
import CircosWrapper from '../components/CircosWrapper';
import ExportButtons from '../components/ExportButtons';
import Sidebar from '../components/Sidebar';
import VariantTable from '../components/VariantTable';
import { getResult } from '../lib/api';
import { demoData } from '../lib/demoData';
import { useMitoStore } from '../lib/store';

export default function ResultPage() {
  const { jobId } = useParams<{ jobId: string }>();
  const isDemo = jobId === 'demo';
  const activeLayers = useMitoStore((state) => state.activeLayers);
  const setSelectedCluster = useMitoStore((state) => state.setSelectedCluster);
  const setSelectedSvId = useMitoStore((state) => state.setSelectedSvId);
  const setData = useMitoStore((state) => state.setData);

  const result = useQuery({
    queryKey: ['result', jobId],
    queryFn: ({ signal }) => getResult(jobId!, signal),
    enabled: Boolean(jobId) && !isDemo
  });

  useEffect(() => {
    setData(isDemo ? demoData : result.data);
  }, [isDemo, result.data, setData]);

  const handleClusterSelect = useCallback(
    (clusterId: number) => {
      setSelectedCluster(clusterId);
    },
    [setSelectedCluster]
  );
  const handleSvSelect = useCallback(
    (svId: string) => {
      setSelectedSvId(svId);
    },
    [setSelectedSvId]
  );

  if (!isDemo && result.isLoading) {
    return <ResultSkeleton />;
  }

  if ((!isDemo && result.isError) || (!isDemo && !result.data) || !jobId) {
    return (
      <main className="mx-auto max-w-[1500px] px-5 py-6 text-red-300">
        {(result.error as Error | undefined)?.message ?? 'Result is unavailable.'}
      </main>
    );
  }

  const data = isDemo ? demoData : result.data!;

  return (
    <main className="surface-grid mx-auto grid max-w-[1540px] gap-5 px-5 py-6">
      <section className="glass-panel flex flex-wrap items-start justify-between gap-4 rounded-lg border border-line p-4 shadow-tool">
        <div>
          <h1 className="text-2xl font-semibold tracking-normal">{data.metadata.sample}</h1>
          <p className="mt-1 text-sm text-muted">
            {data.metadata.reference_accession ?? 'custom reference'} | {data.metadata.reference_length} bp | engine{' '}
            {data.metadata.engine_version}
          </p>
          <div className="mt-3 flex flex-wrap gap-2">
            {(data.metadata.algorithm_notes ?? []).slice(0, 3).map((note) => (
              <span key={note} className="inline-flex items-center gap-1 rounded-md border border-line bg-panel2 px-2.5 py-1 text-xs text-muted">
                <Binary className="h-3 w-3 text-aqua" aria-hidden />
                {note}
              </span>
            ))}
          </div>
        </div>
        <ExportButtons jobId={jobId} jsonData={data} htmlAvailable={!isDemo} />
      </section>

      <section className="grid gap-3 md:grid-cols-4">
        <Metric icon={<Dna className="h-4 w-4" />} label="Passed reads" value={data.filter_stats.passed_reads} />
        <Metric icon={<Filter className="h-4 w-4" />} label="NUMT filtered" value={data.filter_stats.numt_filtered_reads} />
        <Metric icon={<Activity className="h-4 w-4" />} label="SVs" value={data.svs.length} />
        <Metric icon={<Gauge className="h-4 w-4" />} label="Clusters" value={data.clusters.length} />
      </section>

      <div className="grid gap-5 xl:grid-cols-[minmax(0,1fr)_340px]">
        <section className="glass-panel rounded-lg border border-line p-3 shadow-tool">
          <CircosWrapper
            data={data}
            activeLayers={activeLayers}
            onClusterSelect={handleClusterSelect}
            onSvSelect={handleSvSelect}
          />
        </section>
        <Sidebar genes={data.genes} clusters={data.clusters} />
      </div>

      <ClusterMutationPanel data={data} />
      <AnalysisModulesPanel data={data} />
      <VariantTable data={data} />
    </main>
  );
}

function Metric({ icon, label, value }: { icon: JSX.Element; label: string; value: number }) {
  return (
    <div className="glass-panel rounded-lg border border-line p-4 shadow-tool">
      <div className="flex items-center gap-2 text-sm text-muted">
        <span className="text-aqua">{icon}</span>
        {label}
      </div>
      <div className="mt-2 text-3xl font-semibold tracking-normal">{value}</div>
    </div>
  );
}

function ResultSkeleton() {
  return (
    <main className="surface-grid mx-auto grid max-w-[1540px] gap-5 px-5 py-6">
      <div className="glass-panel h-28 animate-pulse rounded-lg border border-line" />
      <section className="grid gap-3 md:grid-cols-4">
        {Array.from({ length: 4 }).map((_, index) => (
          <div key={index} className="glass-panel h-28 animate-pulse rounded-lg border border-line" />
        ))}
      </section>
      <div className="grid gap-5 xl:grid-cols-[minmax(0,1fr)_340px]">
        <div className="glass-panel h-[560px] animate-pulse rounded-lg border border-line" />
        <div className="grid gap-4">
          <div className="glass-panel h-44 animate-pulse rounded-lg border border-line" />
          <div className="glass-panel h-56 animate-pulse rounded-lg border border-line" />
        </div>
      </div>
    </main>
  );
}
