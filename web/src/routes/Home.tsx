import { Activity, DatabaseZap, FileStack, PlayCircle, ShieldCheck } from 'lucide-react';
import { useNavigate } from 'react-router-dom';
import UploadZone from '../components/UploadZone';
import { demoData } from '../lib/demoData';
import { useMitoStore } from '../lib/store';

export default function Home() {
  const navigate = useNavigate();
  const setSelectedFile = useMitoStore((state) => state.setSelectedFile);
  const setData = useMitoStore((state) => state.setData);
  const setJobId = useMitoStore((state) => state.setJobId);

  return (
    <main className="surface-grid mx-auto grid max-w-[1540px] gap-5 px-5 py-6">
      <section className="grid min-h-[calc(100vh-108px)] gap-5 xl:grid-cols-[minmax(0,1fr)_360px]">
        <UploadZone
          onFile={(file) => {
            setSelectedFile(file);
            navigate('/upload');
          }}
        />
        <aside className="grid content-start gap-4">
          <button
            type="button"
            onClick={() => {
              setData(demoData);
              setJobId('demo');
              setSelectedFile(undefined);
              navigate('/result/demo');
            }}
            className="glass-panel group rounded-lg border border-line p-4 text-left shadow-tool hover:border-aqua hover:shadow-focus"
          >
            <span className="flex items-center justify-between gap-3">
              <span>
                <span className="block text-sm font-semibold">Demo result</span>
                <span className="mt-1 block text-xs text-muted">rCRS, SV, clinical, protein mapping</span>
              </span>
              <span className="grid h-10 w-10 place-items-center rounded-md border border-aqua/40 bg-aqua/10 text-aqua">
                <PlayCircle className="h-5 w-5" aria-hidden />
              </span>
            </span>
          </button>

          <div className="glass-panel rounded-lg border border-line p-4">
            <div className="mb-3 text-xs font-semibold uppercase tracking-normal text-muted">Pipeline</div>
            <div className="grid gap-2">
              <SignalRow icon={<FileStack className="h-4 w-4" />} label="BAM/CRAM/SAM" status="htslib path" />
              <SignalRow icon={<ShieldCheck className="h-4 w-4" />} label="NUMT filter" status="tracked" />
              <SignalRow icon={<Activity className="h-4 w-4" />} label="SV + soft clips" status="enabled" />
              <SignalRow icon={<DatabaseZap className="h-4 w-4" />} label="ClinVar cache" status="local/live-ready" />
            </div>
          </div>

          <div className="glass-panel rounded-lg border border-line p-4">
            <div className="mb-3 text-xs font-semibold uppercase tracking-normal text-muted">Accepted input</div>
            <div className="grid grid-cols-5 gap-2 text-center text-xs font-semibold">
              {['FASTQ', 'SAM', 'BAM', 'CRAM', 'FQ'].map((format) => (
                <span key={format} className="rounded-md border border-line bg-panel2 px-2 py-2 text-muted">
                  {format}
                </span>
              ))}
            </div>
          </div>
        </aside>
      </section>
    </main>
  );
}

function SignalRow({ icon, label, status }: { icon: JSX.Element; label: string; status: string }) {
  return (
    <div className="flex items-center justify-between gap-3 rounded-md border border-line bg-panel2 px-3 py-2">
      <span className="flex items-center gap-2 text-sm">
        <span className="text-aqua">{icon}</span>
        {label}
      </span>
      <span className="text-xs text-muted">{status}</span>
    </div>
  );
}
