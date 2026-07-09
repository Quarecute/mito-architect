import type { JobStatusValue } from '../lib/api';

interface ProgressBarProps {
  status?: JobStatusValue;
  progress?: number;
  error?: string | null;
}

export default function ProgressBar({ status = 'queued', progress = 0, error }: ProgressBarProps) {
  const stages = ['queued', 'processing', 'done'] as const;
  return (
    <div className="rounded-lg border border-line bg-panel p-5">
      <div className="flex items-center justify-between gap-4">
        <div>
          <h2 className="text-lg font-semibold tracking-normal">Analysis progress</h2>
          <p className="mt-1 text-sm text-muted">{error ?? statusLabel(status)}</p>
        </div>
        <div className="text-2xl font-semibold">{Math.round(progress)}%</div>
      </div>
      <div className="mt-5 h-2 overflow-hidden rounded-full bg-panel2">
        <div
          className="h-full bg-aqua transition-all duration-300"
          style={{ width: `${Math.max(3, Math.min(100, progress))}%` }}
        />
      </div>
      <div className="mt-5 grid grid-cols-3 gap-2">
        {stages.map((stage) => (
          <div
            key={stage}
            className={[
              'rounded-md border px-3 py-2 text-center text-xs font-semibold uppercase tracking-normal',
              stage === status || (status === 'done' && stage !== 'queued')
                ? 'border-aqua text-aqua'
                : 'border-line text-muted'
            ].join(' ')}
          >
            {stage}
          </div>
        ))}
      </div>
    </div>
  );
}

function statusLabel(status: JobStatusValue): string {
  switch (status) {
    case 'queued':
      return 'Queued for background processing.';
    case 'processing':
      return 'Extracting features, calling SVs, and clustering reads.';
    case 'done':
      return 'Result is ready.';
    case 'error':
      return 'Analysis failed.';
    case 'cancelled':
      return 'Analysis was cancelled.';
  }
}
