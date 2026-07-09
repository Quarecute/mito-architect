import { Download, FileJson } from 'lucide-react';
import { downloadReportUrl } from '../lib/api';

interface ExportButtonsProps {
  jobId: string;
  jsonData: unknown;
  htmlAvailable?: boolean;
}

export default function ExportButtons({ jobId, jsonData, htmlAvailable = true }: ExportButtonsProps) {
  function downloadJson() {
    const blob = new Blob([JSON.stringify(jsonData, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `mito-${jobId}.json`;
    link.click();
    URL.revokeObjectURL(url);
  }

  return (
    <div className="flex flex-wrap gap-2">
      <button
        type="button"
        onClick={downloadJson}
        className="inline-flex items-center gap-2 rounded-md border border-line bg-panel2 px-3 py-2 text-sm font-semibold hover:border-aqua"
      >
        <FileJson className="h-4 w-4" aria-hidden />
        JSON
      </button>
      {htmlAvailable && (
        <a
          href={downloadReportUrl(jobId)}
          className="inline-flex items-center gap-2 rounded-md border border-line bg-panel2 px-3 py-2 text-sm font-semibold hover:border-aqua"
        >
          <Download className="h-4 w-4" aria-hidden />
          HTML
        </a>
      )}
    </div>
  );
}
