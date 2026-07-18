import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { Download, FileJson, Table } from 'lucide-react';
import { downloadReportUrl, downloadResultJsonUrl } from '../lib/api';

interface ExportButtonsProps {
  jobId: string;
  jsonData: MitoAnalysisData;
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

  function downloadVariantTsv() {
    const header = [
      'event_id', 'type', 'position', 'ref', 'alt', 'alt_molecules', 'ref_molecules',
      'other_molecules', 'callable_molecules', 'heteroplasmy', 'ci95_low', 'ci95_high',
      'numt_assessability', 'qc_flags', 'supporting_molecule_ids'
    ];
    const rows = (jsonData.variants ?? []).map((variant) => [
      variant.event_id ?? `legacy:${variant.position}:${variant.ref}:${variant.alt}`,
      variant.type ?? 'SNV',
      variant.position,
      variant.ref,
      variant.alt,
      variant.alt_depth,
      variant.ref_depth,
      variant.other_depth,
      variant.callable_depth,
      variant.heteroplasmy.toFixed(8),
      variant.ci95_low.toFixed(8),
      variant.ci95_high.toFixed(8),
      variant.numt_assessability ?? '',
      (variant.qc_flags ?? []).join(';'),
      (variant.supporting_molecule_ids ?? variant.supporting_reads).join(';')
    ].map(tsvField).join('\t'));
    downloadBlob(
      `${header.join('\t')}\n${rows.join('\n')}\n`,
      'text/tab-separated-values;charset=utf-8',
      `mito-${jobId}-variants.tsv`
    );
  }

  return (
    <div className="flex flex-wrap gap-2">
      {htmlAvailable ? (
        <a
          href={downloadResultJsonUrl(jobId)}
          download={`mito-${jobId}.json`}
          className="inline-flex items-center gap-2 rounded-md border border-line bg-panel2 px-3 py-2 text-sm font-semibold hover:border-aqua focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-aqua"
        >
          <FileJson className="h-4 w-4" aria-hidden />
          Full JSON
        </a>
      ) : (
        <button
          type="button"
          onClick={downloadJson}
          className="inline-flex items-center gap-2 rounded-md border border-line bg-panel2 px-3 py-2 text-sm font-semibold hover:border-aqua focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-aqua"
        >
          <FileJson className="h-4 w-4" aria-hidden />
          JSON
        </button>
      )}
      <button
        type="button"
        onClick={downloadVariantTsv}
        className="inline-flex items-center gap-2 rounded-md border border-line bg-panel2 px-3 py-2 text-sm font-semibold hover:border-aqua focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-aqua"
      >
        <Table className="h-4 w-4" aria-hidden />
        Variant TSV
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

function tsvField(value: string | number): string {
  return String(value).replace(/[\t\r\n]+/g, ' ').trim();
}

function downloadBlob(contents: string, type: string, fileName: string) {
  const blob = new Blob([contents], { type });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = fileName;
  link.click();
  URL.revokeObjectURL(url);
}
