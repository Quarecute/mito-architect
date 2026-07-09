import { FileUp, HardDriveUpload, Upload } from 'lucide-react';
import { DragEvent, useRef, useState } from 'react';

const ACCEPTED = ['.fastq', '.fq', '.sam', '.bam', '.cram'];

interface UploadZoneProps {
  onFile: (file: File) => void;
}

export default function UploadZone({ onFile }: UploadZoneProps) {
  const inputRef = useRef<HTMLInputElement>(null);
  const [dragging, setDragging] = useState(false);
  const [error, setError] = useState<string>();

  function handleFile(file?: File) {
    if (!file) return;
    const lower = file.name.toLowerCase();
    if (!ACCEPTED.some((extension) => lower.endsWith(extension))) {
      setError('Expected FASTQ, SAM, BAM, or CRAM input.');
      return;
    }
    setError(undefined);
    onFile(file);
  }

  function onDrop(event: DragEvent<HTMLDivElement>) {
    event.preventDefault();
    setDragging(false);
    handleFile(event.dataTransfer.files.item(0) ?? undefined);
  }

  return (
    <div
      onDragOver={(event) => {
        event.preventDefault();
        setDragging(true);
      }}
      onDragLeave={() => setDragging(false)}
      onDrop={onDrop}
      className={[
        'glass-panel grid min-h-[560px] place-items-center rounded-lg border border-dashed px-6 py-10 text-center shadow-tool transition',
        dragging ? 'border-aqua bg-aqua/10 shadow-focus' : 'border-line'
      ].join(' ')}
    >
      <div className="max-w-xl">
        <div className="mx-auto mb-6 grid h-16 w-16 place-items-center rounded-md border border-aqua/40 bg-aqua/10 shadow-focus">
          {dragging ? <HardDriveUpload className="h-8 w-8 text-aqua" aria-hidden /> : <FileUp className="h-8 w-8 text-aqua" aria-hidden />}
        </div>
        <h1 className="text-4xl font-semibold tracking-normal">Analyze mtDNA reads</h1>
        <p className="mt-3 text-sm leading-6 text-muted">Drop an alignment or choose a file from disk.</p>
        <div className="mt-7 flex flex-wrap justify-center gap-3">
          <button
            type="button"
            onClick={() => inputRef.current?.click()}
            className="inline-flex items-center gap-2 rounded-md border border-aqua bg-aqua px-5 py-3 text-sm font-semibold text-slate-950 shadow-focus hover:bg-teal-200"
          >
            <Upload className="h-4 w-4" aria-hidden />
            Choose file
          </button>
          <div className="rounded-md border border-line bg-panel2 px-4 py-3 text-sm text-muted">
            {ACCEPTED.join(', ')}
          </div>
        </div>
        {error && <p className="mt-4 text-sm text-red-300">{error}</p>}
        <input
          ref={inputRef}
          type="file"
          className="hidden"
          accept={ACCEPTED.join(',')}
          onChange={(event) => handleFile(event.target.files?.item(0) ?? undefined)}
        />
      </div>
    </div>
  );
}
