import { useMutation, useQuery } from '@tanstack/react-query';
import { useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import ProgressBar from '../components/ProgressBar';
import UploadZone from '../components/UploadZone';
import { getStatus, uploadFile } from '../lib/api';
import { useMitoStore } from '../lib/store';

export default function UploadPage() {
  const navigate = useNavigate();
  const selectedFile = useMitoStore((state) => state.selectedFile);
  const setSelectedFile = useMitoStore((state) => state.setSelectedFile);
  const jobId = useMitoStore((state) => state.jobId);
  const setJobId = useMitoStore((state) => state.setJobId);

  const upload = useMutation({
    mutationFn: uploadFile,
    onSuccess: (response) => setJobId(response.job_id)
  });

  useEffect(() => {
    if (selectedFile && !jobId && !upload.isPending && !upload.isSuccess && !upload.isError) {
      upload.mutate(selectedFile);
    }
  }, [jobId, selectedFile, upload]);

  const status = useQuery({
    queryKey: ['status', jobId],
    queryFn: ({ signal }) => getStatus(jobId!, signal),
    enabled: Boolean(jobId),
    refetchInterval: (query) =>
      query.state.data && ['done', 'error', 'cancelled'].includes(query.state.data.status)
        ? false
        : 1000
  });

  useEffect(() => {
    if (status.data?.status === 'done') {
      navigate(`/result/${status.data.job_id}`);
    }
  }, [navigate, status.data]);

  if (!selectedFile && !jobId) {
    return (
      <main className="mx-auto max-w-[1500px] px-5 py-6">
        <UploadZone
          onFile={(file) => {
            setSelectedFile(file);
            setJobId(undefined);
          }}
        />
      </main>
    );
  }

  return (
    <main className="mx-auto grid max-w-[1100px] gap-5 px-5 py-6">
      <section className="rounded-lg border border-line bg-panel p-5">
        <h1 className="text-2xl font-semibold tracking-normal">{selectedFile?.name ?? 'Upload'}</h1>
        <p className="mt-2 text-sm text-muted">
          The backend stores the file in a per-job temp directory and runs the native engine in a
          blocking worker.
        </p>
      </section>
      <ProgressBar
        status={status.data?.status ?? (upload.isError ? 'error' : upload.isPending ? 'queued' : 'processing')}
        progress={status.data?.progress ?? (upload.isPending ? 5 : 0)}
        error={(upload.error as Error | undefined)?.message ?? status.data?.error}
      />
    </main>
  );
}
