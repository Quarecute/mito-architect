import type { MitoAnalysisData } from '@mito-architect/visualization-lib';

const API_BASE = import.meta.env.VITE_API_URL ?? 'http://127.0.0.1:8080';

export interface UploadResponse {
  job_id: string;
  status: JobStatusValue;
}

export type JobStatusValue = 'queued' | 'processing' | 'done' | 'error' | 'cancelled';

export interface StatusResponse {
  job_id: string;
  status: JobStatusValue;
  progress: number;
  error?: string | null;
}

export async function uploadFile(file: File): Promise<UploadResponse> {
  const body = new FormData();
  body.append('file', file);
  const response = await fetch(`${API_BASE}/upload`, {
    method: 'POST',
    body
  });
  return parseResponse<UploadResponse>(response);
}

export async function getStatus(jobId: string, signal?: AbortSignal): Promise<StatusResponse> {
  const response = await fetch(`${API_BASE}/status/${encodeURIComponent(jobId)}`, { signal });
  return parseResponse<StatusResponse>(response);
}

export async function getResult(jobId: string, signal?: AbortSignal): Promise<MitoAnalysisData> {
  const response = await fetch(`${API_BASE}/result/${encodeURIComponent(jobId)}`, { signal });
  return parseResponse<MitoAnalysisData>(response);
}

export function downloadReportUrl(jobId: string): string {
  return `${API_BASE}/download/${encodeURIComponent(jobId)}`;
}

async function parseResponse<T>(response: Response): Promise<T> {
  const text = await response.text();
  let payload: unknown = {};
  try {
    payload = text ? JSON.parse(text) : {};
  } catch {
    payload = { error: text || response.statusText };
  }
  if (!response.ok) {
    const message =
      typeof payload === 'object' && payload !== null && 'error' in payload
        ? String((payload as { error?: unknown }).error)
        : response.statusText;
    throw new Error(message);
  }
  return payload as T;
}
