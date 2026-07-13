import type { MitoAnalysisData } from '@mito-architect/visualization-lib';

const API_BASE = import.meta.env.VITE_API_URL ?? 'http://127.0.0.1:8080';

export interface UploadResponse {
  job_id: string;
  status: JobStatusValue;
}

export type JobStatusValue = 'queued' | 'processing' | 'done' | 'error' | 'cancelled';

export interface FailureBody {
  schema_version: string;
  code: string;
  message: string;
  retryable: boolean;
}

export interface StatusResponse {
  job_id: string;
  status: JobStatusValue;
  progress: number;
  error?: FailureBody | null;
}

export class ApiClientError extends Error {
  readonly code: string;
  readonly retryable: boolean;
  readonly status: number;

  constructor(failure: FailureBody, status: number) {
    super(`[${failure.code}] ${failure.message}`);
    this.name = 'ApiClientError';
    this.code = failure.code;
    this.retryable = failure.retryable;
    this.status = status;
  }
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
    const error = extractFailure(payload);
    if (error) {
      throw new ApiClientError(error, response.status);
    }
    throw new Error(response.statusText || `HTTP ${response.status}`);
  }
  return payload as T;
}

function extractFailure(payload: unknown): FailureBody | null {
  if (typeof payload !== 'object' || payload === null || !('error' in payload)) {
    return null;
  }
  const error = (payload as { error?: unknown }).error;
  if (typeof error !== 'object' || error === null) {
    return null;
  }
  const candidate = error as Partial<FailureBody>;
  if (
    typeof candidate.schema_version !== 'string' ||
    typeof candidate.code !== 'string' ||
    typeof candidate.message !== 'string' ||
    typeof candidate.retryable !== 'boolean'
  ) {
    return null;
  }
  return candidate as FailureBody;
}
