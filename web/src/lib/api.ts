import type { EvidenceObservation, MitoAnalysisData } from '@mito-architect/visualization-lib';

const API_BASE = import.meta.env.VITE_API_URL ?? 'http://127.0.0.1:8080';

export interface UploadResponse {
  job_id: string;
  status: JobStatusValue;
  result_schema: string;
}

export interface AnalysisUploadOptions {
  evidencePageSize?: number;
  moleculeIdTag?: string;
  umiTag?: string;
  duplexTag?: string;
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

export interface EvidenceSearchFilters {
  moleculeId?: string;
  eventId?: string;
  state?: string;
}

export interface EvidenceSearchRow {
  id: string;
  page_index: number;
  row_index: number;
  molecule_id: string;
  event_id: string;
  alignment_id: string;
  state: EvidenceObservation['state'];
  observed_allele: string | null;
  base_quality: number | null;
  mapping_quality: number;
  strand: string;
  evidence_source: string;
  read_position: number | null;
}

export interface EvidenceSearchResponse {
  schema_version: '1.0';
  filters: {
    molecule_id: string | null;
    event_id: string | null;
    state: string | null;
  };
  cursor: number;
  next_cursor: number | null;
  total_matches: number;
  rows: EvidenceSearchRow[];
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
  return uploadFileWithOptions(file, {});
}

export async function uploadFileWithOptions(
  file: File,
  options: AnalysisUploadOptions
): Promise<UploadResponse> {
  const body = new FormData();
  body.append('evidence_graph', 'true');
  if (options.evidencePageSize !== undefined) {
    body.append('evidence_page_size', options.evidencePageSize.toString());
  }
  if (options.moleculeIdTag) body.append('molecule_id_tag', options.moleculeIdTag);
  if (options.umiTag) body.append('umi_tag', options.umiTag);
  if (options.duplexTag) body.append('duplex_tag', options.duplexTag);
  body.append('file', file);
  const response = await fetch(`${API_BASE}/upload`, {
    method: 'POST',
    body
  });
  return parseResponse<UploadResponse>(response);
}

export async function getEvidencePage(jobId: string, pageIndex: number, signal?: AbortSignal) {
  const response = await fetch(
    `${API_BASE}/result/${encodeURIComponent(jobId)}/evidence/${pageIndex}`,
    { signal }
  );
  return parseResponse<NonNullable<MitoAnalysisData['observation_pages']>[number]>(response);
}

export async function searchEvidence(
  jobId: string,
  filters: EvidenceSearchFilters,
  cursor: number,
  limit: number,
  signal?: AbortSignal
): Promise<EvidenceSearchResponse> {
  const query = new URLSearchParams({ cursor: cursor.toString(), limit: limit.toString() });
  if (filters.moleculeId) query.set('molecule_id', filters.moleculeId);
  if (filters.eventId) query.set('event_id', filters.eventId);
  if (filters.state) query.set('state', filters.state);
  const response = await fetch(
    `${API_BASE}/result/${encodeURIComponent(jobId)}/evidence?${query.toString()}`,
    { signal }
  );
  return parseResponse<EvidenceSearchResponse>(response);
}

export async function getStatus(jobId: string, signal?: AbortSignal): Promise<StatusResponse> {
  const response = await fetch(`${API_BASE}/status/${encodeURIComponent(jobId)}`, { signal });
  return parseResponse<StatusResponse>(response);
}

export async function getResult(jobId: string, signal?: AbortSignal): Promise<MitoAnalysisData> {
  const response = await fetch(`${API_BASE}/result/${encodeURIComponent(jobId)}/summary`, { signal });
  return parseResponse<MitoAnalysisData>(response);
}

export function downloadResultJsonUrl(jobId: string): string {
  return `${API_BASE}/result/${encodeURIComponent(jobId)}`;
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
