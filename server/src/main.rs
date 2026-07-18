use axum::body::{Body, Bytes};
use axum::extract::{DefaultBodyLimit, Multipart, Path as AxumPath, Query, Request, State};
use axum::http::{header, HeaderValue, Method, StatusCode};
use axum::middleware::{self, Next};
use axum::response::{IntoResponse, Response};
use axum::routing::{get, post};
use axum::{Json, Router};
use mito_ffi::{AnalyzeOptions, MitoEngine};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{HashMap, HashSet};
use std::net::SocketAddr;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::{Duration, SystemTime, UNIX_EPOCH};
use tokio::fs;
use tokio::io::AsyncWriteExt;
use tokio::sync::{OwnedSemaphorePermit, RwLock, Semaphore};
use tower_http::cors::CorsLayer;
use uuid::Uuid;

const DEFAULT_MAX_UPLOAD_BYTES: usize = 64 * 1024 * 1024;
const JOB_TTL_SECS: u64 = 24 * 60 * 60;

#[derive(Clone)]
struct AppState {
    jobs: Arc<RwLock<HashMap<Uuid, JobRecord>>>,
    cancel_flags: Arc<RwLock<HashMap<Uuid, Arc<AtomicBool>>>>,
    analysis_slots: Arc<Semaphore>,
    job_slots: Arc<Semaphore>,
    worker_threads: usize,
    min_mapping_quality: u8,
    min_base_quality: u8,
    excluded_snp_flags: u16,
    numt_threshold: f64,
    max_evidence_observations: usize,
    max_phase_links: usize,
    evidence_page_size: usize,
    engine_version: Arc<str>,
    schema_version: Arc<str>,
    error_schema_version: Arc<str>,
    htslib_enabled: bool,
    api_key: Option<Arc<str>>,
    max_upload_bytes: usize,
    tmp_dir: PathBuf,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
enum JobStatus {
    Queued,
    Processing,
    Done,
    Error,
    Cancelled,
}

#[derive(Debug)]
struct JobRecord {
    status: JobStatus,
    progress: u8,
    input_path: PathBuf,
    result: Option<Bytes>,
    result_summary: Option<Bytes>,
    evidence_pages: Vec<Bytes>,
    evidence_search_index: Option<Arc<EvidenceSearchIndex>>,
    html_report: Option<Bytes>,
    error: Option<FailureBody>,
    cancel_requested: bool,
    created_at: u64,
}

#[derive(Debug)]
struct CompletedArtifacts {
    result: Bytes,
    result_summary: Bytes,
    evidence_pages: Vec<Bytes>,
    evidence_search_index: Option<Arc<EvidenceSearchIndex>>,
    html_report: Bytes,
}

#[derive(Debug)]
struct EvidenceSearchIndex {
    strings: Vec<Arc<str>>,
    string_ids: HashMap<Arc<str>, u32>,
    rows: Vec<EvidenceSearchRecord>,
    by_molecule: HashMap<u32, Vec<usize>>,
    by_event: HashMap<u32, Vec<usize>>,
    by_state: HashMap<u32, Vec<usize>>,
}

#[derive(Clone, Copy, Debug)]
struct EvidenceSearchRecord {
    global_index: usize,
    page_index: usize,
    row_index: usize,
    molecule_id: u32,
    event_id: u32,
    alignment_id: u32,
    state: u32,
    observed_allele: Option<u32>,
    base_quality: Option<u64>,
    mapping_quality: u64,
    strand: u32,
    evidence_source: u32,
    read_position: Option<f64>,
}

#[derive(Debug, Deserialize)]
struct EvidenceSearchQuery {
    molecule_id: Option<String>,
    event_id: Option<String>,
    state: Option<String>,
    cursor: Option<usize>,
    limit: Option<usize>,
}

#[derive(Debug, Serialize)]
struct EvidenceSearchFilters {
    molecule_id: Option<String>,
    event_id: Option<String>,
    state: Option<String>,
}

#[derive(Debug, Serialize)]
struct EvidenceSearchResponse {
    schema_version: &'static str,
    filters: EvidenceSearchFilters,
    cursor: usize,
    next_cursor: Option<usize>,
    total_matches: usize,
    rows: Vec<EvidenceSearchResult>,
}

#[derive(Debug, Serialize)]
struct EvidenceSearchResult {
    id: String,
    page_index: usize,
    row_index: usize,
    molecule_id: String,
    event_id: String,
    alignment_id: String,
    state: String,
    observed_allele: Option<String>,
    base_quality: Option<u64>,
    mapping_quality: u64,
    strand: String,
    evidence_source: String,
    read_position: Option<f64>,
}

#[derive(Debug, Serialize)]
struct UploadResponse {
    job_id: Uuid,
    status: JobStatus,
    result_schema: &'static str,
}

#[derive(Clone, Debug)]
struct JobAnalysisOptions {
    emit_evidence_graph: bool,
    evidence_page_size: usize,
    molecule_id_tag: String,
    umi_tag: String,
    duplex_tag: String,
}

#[derive(Debug, Serialize)]
struct StatusResponse {
    job_id: Uuid,
    status: JobStatus,
    progress: u8,
    error: Option<FailureBody>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct FailureBody {
    schema_version: &'static str,
    code: String,
    message: String,
    retryable: bool,
}

impl FailureBody {
    fn new(code: impl Into<String>, message: impl Into<String>, retryable: bool) -> Self {
        Self {
            schema_version: "1.0",
            code: code.into(),
            message: message.into(),
            retryable,
        }
    }
}

#[derive(Debug)]
struct ApiError {
    status: StatusCode,
    body: FailureBody,
}

impl ApiError {
    fn new(status: StatusCode, message: impl Into<String>) -> Self {
        let (code, retryable) = match status {
            StatusCode::BAD_REQUEST => ("MITO-API-E1001", false),
            StatusCode::UNAUTHORIZED => ("MITO-API-E1002", false),
            StatusCode::NOT_FOUND => ("MITO-API-E1003", false),
            StatusCode::CONFLICT => ("MITO-API-E1004", false),
            StatusCode::PAYLOAD_TOO_LARGE => ("MITO-API-E1005", false),
            StatusCode::TOO_MANY_REQUESTS => ("MITO-API-E1006", true),
            StatusCode::ACCEPTED => ("MITO-API-E1007", true),
            StatusCode::SERVICE_UNAVAILABLE => ("MITO-API-E2001", true),
            _ => ("MITO-API-E9001", false),
        };
        Self::coded(status, code, message, retryable)
    }

    fn coded(
        status: StatusCode,
        code: impl Into<String>,
        message: impl Into<String>,
        retryable: bool,
    ) -> Self {
        Self {
            status,
            body: FailureBody::new(code, message, retryable),
        }
    }
}

impl IntoResponse for ApiError {
    fn into_response(self) -> Response {
        (self.status, Json(serde_json::json!({ "error": self.body }))).into_response()
    }
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let tmp_dir = PathBuf::from("./tmp");
    fs::create_dir_all(&tmp_dir).await?;

    let available_threads = std::thread::available_parallelism().map_or(1, usize::from);
    let worker_threads = positive_env_usize("MITO_JOB_THREADS")?.unwrap_or(available_threads);
    let max_concurrent_jobs = positive_env_usize("MITO_MAX_CONCURRENT_JOBS")?.unwrap_or(1);
    let max_active_jobs = positive_env_usize("MITO_MAX_ACTIVE_JOBS")?.unwrap_or(16);
    let max_upload_bytes =
        positive_env_usize("MITO_MAX_UPLOAD_BYTES")?.unwrap_or(DEFAULT_MAX_UPLOAD_BYTES);
    let max_evidence_observations =
        positive_env_usize("MITO_MAX_EVIDENCE_OBSERVATIONS")?.unwrap_or(5_000_000);
    let max_phase_links = positive_env_usize("MITO_MAX_PHASE_LINKS")?.unwrap_or(1_000_000);
    let evidence_page_size = positive_env_usize("MITO_EVIDENCE_PAGE_SIZE")?.unwrap_or(4096);
    if evidence_page_size > 1_000_000 {
        anyhow::bail!("MITO_EVIDENCE_PAGE_SIZE must not exceed 1000000");
    }
    let min_mapping_quality = env_u8("MITO_MIN_MAPQ")?.unwrap_or(20);
    let min_base_quality = env_u8("MITO_MIN_BASE_QUALITY")?.unwrap_or(10);
    let excluded_snp_flags = env_u16("MITO_EXCLUDED_SNP_FLAGS")?.unwrap_or(0xF00);
    let numt_threshold = optional_env_number("MITO_NUMT_THRESHOLD")?.unwrap_or(0.30);
    if !(0.0..=1.0).contains(&numt_threshold) {
        anyhow::bail!("MITO_NUMT_THRESHOLD must be between 0 and 1");
    }
    let capabilities = MitoEngine::capabilities();
    if !capabilities.htslib {
        anyhow::bail!("mito-server requires a native core built with htslib for BAM/CRAM support");
    }
    let api_key = std::env::var("MITO_API_KEY")
        .ok()
        .filter(|value| !value.is_empty())
        .map(Arc::<str>::from);
    let addr = std::env::var("MITO_SERVER_ADDR")
        .ok()
        .and_then(|value| value.parse::<SocketAddr>().ok())
        .unwrap_or_else(|| SocketAddr::from(([127, 0, 0, 1], 8080)));
    if !addr.ip().is_loopback() && api_key.is_none() {
        anyhow::bail!("MITO_API_KEY is required when MITO_SERVER_ADDR is not loopback");
    }
    let state = AppState {
        jobs: Arc::new(RwLock::new(HashMap::new())),
        cancel_flags: Arc::new(RwLock::new(HashMap::new())),
        analysis_slots: Arc::new(Semaphore::new(max_concurrent_jobs)),
        job_slots: Arc::new(Semaphore::new(max_active_jobs)),
        worker_threads,
        min_mapping_quality,
        min_base_quality,
        excluded_snp_flags,
        numt_threshold,
        max_evidence_observations,
        max_phase_links,
        evidence_page_size,
        engine_version: Arc::from(capabilities.engine_version),
        schema_version: Arc::from(capabilities.schema_version),
        error_schema_version: Arc::from(capabilities.error_schema_version),
        htslib_enabled: capabilities.htslib,
        api_key,
        max_upload_bytes,
        tmp_dir,
    };

    spawn_cleanup(state.clone());
    let cors = cors_layer()?;
    let app = build_app(state, cors);

    let listener = tokio::net::TcpListener::bind(addr).await?;
    println!("mito-server listening on http://{addr}");
    axum::serve(listener, app).await?;
    Ok(())
}

fn build_app(state: AppState, cors: CorsLayer) -> Router {
    let max_upload_bytes = state.max_upload_bytes;
    let protected = Router::new()
        .route("/upload", post(upload))
        .route("/status/:job_id", get(status))
        .route("/result/:job_id", get(result))
        .route("/result/:job_id/summary", get(result_summary))
        .route("/result/:job_id/evidence", get(search_evidence))
        .route("/result/:job_id/evidence/:page_index", get(evidence_page))
        .route("/download/:job_id", get(download))
        .route("/cancel/:job_id", post(cancel))
        .route_layer(middleware::from_fn_with_state(
            state.clone(),
            require_api_key,
        ));
    Router::new()
        .route("/healthz", get(health))
        .route("/readyz", get(ready))
        .merge(protected)
        .layer(DefaultBodyLimit::max(max_upload_bytes))
        .layer(cors)
        .with_state(state)
}

async fn health() -> Json<Value> {
    Json(serde_json::json!({ "status": "ok" }))
}

async fn ready(State(state): State<AppState>) -> Result<Json<Value>, ApiError> {
    fs::metadata(&state.tmp_dir)
        .await
        .map_err(|error| ApiError::new(StatusCode::SERVICE_UNAVAILABLE, error.to_string()))?;
    Ok(Json(serde_json::json!({
        "status": "ready",
        "engine_version": state.engine_version.as_ref(),
        "schema_version": state.schema_version.as_ref(),
        "error_schema_version": state.error_schema_version.as_ref(),
        "htslib_enabled": state.htslib_enabled,
        "available_analysis_slots": state.analysis_slots.available_permits(),
        "available_job_slots": state.job_slots.available_permits()
    })))
}

async fn require_api_key(
    State(state): State<AppState>,
    request: Request,
    next: Next,
) -> Result<Response, ApiError> {
    let Some(expected) = state.api_key.as_deref() else {
        return Ok(next.run(request).await);
    };
    let supplied = request
        .headers()
        .get(header::AUTHORIZATION)
        .and_then(|value| value.to_str().ok())
        .and_then(|value| value.strip_prefix("Bearer "));
    if !supplied.is_some_and(|value| constant_time_eq(value.as_bytes(), expected.as_bytes())) {
        return Err(ApiError::new(
            StatusCode::UNAUTHORIZED,
            "invalid or missing API key",
        ));
    }
    Ok(next.run(request).await)
}

fn constant_time_eq(left: &[u8], right: &[u8]) -> bool {
    let mut difference = left.len() ^ right.len();
    let maximum = left.len().max(right.len());
    for index in 0..maximum {
        let lhs = left.get(index).copied().unwrap_or(0);
        let rhs = right.get(index).copied().unwrap_or(0);
        difference |= usize::from(lhs ^ rhs);
    }
    difference == 0
}

async fn upload(
    State(state): State<AppState>,
    multipart: Result<Multipart, axum::extract::multipart::MultipartRejection>,
) -> Result<Json<UploadResponse>, ApiError> {
    let mut multipart = multipart.map_err(|error| {
        eprintln!("multipart extraction failed: {error}");
        ApiError::new(StatusCode::BAD_REQUEST, "invalid multipart upload")
    })?;
    let job_permit = state
        .job_slots
        .clone()
        .try_acquire_owned()
        .map_err(|_| ApiError::new(StatusCode::TOO_MANY_REQUESTS, "analysis queue is full"))?;
    let mut analysis_options = JobAnalysisOptions {
        emit_evidence_graph: false,
        evidence_page_size: state.evidence_page_size,
        molecule_id_tag: String::new(),
        umi_tag: String::new(),
        duplex_tag: String::new(),
    };
    let mut uploaded: Option<(Uuid, PathBuf, Arc<AtomicBool>)> = None;
    loop {
        let next_field = match multipart.next_field().await {
            Ok(field) => field,
            Err(error) => {
                if let Some((job_id, _, _)) = &uploaded {
                    let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
                }
                return Err(ApiError::new(StatusCode::BAD_REQUEST, error.to_string()));
            }
        };
        let Some(mut field) = next_field else {
            break;
        };
        let field_name = field.name().unwrap_or_default().to_owned();
        if matches!(
            field_name.as_str(),
            "evidence_graph" | "evidence_page_size" | "molecule_id_tag" | "umi_tag" | "duplex_tag"
        ) {
            let value = match field.text().await {
                Ok(value) => value,
                Err(error) => {
                    if let Some((job_id, _, _)) = &uploaded {
                        let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
                    }
                    return Err(ApiError::new(StatusCode::BAD_REQUEST, error.to_string()));
                }
            };
            match field_name.as_str() {
                "evidence_graph" => {
                    analysis_options.emit_evidence_graph = match value.trim() {
                        "true" | "1" => true,
                        "false" | "0" => false,
                        _ => {
                            if let Some((job_id, _, _)) = &uploaded {
                                let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string()))
                                    .await;
                            }
                            return Err(ApiError::new(
                                StatusCode::BAD_REQUEST,
                                "evidence_graph must be true, false, 1, or 0",
                            ));
                        }
                    };
                }
                "evidence_page_size" => {
                    let parsed = value
                        .trim()
                        .parse::<usize>()
                        .ok()
                        .filter(|value| (1..=1_000_000).contains(value));
                    let Some(parsed) = parsed else {
                        if let Some((job_id, _, _)) = &uploaded {
                            let _ =
                                fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
                        }
                        return Err(ApiError::new(
                            StatusCode::BAD_REQUEST,
                            "evidence_page_size must be between 1 and 1000000",
                        ));
                    };
                    analysis_options.evidence_page_size = parsed;
                }
                "molecule_id_tag" => analysis_options.molecule_id_tag = value.trim().to_owned(),
                "umi_tag" => analysis_options.umi_tag = value.trim().to_owned(),
                "duplex_tag" => analysis_options.duplex_tag = value.trim().to_owned(),
                _ => unreachable!("multipart protocol field was matched above"),
            }
            for tag in [
                &analysis_options.molecule_id_tag,
                &analysis_options.umi_tag,
                &analysis_options.duplex_tag,
            ] {
                if !valid_optional_sam_tag(tag) {
                    if let Some((job_id, _, _)) = &uploaded {
                        let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
                    }
                    return Err(ApiError::new(
                        StatusCode::BAD_REQUEST,
                        "protocol SAM tags must match [A-Za-z][A-Za-z0-9]",
                    ));
                }
            }
            continue;
        }
        if field.name() != Some("file") {
            continue;
        }
        if let Some((job_id, _, _)) = &uploaded {
            let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
            return Err(ApiError::new(
                StatusCode::BAD_REQUEST,
                "multipart upload must include exactly one file field",
            ));
        }

        let original_name = field.file_name().unwrap_or("input.fastq").to_string();
        validate_upload_name(&original_name)?;
        let job_id = Uuid::new_v4();
        let cancel_flag = Arc::new(AtomicBool::new(false));
        let job_dir = state.tmp_dir.join(job_id.to_string());
        fs::create_dir_all(&job_dir).await.map_err(internal_error)?;
        let input_path = job_dir.join(sanitize_file_name(&original_name));
        let mut output = match fs::File::create(&input_path).await {
            Ok(output) => output,
            Err(error) => {
                let _ = fs::remove_dir_all(&job_dir).await;
                return Err(internal_error(error));
            }
        };
        let mut written = 0usize;

        loop {
            let chunk = match field.chunk().await {
                Ok(Some(chunk)) => chunk,
                Ok(None) => break,
                Err(error) => {
                    drop(output);
                    let _ = fs::remove_dir_all(&job_dir).await;
                    return Err(ApiError::new(StatusCode::BAD_REQUEST, error.to_string()));
                }
            };
            written = match written.checked_add(chunk.len()) {
                Some(total) => total,
                None => {
                    drop(output);
                    let _ = fs::remove_dir_all(&job_dir).await;
                    return Err(ApiError::new(
                        StatusCode::PAYLOAD_TOO_LARGE,
                        "upload size overflow",
                    ));
                }
            };
            if written > state.max_upload_bytes {
                drop(output);
                let _ = fs::remove_dir_all(&job_dir).await;
                return Err(ApiError::new(
                    StatusCode::PAYLOAD_TOO_LARGE,
                    format!("upload exceeds the {} byte limit", state.max_upload_bytes),
                ));
            }
            if let Err(error) = output.write_all(&chunk).await {
                drop(output);
                let _ = fs::remove_dir_all(&job_dir).await;
                return Err(internal_error(error));
            }
        }

        if written == 0 {
            drop(output);
            let _ = fs::remove_dir_all(&job_dir).await;
            return Err(ApiError::new(
                StatusCode::BAD_REQUEST,
                "uploaded file is empty",
            ));
        }
        if let Err(error) = output.flush().await {
            drop(output);
            let _ = fs::remove_dir_all(&job_dir).await;
            return Err(internal_error(error));
        }
        drop(output);
        uploaded = Some((job_id, input_path, cancel_flag));
    }

    let Some((job_id, input_path, cancel_flag)) = uploaded else {
        return Err(ApiError::new(
            StatusCode::BAD_REQUEST,
            "multipart upload must include a file field named 'file'",
        ));
    };
    if !analysis_options.emit_evidence_graph
        && (!analysis_options.molecule_id_tag.is_empty()
            || !analysis_options.umi_tag.is_empty()
            || !analysis_options.duplex_tag.is_empty())
    {
        let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
        return Err(ApiError::new(
            StatusCode::BAD_REQUEST,
            "protocol SAM tags require evidence_graph=true",
        ));
    }
    let record = JobRecord {
        status: JobStatus::Queued,
        progress: 0,
        input_path: input_path.clone(),
        result: None,
        result_summary: None,
        evidence_pages: Vec::new(),
        evidence_search_index: None,
        html_report: None,
        error: None,
        cancel_requested: false,
        created_at: now_epoch_secs(),
    };
    state.jobs.write().await.insert(job_id, record);
    state
        .cancel_flags
        .write()
        .await
        .insert(job_id, cancel_flag.clone());

    let job_state = state.clone();
    let job_options = analysis_options.clone();
    tokio::spawn(async move {
        run_job(
            job_state,
            job_id,
            input_path,
            cancel_flag,
            job_options,
            job_permit,
        )
        .await;
    });

    Ok(Json(UploadResponse {
        job_id,
        status: JobStatus::Queued,
        result_schema: if analysis_options.emit_evidence_graph {
            "0.6"
        } else {
            "0.5"
        },
    }))
}

async fn status(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
) -> Result<Json<StatusResponse>, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let jobs = state.jobs.read().await;
    let job = jobs
        .get(&job_id)
        .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;
    Ok(Json(StatusResponse {
        job_id,
        status: job.status.clone(),
        progress: job.progress,
        error: job.error.clone(),
    }))
}

async fn cancel(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
) -> Result<Json<StatusResponse>, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let response = {
        let mut jobs = state.jobs.write().await;
        let job = jobs
            .get_mut(&job_id)
            .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;

        match job.status {
            JobStatus::Done => {
                return Err(ApiError::new(
                    StatusCode::CONFLICT,
                    "completed jobs cannot be cancelled",
                ));
            }
            JobStatus::Error | JobStatus::Cancelled => {}
            JobStatus::Queued | JobStatus::Processing => {
                job.cancel_requested = true;
                job.status = JobStatus::Cancelled;
                job.progress = 100;
                job.error = Some(FailureBody::new("MITO-E1501", "cancelled by user", false));
            }
        }

        StatusResponse {
            job_id,
            status: job.status.clone(),
            progress: job.progress,
            error: job.error.clone(),
        }
    };
    if let Some(flag) = state.cancel_flags.read().await.get(&job_id) {
        flag.store(true, Ordering::Relaxed);
    }

    Ok(Json(response))
}

async fn result(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
) -> Result<Response, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let jobs = state.jobs.read().await;
    let job = jobs
        .get(&job_id)
        .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;

    match (&job.status, &job.result) {
        (JobStatus::Done, Some(value)) => {
            let mut response = Response::new(Body::from(value.clone()));
            response.headers_mut().insert(
                header::CONTENT_TYPE,
                HeaderValue::from_static("application/json; charset=utf-8"),
            );
            response.headers_mut().insert(
                header::CONTENT_DISPOSITION,
                HeaderValue::from_str(&format!("attachment; filename=\"mito-{job_id}.json\""))
                    .map_err(internal_error)?,
            );
            Ok(response)
        }
        (JobStatus::Error, _) => {
            let failure = job
                .error
                .clone()
                .unwrap_or_else(|| FailureBody::new("MITO-E9001", "analysis failed", false));
            Err(ApiError::coded(
                StatusCode::CONFLICT,
                failure.code,
                failure.message,
                failure.retryable,
            ))
        }
        (JobStatus::Cancelled, _) => Err(ApiError::coded(
            StatusCode::CONFLICT,
            "MITO-E1501",
            "job was cancelled",
            false,
        )),
        _ => Err(ApiError::new(
            StatusCode::ACCEPTED,
            "analysis result is not ready yet",
        )),
    }
}

async fn result_summary(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
) -> Result<Response, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let jobs = state.jobs.read().await;
    let job = jobs
        .get(&job_id)
        .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;

    match (&job.status, &job.result_summary) {
        (JobStatus::Done, Some(value)) => {
            let mut response = Response::new(Body::from(value.clone()));
            response.headers_mut().insert(
                header::CONTENT_TYPE,
                HeaderValue::from_static("application/json; charset=utf-8"),
            );
            response.headers_mut().insert(
                header::CACHE_CONTROL,
                HeaderValue::from_static("private, immutable, max-age=86400"),
            );
            Ok(response)
        }
        (JobStatus::Error, _) => {
            let failure = job
                .error
                .clone()
                .unwrap_or_else(|| FailureBody::new("MITO-E9001", "analysis failed", false));
            Err(ApiError::coded(
                StatusCode::CONFLICT,
                failure.code,
                failure.message,
                failure.retryable,
            ))
        }
        (JobStatus::Cancelled, _) => Err(ApiError::coded(
            StatusCode::CONFLICT,
            "MITO-E1501",
            "job was cancelled",
            false,
        )),
        _ => Err(ApiError::new(
            StatusCode::ACCEPTED,
            "analysis result is not ready yet",
        )),
    }
}

async fn search_evidence(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
    Query(query): Query<EvidenceSearchQuery>,
) -> Result<Response, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let limit = query.limit.unwrap_or(100);
    if !(1..=500).contains(&limit) {
        return Err(ApiError::new(
            StatusCode::BAD_REQUEST,
            "evidence search limit must be between 1 and 500",
        ));
    }
    let filters = EvidenceSearchFilters {
        molecule_id: normalized_search_filter(query.molecule_id),
        event_id: normalized_search_filter(query.event_id),
        state: normalized_search_filter(query.state).map(|state| state.to_ascii_uppercase()),
    };
    if let Some(state) = &filters.state {
        if !matches!(
            state.as_str(),
            "REFERENCE" | "ALTERNATE" | "EVENT_ABSENT" | "LOW_QUALITY" | "CONFLICT"
        ) {
            return Err(ApiError::new(
                StatusCode::BAD_REQUEST,
                "evidence search state is not a stored observation state",
            ));
        }
    }
    let search_index = {
        let jobs = state.jobs.read().await;
        let job = jobs
            .get(&job_id)
            .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;
        match job.status {
            JobStatus::Done => job.evidence_search_index.clone().ok_or_else(|| {
                ApiError::new(
                    StatusCode::CONFLICT,
                    "analysis result has no searchable evidence graph",
                )
            })?,
            JobStatus::Error | JobStatus::Cancelled => {
                return Err(ApiError::new(
                    StatusCode::CONFLICT,
                    "analysis has no searchable evidence graph",
                ));
            }
            JobStatus::Queued | JobStatus::Processing => {
                return Err(ApiError::new(
                    StatusCode::ACCEPTED,
                    "analysis result is not ready yet",
                ));
            }
        }
    };
    let mut response =
        Json(search_index.search(filters, query.cursor.unwrap_or(0), limit)).into_response();
    response.headers_mut().insert(
        header::CACHE_CONTROL,
        HeaderValue::from_static("private, immutable, max-age=86400"),
    );
    Ok(response)
}

fn normalized_search_filter(value: Option<String>) -> Option<String> {
    value.and_then(|value| {
        let trimmed = value.trim();
        (!trimmed.is_empty()).then(|| trimmed.to_string())
    })
}

impl EvidenceSearchIndex {
    fn search(
        &self,
        filters: EvidenceSearchFilters,
        cursor: usize,
        limit: usize,
    ) -> EvidenceSearchResponse {
        let molecule_id = filters
            .molecule_id
            .as_ref()
            .and_then(|value| self.string_ids.get(value.as_str()).copied());
        let event_id = filters
            .event_id
            .as_ref()
            .and_then(|value| self.string_ids.get(value.as_str()).copied());
        let state = filters
            .state
            .as_ref()
            .and_then(|value| self.string_ids.get(value.as_str()).copied());
        let unknown_filter = (filters.molecule_id.is_some() && molecule_id.is_none())
            || (filters.event_id.is_some() && event_id.is_none())
            || (filters.state.is_some() && state.is_none());
        if unknown_filter {
            return EvidenceSearchResponse {
                schema_version: "1.0",
                filters,
                cursor,
                next_cursor: None,
                total_matches: 0,
                rows: Vec::new(),
            };
        }

        let mut candidate_postings: Option<&[usize]> = None;
        for postings in [
            molecule_id.and_then(|id| self.by_molecule.get(&id)),
            event_id.and_then(|id| self.by_event.get(&id)),
            state.and_then(|id| self.by_state.get(&id)),
        ]
        .into_iter()
        .flatten()
        {
            if match candidate_postings {
                Some(current) => postings.len() < current.len(),
                None => true,
            } {
                candidate_postings = Some(postings);
            }
        }

        let mut total_matches = 0usize;
        let mut response_rows = Vec::with_capacity(limit);
        let mut visit = |record_index: usize| {
            let record = &self.rows[record_index];
            if molecule_id.is_some_and(|id| record.molecule_id != id)
                || event_id.is_some_and(|id| record.event_id != id)
                || state.is_some_and(|id| record.state != id)
            {
                return;
            }
            if total_matches >= cursor && response_rows.len() < limit {
                response_rows.push(self.response_row(record));
            }
            total_matches += 1;
        };
        if let Some(postings) = candidate_postings {
            for &record_index in postings {
                visit(record_index);
            }
        } else {
            for record_index in 0..self.rows.len() {
                visit(record_index);
            }
        }
        let returned_end = cursor.saturating_add(response_rows.len());
        EvidenceSearchResponse {
            schema_version: "1.0",
            filters,
            cursor,
            next_cursor: (returned_end < total_matches).then_some(returned_end),
            total_matches,
            rows: response_rows,
        }
    }

    fn response_row(&self, record: &EvidenceSearchRecord) -> EvidenceSearchResult {
        let resolve = |id: u32| self.strings[id as usize].to_string();
        EvidenceSearchResult {
            id: format!("observation:{}", record.global_index),
            page_index: record.page_index,
            row_index: record.row_index,
            molecule_id: resolve(record.molecule_id),
            event_id: resolve(record.event_id),
            alignment_id: resolve(record.alignment_id),
            state: resolve(record.state),
            observed_allele: record.observed_allele.map(resolve),
            base_quality: record.base_quality,
            mapping_quality: record.mapping_quality,
            strand: resolve(record.strand),
            evidence_source: resolve(record.evidence_source),
            read_position: record.read_position,
        }
    }
}

async fn evidence_page(
    State(state): State<AppState>,
    AxumPath((job_id, page_index)): AxumPath<(String, String)>,
) -> Result<Response, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let page_index = page_index.parse::<usize>().map_err(|_| {
        ApiError::new(
            StatusCode::BAD_REQUEST,
            "evidence page index must be a non-negative integer",
        )
    })?;
    let page = {
        let jobs = state.jobs.read().await;
        let job = jobs
            .get(&job_id)
            .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;
        match job.status {
            JobStatus::Done => {
                job.evidence_pages.get(page_index).cloned().ok_or_else(|| {
                    ApiError::new(StatusCode::NOT_FOUND, "evidence page not found")
                })?
            }
            JobStatus::Error | JobStatus::Cancelled => {
                return Err(ApiError::new(
                    StatusCode::CONFLICT,
                    "analysis has no evidence pages",
                ));
            }
            JobStatus::Queued | JobStatus::Processing => {
                return Err(ApiError::new(
                    StatusCode::ACCEPTED,
                    "analysis result is not ready yet",
                ));
            }
        }
    };
    let mut response = Response::new(Body::from(page));
    response.headers_mut().insert(
        header::CONTENT_TYPE,
        HeaderValue::from_static("application/json; charset=utf-8"),
    );
    response.headers_mut().insert(
        header::CACHE_CONTROL,
        HeaderValue::from_static("private, immutable, max-age=86400"),
    );
    Ok(response)
}

async fn download(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
) -> Result<Response, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let jobs = state.jobs.read().await;
    let job = jobs
        .get(&job_id)
        .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;
    let html = job
        .html_report
        .clone()
        .ok_or_else(|| ApiError::new(StatusCode::ACCEPTED, "HTML report is not ready yet"))?;

    let mut response = Response::new(Body::from(html));
    response.headers_mut().insert(
        header::CONTENT_TYPE,
        HeaderValue::from_static("text/html; charset=utf-8"),
    );
    response.headers_mut().insert(
        header::CONTENT_DISPOSITION,
        HeaderValue::from_str(&format!("attachment; filename=\"mito-{job_id}.html\""))
            .map_err(internal_error)?,
    );
    Ok(response)
}

async fn run_job(
    state: AppState,
    job_id: Uuid,
    input_path: PathBuf,
    cancel_flag: Arc<AtomicBool>,
    analysis_options: JobAnalysisOptions,
    _job_permit: OwnedSemaphorePermit,
) {
    let Ok(_permit) = state.analysis_slots.clone().acquire_owned().await else {
        update_job(
            &state,
            job_id,
            JobStatus::Error,
            100,
            None,
            Some(FailureBody::new(
                "MITO-API-E2001",
                "analysis scheduler is unavailable",
                true,
            )),
        )
        .await;
        return;
    };
    if is_cancelled(&state, job_id).await {
        state.cancel_flags.write().await.remove(&job_id);
        return;
    }
    update_job(&state, job_id, JobStatus::Processing, 12, None, None).await;

    let input_for_worker = input_path.clone();
    let threads = state.worker_threads;
    let min_mapping_quality = state.min_mapping_quality;
    let min_base_quality = state.min_base_quality;
    let excluded_snp_flags = state.excluded_snp_flags;
    let numt_threshold = state.numt_threshold;
    let max_evidence_observations = state.max_evidence_observations;
    let max_phase_links = state.max_phase_links;
    let analysis = tokio::task::spawn_blocking(move || {
        let engine = MitoEngine::new()?;
        engine.analyze_with_cancel_flag(
            &input_for_worker,
            None,
            AnalyzeOptions {
                filter_numt: true,
                threads,
                min_mapping_quality,
                min_base_quality,
                excluded_snp_flags,
                numt_threshold,
                allow_development_tags: false,
                emit_evidence_graph: analysis_options.emit_evidence_graph,
                max_evidence_observations,
                max_phase_links,
                evidence_page_size: analysis_options.evidence_page_size,
                molecule_id_tag: analysis_options.molecule_id_tag,
                umi_tag: analysis_options.umi_tag,
                duplex_tag: analysis_options.duplex_tag,
            },
            &cancel_flag,
        )
    })
    .await;

    if is_cancelled(&state, job_id).await {
        state.cancel_flags.write().await.remove(&job_id);
        return;
    }

    match analysis {
        Ok(Ok(json_text)) => match serde_json::from_str::<Value>(&json_text) {
            Ok(value) => {
                if let Err(error) = validate_result_contract(&value) {
                    eprintln!("job {job_id} result contract violation: {error}");
                    update_job(
                        &state,
                        job_id,
                        JobStatus::Error,
                        100,
                        None,
                        Some(FailureBody::new(
                            "MITO-E9001",
                            "internal analysis error",
                            false,
                        )),
                    )
                    .await;
                } else {
                    let html = render_html_report(&value);
                    match prepare_completed_artifacts(job_id, &value, json_text, html) {
                        Ok(artifacts) => {
                            update_job(&state, job_id, JobStatus::Done, 100, Some(artifacts), None)
                                .await;
                        }
                        Err(error) => {
                            eprintln!("job {job_id} transport preparation failed: {error}");
                            update_job(
                                &state,
                                job_id,
                                JobStatus::Error,
                                100,
                                None,
                                Some(FailureBody::new(
                                    "MITO-E9001",
                                    "internal analysis error",
                                    false,
                                )),
                            )
                            .await;
                        }
                    }
                }
            }
            Err(error) => {
                eprintln!("job {job_id} returned invalid JSON: {error}");
                update_job(
                    &state,
                    job_id,
                    JobStatus::Error,
                    100,
                    None,
                    Some(FailureBody::new(
                        "MITO-E9001",
                        "internal analysis error",
                        false,
                    )),
                )
                .await;
            }
        },
        Ok(Err(error)) => {
            let failure = public_engine_failure(error, job_id);
            update_job(&state, job_id, JobStatus::Error, 100, None, Some(failure)).await;
        }
        Err(error) => {
            eprintln!("job {job_id} analysis task failed: {error}");
            update_job(
                &state,
                job_id,
                JobStatus::Error,
                100,
                None,
                Some(FailureBody::new(
                    "MITO-E9001",
                    "internal analysis error",
                    false,
                )),
            )
            .await;
        }
    }

    state.cancel_flags.write().await.remove(&job_id);
}

fn prepare_completed_artifacts(
    job_id: Uuid,
    value: &Value,
    json_text: String,
    html: String,
) -> anyhow::Result<CompletedArtifacts> {
    let root = value
        .as_object()
        .ok_or_else(|| anyhow::anyhow!("analysis result root is not an object"))?;
    let evidence_pages = value
        .pointer("/observation_pages")
        .and_then(Value::as_array)
        .map(|pages| {
            pages
                .iter()
                .map(serde_json::to_vec)
                .map(|page| page.map(Bytes::from))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?
        .unwrap_or_default();

    let is_schema_0_6 = value
        .pointer("/metadata/schema_version")
        .and_then(Value::as_str)
        == Some("0.6");
    let evidence_search_index = if is_schema_0_6 {
        Some(Arc::new(build_evidence_search_index(value)?))
    } else {
        None
    };
    let mut summary_root = serde_json::Map::with_capacity(root.len());
    for (name, field) in root {
        if is_schema_0_6 && name == "observation_pages" {
            summary_root.insert(name.clone(), Value::Array(Vec::new()));
        } else {
            summary_root.insert(name.clone(), field.clone());
        }
    }
    let mut summary = Value::Object(summary_root);
    if is_schema_0_6 {
        if !root.contains_key("observation_pages") {
            anyhow::bail!("schema 0.6 summary has no observation_pages");
        }
        let encoding = summary
            .pointer_mut("/evidence_encoding")
            .and_then(Value::as_object_mut)
            .ok_or_else(|| anyhow::anyhow!("schema 0.6 summary has no evidence_encoding"))?;
        encoding.insert(
            "observation_storage".to_string(),
            Value::String("remote_http_pages".to_string()),
        );
        encoding.insert(
            "observation_page_endpoint".to_string(),
            Value::String(format!("/result/{job_id}/evidence/{{page_index}}")),
        );
        encoding.insert(
            "observation_search_endpoint".to_string(),
            Value::String(format!("/result/{job_id}/evidence")),
        );
    }
    Ok(CompletedArtifacts {
        result: Bytes::from(json_text),
        result_summary: Bytes::from(serde_json::to_vec(&summary)?),
        evidence_pages,
        evidence_search_index,
        html_report: Bytes::from(html),
    })
}

fn build_evidence_search_index(value: &Value) -> anyhow::Result<EvidenceSearchIndex> {
    let pages = value
        .pointer("/observation_pages")
        .and_then(Value::as_array)
        .ok_or_else(|| anyhow::anyhow!("schema 0.6 result has no observation_pages"))?;
    let expected_count = value
        .pointer("/evidence_encoding/observation_count")
        .and_then(Value::as_u64)
        .and_then(|count| usize::try_from(count).ok())
        .ok_or_else(|| anyhow::anyhow!("schema 0.6 result has invalid observation_count"))?;
    let mut strings = Vec::<Arc<str>>::new();
    let mut string_ids = HashMap::<Arc<str>, u32>::new();
    let mut rows = Vec::<EvidenceSearchRecord>::with_capacity(expected_count);
    let mut by_molecule = HashMap::<u32, Vec<usize>>::new();
    let mut by_event = HashMap::<u32, Vec<usize>>::new();
    let mut by_state = HashMap::<u32, Vec<usize>>::new();

    for page in pages {
        let page_index = required_page_usize(page, "index")?;
        let offset = required_page_usize(page, "offset")?;
        let count = required_page_usize(page, "count")?;
        let columns = page
            .get("columns")
            .and_then(Value::as_object)
            .ok_or_else(|| anyhow::anyhow!("evidence page {page_index} has no columns"))?;
        let molecule_ids = required_page_column(columns, "molecule_id", page_index, count)?;
        let event_ids = required_page_column(columns, "event_id", page_index, count)?;
        let alignment_ids = required_page_column(columns, "alignment_id", page_index, count)?;
        let states = required_page_column(columns, "state", page_index, count)?;
        let observed_alleles = required_page_column(columns, "observed_allele", page_index, count)?;
        let base_qualities = required_page_column(columns, "base_quality", page_index, count)?;
        let mapping_qualities =
            required_page_column(columns, "mapping_quality", page_index, count)?;
        let strands = required_page_column(columns, "strand", page_index, count)?;
        let evidence_sources = required_page_column(columns, "evidence_source", page_index, count)?;
        let read_positions = required_page_column(columns, "read_position", page_index, count)?;

        for row_index in 0..count {
            let molecule_id = intern_required_column_string(
                molecule_ids,
                row_index,
                "molecule_id",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let event_id = intern_required_column_string(
                event_ids,
                row_index,
                "event_id",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let alignment_id = intern_required_column_string(
                alignment_ids,
                row_index,
                "alignment_id",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let state = intern_required_column_string(
                states,
                row_index,
                "state",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let strand = intern_required_column_string(
                strands,
                row_index,
                "strand",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let evidence_source = intern_required_column_string(
                evidence_sources,
                row_index,
                "evidence_source",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let observed_allele = match &observed_alleles[row_index] {
                Value::Null => None,
                Value::String(value) => Some(intern_evidence_string(
                    value,
                    &mut strings,
                    &mut string_ids,
                )?),
                _ => anyhow::bail!(
                    "evidence page {page_index} observed_allele[{row_index}] is invalid"
                ),
            };
            let base_quality =
                optional_u64_column(base_qualities, row_index, "base_quality", page_index)?;
            let mapping_quality = mapping_qualities[row_index].as_u64().ok_or_else(|| {
                anyhow::anyhow!(
                    "evidence page {page_index} mapping_quality[{row_index}] is invalid"
                )
            })?;
            let read_position =
                optional_f64_column(read_positions, row_index, "read_position", page_index)?;
            let record_index = rows.len();
            rows.push(EvidenceSearchRecord {
                global_index: offset + row_index,
                page_index,
                row_index,
                molecule_id,
                event_id,
                alignment_id,
                state,
                observed_allele,
                base_quality,
                mapping_quality,
                strand,
                evidence_source,
                read_position,
            });
            by_molecule
                .entry(molecule_id)
                .or_default()
                .push(record_index);
            by_event.entry(event_id).or_default().push(record_index);
            by_state.entry(state).or_default().push(record_index);
        }
    }
    if rows.len() != expected_count {
        anyhow::bail!(
            "evidence search index row count {} differs from declared observation_count {expected_count}",
            rows.len()
        );
    }
    Ok(EvidenceSearchIndex {
        strings,
        string_ids,
        rows,
        by_molecule,
        by_event,
        by_state,
    })
}

fn required_page_usize(page: &Value, name: &str) -> anyhow::Result<usize> {
    page.get(name)
        .and_then(Value::as_u64)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or_else(|| anyhow::anyhow!("evidence page has invalid {name}"))
}

fn required_page_column<'a>(
    columns: &'a serde_json::Map<String, Value>,
    name: &str,
    page_index: usize,
    count: usize,
) -> anyhow::Result<&'a Vec<Value>> {
    let values = columns
        .get(name)
        .and_then(Value::as_array)
        .ok_or_else(|| anyhow::anyhow!("evidence page {page_index} has invalid {name} column"))?;
    if values.len() != count {
        anyhow::bail!(
            "evidence page {page_index} {name} column has {} rows, expected {count}",
            values.len()
        );
    }
    Ok(values)
}

fn intern_required_column_string(
    values: &[Value],
    row_index: usize,
    name: &str,
    page_index: usize,
    strings: &mut Vec<Arc<str>>,
    string_ids: &mut HashMap<Arc<str>, u32>,
) -> anyhow::Result<u32> {
    let value = values[row_index].as_str().ok_or_else(|| {
        anyhow::anyhow!("evidence page {page_index} {name}[{row_index}] is invalid")
    })?;
    intern_evidence_string(value, strings, string_ids)
}

fn intern_evidence_string(
    value: &str,
    strings: &mut Vec<Arc<str>>,
    string_ids: &mut HashMap<Arc<str>, u32>,
) -> anyhow::Result<u32> {
    if let Some(id) = string_ids.get(value) {
        return Ok(*id);
    }
    let id = u32::try_from(strings.len())
        .map_err(|_| anyhow::anyhow!("evidence string dictionary exceeds u32"))?;
    let interned: Arc<str> = Arc::from(value);
    strings.push(interned.clone());
    string_ids.insert(interned, id);
    Ok(id)
}

fn optional_u64_column(
    values: &[Value],
    row_index: usize,
    name: &str,
    page_index: usize,
) -> anyhow::Result<Option<u64>> {
    match &values[row_index] {
        Value::Null => Ok(None),
        value => value.as_u64().map(Some).ok_or_else(|| {
            anyhow::anyhow!("evidence page {page_index} {name}[{row_index}] is invalid")
        }),
    }
}

fn optional_f64_column(
    values: &[Value],
    row_index: usize,
    name: &str,
    page_index: usize,
) -> anyhow::Result<Option<f64>> {
    match &values[row_index] {
        Value::Null => Ok(None),
        value => value
            .as_f64()
            .filter(|number| number.is_finite())
            .map(Some)
            .ok_or_else(|| {
                anyhow::anyhow!("evidence page {page_index} {name}[{row_index}] is invalid")
            }),
    }
}

fn public_engine_failure(error: mito_ffi::MitoError, job_id: Uuid) -> FailureBody {
    let message = match error.code.as_str() {
        "MITO-E1101" => "input could not be opened".to_string(),
        "MITO-E1201" => "reference could not be opened".to_string(),
        "MITO-E1301" => "required analysis resource could not be opened".to_string(),
        "MITO-E9001" => {
            eprintln!("job {job_id} internal engine error: {}", error.message);
            "internal analysis error".to_string()
        }
        _ => error.message,
    };
    FailureBody::new(error.code, message, false)
}

async fn is_cancelled(state: &AppState, job_id: Uuid) -> bool {
    state
        .jobs
        .read()
        .await
        .get(&job_id)
        .is_some_and(|job| job.cancel_requested || job.status == JobStatus::Cancelled)
}

async fn update_job(
    state: &AppState,
    job_id: Uuid,
    status: JobStatus,
    progress: u8,
    artifacts: Option<CompletedArtifacts>,
    error: Option<FailureBody>,
) {
    if let Some(job) = state.jobs.write().await.get_mut(&job_id) {
        if job.status == JobStatus::Cancelled && status != JobStatus::Cancelled {
            return;
        }
        job.status = status;
        job.progress = progress.min(100);
        if let Some(artifacts) = artifacts {
            job.result = Some(artifacts.result);
            job.result_summary = Some(artifacts.result_summary);
            job.evidence_pages = artifacts.evidence_pages;
            job.evidence_search_index = artifacts.evidence_search_index;
            job.html_report = Some(artifacts.html_report);
        }
        if error.is_some() {
            job.error = error;
        }
    }
}

fn validate_upload_name(file_name: &str) -> Result<(), ApiError> {
    let path = Path::new(file_name);
    let extension = path
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    match extension.as_str() {
        "fastq" | "fq" | "sam" | "bam" | "cram" => Ok(()),
        _ => Err(ApiError::new(
            StatusCode::BAD_REQUEST,
            "expected .fastq, .fq, .sam, .bam, or .cram input",
        )),
    }
}

fn valid_optional_sam_tag(tag: &str) -> bool {
    if tag.is_empty() {
        return true;
    }
    let bytes = tag.as_bytes();
    bytes.len() == 2 && bytes[0].is_ascii_alphabetic() && bytes[1].is_ascii_alphanumeric()
}

fn parse_job_id(value: &str) -> Result<Uuid, ApiError> {
    value
        .parse::<Uuid>()
        .map_err(|_| ApiError::new(StatusCode::BAD_REQUEST, "job ID must be a UUID"))
}

fn sanitize_file_name(file_name: &str) -> String {
    let sanitized: String = file_name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-') {
                c
            } else {
                '_'
            }
        })
        .collect();
    if sanitized.is_empty() {
        "input.fastq".to_string()
    } else {
        sanitized
    }
}

fn render_html_report(value: &Value) -> String {
    let sample = value
        .pointer("/metadata/sample")
        .and_then(Value::as_str)
        .unwrap_or("mtDNA");
    let pretty = serde_json::to_string_pretty(value).unwrap_or_else(|_| "{}".to_string());
    format!(
        r#"<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Mito-Architect {sample}</title>
<style>
body {{ margin: 0; font-family: ui-sans-serif, system-ui, sans-serif; background: #08111f; color: #e5edf8; }}
header {{ padding: 24px 28px; border-bottom: 1px solid #263244; }}
main {{ padding: 24px 28px; display: grid; gap: 18px; }}
h1 {{ margin: 0; font-size: 32px; letter-spacing: 0; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 12px; }}
.metric {{ border: 1px solid #263244; border-radius: 8px; padding: 16px; background: #101a2b; }}
.metric span {{ color: #94a3b8; font-size: 13px; }}
.metric strong {{ display: block; font-size: 28px; margin-top: 4px; }}
pre {{ overflow: auto; border: 1px solid #263244; border-radius: 8px; padding: 16px; background: #101a2b; }}
</style>
</head>
<body>
<header><h1>Mito-Architect Report</h1></header>
<main>
<section class="grid">
<div class="metric"><span>Sample</span><strong>{sample}</strong></div>
<div class="metric"><span>Clusters</span><strong>{clusters}</strong></div>
<div class="metric"><span>SVs</span><strong>{svs}</strong></div>
<div class="metric"><span>Complex paths</span><strong>{complex_events}</strong></div>
<div class="metric"><span>Passed reads</span><strong>{passed}</strong></div>
</section>
<pre>{pretty}</pre>
</main>
</body>
</html>"#,
        sample = escape_html(sample),
        clusters = value
            .pointer("/clusters")
            .and_then(Value::as_array)
            .map_or(0, Vec::len),
        svs = value
            .pointer("/svs")
            .and_then(Value::as_array)
            .map_or(0, Vec::len),
        complex_events = value
            .pointer("/complex_events")
            .and_then(Value::as_array)
            .map_or(0, Vec::len),
        passed = value
            .pointer("/filter_stats/passed_reads")
            .and_then(Value::as_u64)
            .unwrap_or(0),
        pretty = escape_html(&pretty)
    )
}

fn validate_result_contract(value: &Value) -> Result<(), String> {
    require_object(value, "/metadata")?;
    require_object(value, "/filter_stats")?;
    require_string(value, "/metadata/schema_version")?;
    let schema_version = value
        .pointer("/metadata/schema_version")
        .and_then(Value::as_str)
        .ok_or_else(|| "/metadata/schema_version must be a string".to_string())?;
    if schema_version != "0.5" && schema_version != "0.6" {
        return Err(format!("unsupported result schema {schema_version}"));
    }
    require_string(value, "/metadata/sv_event_schema_version")?;
    require_string(value, "/metadata/complex_sv_event_schema_version")?;
    require_string(value, "/metadata/clinical_annotation_schema_version")?;
    require_string(value, "/metadata/engine_version")?;
    require_u64(value, "/metadata/reference_length")?;
    require_array(value, "/reads")?;
    require_array(value, "/variants")?;
    require_array(value, "/svs")?;
    require_array(value, "/complex_events")?;
    require_array(value, "/clusters")?;
    require_array(value, "/metadata/resources")?;
    require_object(value, "/filter_stats/numt_assessment")?;
    require_u64(value, "/filter_stats/passed_reads")?;
    if schema_version == "0.6" {
        require_object(value, "/evidence_encoding")?;
        require_array(value, "/alignments")?;
        require_array(value, "/molecules")?;
        require_array(value, "/callability")?;
        require_array(value, "/events")?;
        require_array(value, "/observation_pages")?;
        require_array(value, "/phase_links")?;
        require_array(value, "/architectures")?;
        validate_observation_pages(value)?;
        validate_unified_variant_projection(value)?;
    }
    Ok(())
}

fn validate_observation_pages(value: &Value) -> Result<(), String> {
    if value.get("observations").is_some() {
        return Err("schema 0.6 row observations are forbidden; use observation_pages".to_string());
    }
    let encoding = value
        .pointer("/evidence_encoding")
        .and_then(Value::as_object)
        .ok_or_else(|| "/evidence_encoding must be an object".to_string())?;
    if encoding.get("layout").and_then(Value::as_str) != Some("paged_columnar_molecule_event")
        || encoding.get("observation_storage").and_then(Value::as_str)
            != Some("embedded_columnar_pages")
        || encoding
            .get("phase_molecule_policy")
            .and_then(Value::as_str)
            != Some("evidence_eligible_only")
        || encoding
            .get("phase_molecule_reference")
            .and_then(Value::as_str)
            != Some("molecules[].index")
        || encoding.get("phase_null_model").and_then(Value::as_str)
            != Some("independent_marginals_within_jointly_callable")
    {
        return Err("unsupported schema 0.6 evidence page encoding".to_string());
    }
    let expected_count = encoding
        .get("observation_count")
        .and_then(Value::as_u64)
        .ok_or_else(|| "/evidence_encoding/observation_count must be an integer".to_string())?;
    let expected_page_count = encoding
        .get("observation_page_count")
        .and_then(Value::as_u64)
        .ok_or_else(|| {
            "/evidence_encoding/observation_page_count must be an integer".to_string()
        })?;
    let page_size = encoding
        .get("observation_page_size")
        .and_then(Value::as_u64)
        .filter(|value| (1..=1_000_000).contains(value))
        .ok_or_else(|| "/evidence_encoding/observation_page_size is invalid".to_string())?;
    let pages = value
        .pointer("/observation_pages")
        .and_then(Value::as_array)
        .ok_or_else(|| "/observation_pages must be an array".to_string())?;
    if pages.len() as u64 != expected_page_count {
        return Err("observation page cardinality does not match metadata".to_string());
    }
    let column_names = [
        "molecule_id",
        "event_id",
        "alignment_id",
        "state",
        "observed_allele",
        "base_quality",
        "mapping_quality",
        "strand",
        "evidence_source",
        "read_position",
    ];
    let mut offset = 0_u64;
    for (page_index, page) in pages.iter().enumerate() {
        if page.get("index").and_then(Value::as_u64) != Some(page_index as u64)
            || page.get("offset").and_then(Value::as_u64) != Some(offset)
        {
            return Err(format!("observation page {page_index} is not contiguous"));
        }
        let count = page
            .get("count")
            .and_then(Value::as_u64)
            .ok_or_else(|| format!("observation page {page_index} count must be an integer"))?;
        if count > page_size {
            return Err(format!(
                "observation page {page_index} exceeds the page-size limit"
            ));
        }
        let count_usize = usize::try_from(count)
            .map_err(|_| format!("observation page {page_index} count is too large"))?;
        let columns = page
            .get("columns")
            .and_then(Value::as_object)
            .ok_or_else(|| format!("observation page {page_index} columns must be an object"))?;
        for name in column_names {
            if columns.get(name).and_then(Value::as_array).map(Vec::len) != Some(count_usize) {
                return Err(format!(
                    "observation page {page_index} column {name} has the wrong length"
                ));
            }
        }
        offset = offset
            .checked_add(count)
            .ok_or_else(|| "observation page count overflow".to_string())?;
    }
    if offset != expected_count {
        return Err("observation page count does not match observation_count".to_string());
    }
    Ok(())
}

fn validate_unified_variant_projection(value: &Value) -> Result<(), String> {
    #[derive(Clone)]
    struct EventDefinition {
        event_type: String,
        start: u64,
        end: u64,
        reference: String,
        alternate: String,
        support: HashSet<String>,
    }

    let events = value
        .pointer("/events")
        .and_then(Value::as_array)
        .ok_or_else(|| "/events must be an array".to_string())?;
    let mut definitions = HashMap::<String, EventDefinition>::new();
    let mut event_order = HashMap::<String, usize>::new();
    let mut event_complete_callability = HashMap::<String, bool>::new();
    for (index, event) in events.iter().enumerate() {
        let event_id = event
            .get("id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("/events/{index}/id must be a string"))?;
        let event_type = event
            .get("type")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("/events/{index}/type must be a string"))?;
        if event_order.insert(event_id.to_owned(), index).is_some() {
            return Err(format!("duplicate evidence event {event_id}"));
        }
        let assessability = event
            .get("assessability")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("/events/{index}/assessability must be a string"))?;
        event_complete_callability.insert(
            event_id.to_owned(),
            assessability == "REFERENCE_AND_ALTERNATE",
        );
        if !matches!(event_type, "SNV" | "SMALL_INSERTION" | "SMALL_DELETION") {
            continue;
        }
        let support = event
            .get("supporting_molecule_ids")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("/events/{index}/supporting_molecule_ids must be an array"))?
            .iter()
            .map(|item| {
                item.as_str()
                    .map(str::to_owned)
                    .ok_or_else(|| "supporting molecule ID must be a string".to_string())
            })
            .collect::<Result<HashSet<_>, _>>()?;
        let definition = EventDefinition {
            event_type: event_type.to_owned(),
            start: event
                .get("start")
                .and_then(Value::as_u64)
                .ok_or_else(|| format!("/events/{index}/start must be an integer"))?,
            end: event
                .get("end")
                .and_then(Value::as_u64)
                .ok_or_else(|| format!("/events/{index}/end must be an integer"))?,
            reference: event
                .get("ref")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("/events/{index}/ref must be a string"))?
                .to_owned(),
            alternate: event
                .get("alt")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("/events/{index}/alt must be a string"))?
                .to_owned(),
            support,
        };
        if definitions
            .insert(event_id.to_owned(), definition)
            .is_some()
        {
            return Err(format!("duplicate unified variant event {event_id}"));
        }
    }

    let molecules = value
        .pointer("/molecules")
        .and_then(Value::as_array)
        .ok_or_else(|| "/molecules must be an array".to_string())?;
    let mut evidence_eligible = HashSet::new();
    let mut evidence_eligible_molecules = Vec::<(u64, String)>::new();
    for (index, molecule) in molecules.iter().enumerate() {
        let id = molecule
            .get("id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("/molecules/{index}/id must be a string"))?;
        if molecule.get("index").and_then(Value::as_u64) != Some(index as u64) {
            return Err(format!("/molecules/{index}/index is not contiguous"));
        }
        let eligible = molecule
            .get("evidence_eligible")
            .and_then(Value::as_bool)
            .ok_or_else(|| format!("/molecules/{index}/evidence_eligible must be boolean"))?;
        if eligible {
            if !evidence_eligible.insert(id.to_owned()) {
                return Err(format!("duplicate evidence-eligible molecule {id}"));
            }
            evidence_eligible_molecules.push((index as u64, id.to_owned()));
        }
    }
    let mut counts = HashMap::<String, [u64; 6]>::new();
    let mut phase_observations = HashMap::<String, HashMap<String, String>>::new();
    for (page_index, page) in value
        .pointer("/observation_pages")
        .and_then(Value::as_array)
        .ok_or_else(|| "/observation_pages must be an array".to_string())?
        .iter()
        .enumerate()
    {
        let columns = page
            .get("columns")
            .and_then(Value::as_object)
            .ok_or_else(|| format!("observation page {page_index} columns must be an object"))?;
        let event_ids = columns
            .get("event_id")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("observation page {page_index} event_id must be an array"))?;
        let molecule_ids = columns
            .get("molecule_id")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("observation page {page_index} molecule_id must be an array"))?;
        let states = columns
            .get("state")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("observation page {page_index} state must be an array"))?;
        let sources = columns
            .get("evidence_source")
            .and_then(Value::as_array)
            .ok_or_else(|| {
                format!("observation page {page_index} evidence_source must be an array")
            })?;
        for row in 0..event_ids.len() {
            let molecule_id = molecule_ids[row]
                .as_str()
                .ok_or_else(|| "observation molecule ID must be a string".to_string())?;
            let event_id = event_ids[row]
                .as_str()
                .ok_or_else(|| "observation event ID must be a string".to_string())?;
            let Some(definition) = definitions.get(event_id) else {
                continue;
            };
            let state = states[row]
                .as_str()
                .ok_or_else(|| "observation state must be a string".to_string())?;
            let source = sources[row]
                .as_str()
                .ok_or_else(|| "observation evidence source must be a string".to_string())?;
            if !evidence_eligible.contains(molecule_id) {
                return Err(format!(
                    "observation for event {event_id} belongs to an evidence-ineligible molecule"
                ));
            }
            if phase_observations
                .entry(molecule_id.to_owned())
                .or_default()
                .insert(event_id.to_owned(), state.to_owned())
                .is_some()
            {
                return Err(format!(
                    "duplicate molecule/event observation {molecule_id}/{event_id}"
                ));
            }
            let projected = counts.entry(event_id.to_owned()).or_default();
            match state {
                "ALTERNATE" => projected[0] += 1,
                "REFERENCE" => projected[1] += 1,
                "EVENT_ABSENT" => {
                    projected[3] += 1;
                    if definition.event_type == "SNV" || source == "cigar_alternative_small_indel" {
                        projected[2] += 1;
                    } else {
                        projected[1] += 1;
                    }
                }
                "LOW_QUALITY" => projected[4] += 1,
                "CONFLICT" => projected[5] += 1,
                _ => return Err(format!("unsupported observation state {state}")),
            }
        }
    }

    let variants = value
        .pointer("/variants")
        .and_then(Value::as_array)
        .ok_or_else(|| "/variants must be an array".to_string())?;
    if variants.len() != definitions.len() {
        return Err("unified variant/event cardinality mismatch".to_string());
    }
    let mut seen = HashSet::new();
    for (index, variant) in variants.iter().enumerate() {
        let event_id = variant
            .get("event_id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("/variants/{index}/event_id must be a string"))?;
        if !seen.insert(event_id.to_owned()) {
            return Err(format!("duplicate unified variant {event_id}"));
        }
        let definition = definitions
            .get(event_id)
            .ok_or_else(|| format!("variant {event_id} has no source event"))?;
        if variant.get("type").and_then(Value::as_str) != Some(definition.event_type.as_str())
            || variant.get("position").and_then(Value::as_u64) != Some(definition.start)
            || variant.get("start").and_then(Value::as_u64) != Some(definition.start)
            || variant.get("end").and_then(Value::as_u64) != Some(definition.end)
            || variant.get("ref").and_then(Value::as_str) != Some(definition.reference.as_str())
            || variant.get("alt").and_then(Value::as_str) != Some(definition.alternate.as_str())
        {
            return Err(format!("variant/event identity mismatch for {event_id}"));
        }
        let projected = counts.get(event_id).copied().unwrap_or_default();
        let callable = projected[0]
            .checked_add(projected[1])
            .and_then(|total| total.checked_add(projected[2]))
            .ok_or_else(|| "variant count overflow".to_string())?;
        for (field, expected) in [
            ("alt_depth", projected[0]),
            ("ref_depth", projected[1]),
            ("other_depth", projected[2]),
            ("event_absent_depth", projected[3]),
            ("low_quality_depth", projected[4]),
            ("conflict_depth", projected[5]),
            ("callable_depth", callable),
        ] {
            if variant.get(field).and_then(Value::as_u64) != Some(expected) {
                return Err(format!("variant count mismatch for {event_id}/{field}"));
            }
        }
        let support = variant
            .get("supporting_molecule_ids")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("/variants/{index}/supporting_molecule_ids must be an array"))?
            .iter()
            .map(|item| {
                item.as_str()
                    .map(str::to_owned)
                    .ok_or_else(|| "supporting molecule ID must be a string".to_string())
            })
            .collect::<Result<HashSet<_>, _>>()?;
        if support != definition.support {
            return Err(format!("variant support mismatch for {event_id}"));
        }
    }
    if seen.len() != definitions.len() {
        return Err("variant projection does not cover every SNV/small-indel event".to_string());
    }
    let mut expected_phase_pairs = HashSet::<(String, String)>::new();
    for observations in phase_observations.values() {
        for (alternate_event, state) in observations {
            if state != "ALTERNATE" {
                continue;
            }
            for other_event in observations.keys() {
                if alternate_event == other_event {
                    continue;
                }
                let alternate_order = event_order
                    .get(alternate_event)
                    .ok_or_else(|| "phase alternate event is unresolved".to_string())?;
                let other_order = event_order
                    .get(other_event)
                    .ok_or_else(|| "phase neighbor event is unresolved".to_string())?;
                expected_phase_pairs.insert(if alternate_order < other_order {
                    (alternate_event.clone(), other_event.clone())
                } else {
                    (other_event.clone(), alternate_event.clone())
                });
            }
        }
    }
    for (index, link) in value
        .pointer("/phase_links")
        .and_then(Value::as_array)
        .ok_or_else(|| "/phase_links must be an array".to_string())?
        .iter()
        .enumerate()
    {
        let event_a = link
            .get("event_a_id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("phase link {index} event_a_id must be a string"))?;
        let event_b = link
            .get("event_b_id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("phase link {index} event_b_id must be a string"))?;
        let phase_id = link
            .get("id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("phase link {index} ID must be a string"))?;
        if event_order.get(event_a) >= event_order.get(event_b)
            || phase_id != format!("phase:{event_a}|{event_b}")
            || !expected_phase_pairs.remove(&(event_a.to_owned(), event_b.to_owned()))
        {
            return Err(format!(
                "phase link {index} is not a canonical evidence-derived pair"
            ));
        }
        let complete = event_complete_callability.get(event_a) == Some(&true)
            && event_complete_callability.get(event_b) == Some(&true);
        let expected_assessability = if complete {
            "COMPLETE_FOR_BOTH_EVENTS"
        } else {
            "SUPPORT_CONDITIONED"
        };
        if link.get("assessability").and_then(Value::as_str) != Some(expected_assessability) {
            return Err(format!("phase link {index} assessability is inconsistent"));
        }

        let mut projected = [0_u64; 5];
        let mut support = Vec::new();
        let mut uncertain = Vec::new();
        for (molecule_index, molecule_id) in &evidence_eligible_molecules {
            let Some(observations) = phase_observations.get(molecule_id) else {
                continue;
            };
            let (Some(state_a), Some(state_b)) =
                (observations.get(event_a), observations.get(event_b))
            else {
                continue;
            };
            let callable_a = matches!(state_a.as_str(), "REFERENCE" | "ALTERNATE" | "EVENT_ABSENT");
            let callable_b = matches!(state_b.as_str(), "REFERENCE" | "ALTERNATE" | "EVENT_ABSENT");
            if !callable_a || !callable_b {
                projected[4] += 1;
                uncertain.push(*molecule_index);
                continue;
            }
            match (state_a == "ALTERNATE", state_b == "ALTERNATE") {
                (true, true) => {
                    projected[0] += 1;
                    support.push(*molecule_index);
                }
                (true, false) => projected[1] += 1,
                (false, true) => projected[2] += 1,
                (false, false) => projected[3] += 1,
            }
        }
        let jointly_callable = projected[..4].iter().try_fold(0_u64, |total, value| {
            total
                .checked_add(*value)
                .ok_or_else(|| "phase count overflow".to_string())
        })?;
        for (field, expected_count) in [
            ("both_alternate", projected[0]),
            ("a_alternate_b_absent", projected[1]),
            ("a_absent_b_alternate", projected[2]),
            ("neither_alternate", projected[3]),
            ("jointly_uncertain", projected[4]),
            ("jointly_callable", jointly_callable),
        ] {
            if link.get(field).and_then(Value::as_u64) != Some(expected_count) {
                return Err(format!(
                    "phase link {index}/{field} is not evidence-derived"
                ));
            }
        }
        if server_phase_u64_array(link, "supporting_molecule_indices", index)? != support
            || server_phase_u64_array(link, "uncertain_molecule_indices", index)? != uncertain
        {
            return Err(format!(
                "phase link {index} molecule traceability is inconsistent"
            ));
        }
        let expected_qc = [
            (!complete).then_some("SUPPORT_CONDITIONED"),
            (projected[4] != 0).then_some("UNCERTAIN_COOCCURRENCE_EXCLUDED"),
        ]
        .into_iter()
        .flatten()
        .map(str::to_owned)
        .collect::<Vec<_>>();
        if server_phase_string_array(link, "qc_flags", index)? != expected_qc {
            return Err(format!("phase link {index} QC facts are inconsistent"));
        }
        let observed = if jointly_callable == 0 {
            0.0
        } else {
            projected[0] as f64 / jointly_callable as f64
        };
        let expected = if jointly_callable == 0 {
            0.0
        } else {
            let n = jointly_callable as f64;
            ((projected[0] + projected[1]) as f64 / n) * ((projected[0] + projected[2]) as f64 / n)
        };
        let (ci_low, ci_high) = server_phase_wilson_interval(projected[0], jointly_callable);
        for (field, expected_value) in [
            ("co_alternate_fraction", observed),
            ("co_alternate_ci95_low", ci_low),
            ("co_alternate_ci95_high", ci_high),
            ("expected_co_alternate_fraction", expected),
            ("linkage_delta", observed - expected),
        ] {
            let actual = link
                .get(field)
                .and_then(Value::as_f64)
                .ok_or_else(|| format!("phase link {index}/{field} must be a number"))?;
            if (actual - expected_value).abs() > 2e-9 {
                return Err(format!("phase link {index}/{field} is inconsistent"));
            }
        }
    }
    if !expected_phase_pairs.is_empty() {
        return Err("phase projection omits evidence-derived candidate pairs".to_string());
    }
    Ok(())
}

fn server_phase_string_array(
    link: &Value,
    field: &str,
    index: usize,
) -> Result<Vec<String>, String> {
    link.get(field)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("phase link {index}/{field} must be an array"))?
        .iter()
        .map(|value| {
            value
                .as_str()
                .map(str::to_owned)
                .ok_or_else(|| format!("phase link {index}/{field} must contain strings"))
        })
        .collect()
}

fn server_phase_u64_array(link: &Value, field: &str, index: usize) -> Result<Vec<u64>, String> {
    link.get(field)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("phase link {index}/{field} must be an array"))?
        .iter()
        .map(|value| {
            value
                .as_u64()
                .ok_or_else(|| format!("phase link {index}/{field} must contain integers"))
        })
        .collect()
}

fn server_phase_wilson_interval(successes: u64, total: u64) -> (f64, f64) {
    if total == 0 {
        return (0.0, 0.0);
    }
    const Z: f64 = 1.959_963_984_540_054;
    let n = total as f64;
    let p = successes as f64 / n;
    let z2 = Z * Z;
    let denominator = 1.0 + z2 / n;
    let center = (p + z2 / (2.0 * n)) / denominator;
    let margin = Z * ((p * (1.0 - p) / n) + z2 / (4.0 * n * n)).sqrt() / denominator;
    ((center - margin).max(0.0), (center + margin).min(1.0))
}

fn require_object(value: &Value, pointer: &str) -> Result<(), String> {
    value
        .pointer(pointer)
        .and_then(Value::as_object)
        .map(|_| ())
        .ok_or_else(|| format!("{pointer} must be an object"))
}

fn require_array(value: &Value, pointer: &str) -> Result<(), String> {
    value
        .pointer(pointer)
        .and_then(Value::as_array)
        .map(|_| ())
        .ok_or_else(|| format!("{pointer} must be an array"))
}

fn require_string(value: &Value, pointer: &str) -> Result<(), String> {
    value
        .pointer(pointer)
        .and_then(Value::as_str)
        .map(|_| ())
        .ok_or_else(|| format!("{pointer} must be a string"))
}

fn require_u64(value: &Value, pointer: &str) -> Result<(), String> {
    value
        .pointer(pointer)
        .and_then(Value::as_u64)
        .map(|_| ())
        .ok_or_else(|| format!("{pointer} must be an unsigned integer"))
}

fn escape_html(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}

fn spawn_cleanup(state: AppState) {
    tokio::spawn(async move {
        loop {
            tokio::time::sleep(Duration::from_secs(60 * 60)).await;
            let expired = {
                let jobs = state.jobs.read().await;
                jobs.iter()
                    .filter_map(|(job_id, job)| {
                        (now_epoch_secs().saturating_sub(job.created_at) > JOB_TTL_SECS)
                            .then_some((*job_id, job.input_path.clone()))
                    })
                    .collect::<Vec<_>>()
            };

            if expired.is_empty() {
                continue;
            }

            {
                let cancel_flags = state.cancel_flags.read().await;
                for (job_id, _) in &expired {
                    if let Some(flag) = cancel_flags.get(job_id) {
                        flag.store(true, Ordering::Relaxed);
                    }
                }
            }
            {
                let mut jobs = state.jobs.write().await;
                for (job_id, _) in &expired {
                    jobs.remove(job_id);
                }
            }
            {
                let mut cancel_flags = state.cancel_flags.write().await;
                for (job_id, _) in &expired {
                    cancel_flags.remove(job_id);
                }
            }
            for (_, input_path) in expired {
                if let Some(parent) = input_path.parent() {
                    let _ = fs::remove_dir_all(parent).await;
                }
            }
        }
    });
}

fn now_epoch_secs() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |duration| duration.as_secs())
}

fn positive_env_usize(name: &str) -> anyhow::Result<Option<usize>> {
    let value = match std::env::var(name) {
        Ok(value) => value,
        Err(std::env::VarError::NotPresent) => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    let parsed = value
        .parse::<usize>()
        .map_err(|_| anyhow::anyhow!("{name} must be a positive integer"))?;
    if parsed == 0 {
        anyhow::bail!("{name} must be a positive integer");
    }
    Ok(Some(parsed))
}

fn env_u8(name: &str) -> anyhow::Result<Option<u8>> {
    optional_env_integer(name, |value| value.parse::<u8>())
}

fn env_u16(name: &str) -> anyhow::Result<Option<u16>> {
    let value = match std::env::var(name) {
        Ok(value) => value,
        Err(std::env::VarError::NotPresent) => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    let parsed = if let Some(hex) = value
        .strip_prefix("0x")
        .or_else(|| value.strip_prefix("0X"))
    {
        u16::from_str_radix(hex, 16)
    } else {
        value.parse::<u16>()
    }
    .map_err(|_| anyhow::anyhow!("{name} must be an unsigned 16-bit integer"))?;
    Ok(Some(parsed))
}

fn optional_env_number(name: &str) -> anyhow::Result<Option<f64>> {
    let value = match std::env::var(name) {
        Ok(value) => value,
        Err(std::env::VarError::NotPresent) => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    let parsed = value
        .parse::<f64>()
        .map_err(|_| anyhow::anyhow!("{name} must be a finite number"))?;
    if !parsed.is_finite() {
        anyhow::bail!("{name} must be a finite number");
    }
    Ok(Some(parsed))
}

fn optional_env_integer<T, E>(
    name: &str,
    parse: impl FnOnce(&str) -> Result<T, E>,
) -> anyhow::Result<Option<T>> {
    let value = match std::env::var(name) {
        Ok(value) => value,
        Err(std::env::VarError::NotPresent) => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    parse(&value)
        .map(Some)
        .map_err(|_| anyhow::anyhow!("{name} has an invalid integer value"))
}

fn cors_layer() -> anyhow::Result<CorsLayer> {
    let configured = std::env::var("MITO_CORS_ORIGINS")
        .unwrap_or_else(|_| "http://127.0.0.1:5173,http://localhost:5173".to_string());
    let origins = configured
        .split(',')
        .map(str::trim)
        .filter(|origin| !origin.is_empty())
        .map(str::parse::<HeaderValue>)
        .collect::<Result<Vec<_>, _>>()?;
    if origins.is_empty() {
        anyhow::bail!("MITO_CORS_ORIGINS must contain at least one valid origin");
    }
    Ok(CorsLayer::new()
        .allow_origin(origins)
        .allow_methods([Method::GET, Method::POST])
        .allow_headers([header::CONTENT_TYPE, header::AUTHORIZATION]))
}

fn internal_error(error: impl std::fmt::Display) -> ApiError {
    eprintln!("internal server error: {error}");
    ApiError::coded(
        StatusCode::INTERNAL_SERVER_ERROR,
        "MITO-API-E9001",
        "internal server error",
        false,
    )
}

#[cfg(test)]
mod tests {
    use super::{
        build_app, constant_time_eq, public_engine_failure, sanitize_file_name,
        validate_result_contract, validate_upload_name, ApiError, AppState, JobRecord, JobStatus,
    };
    use axum::body::{to_bytes, Body, Bytes};
    use axum::http::{header, Request, StatusCode};
    use serde_json::{json, Value};
    use std::collections::HashMap;
    use std::sync::Arc;
    use tokio::sync::{RwLock, Semaphore};
    use tower::ServiceExt;
    use tower_http::cors::CorsLayer;
    use uuid::Uuid;

    fn test_state(api_key: Option<&str>) -> AppState {
        let tmp_dir =
            std::env::temp_dir().join(format!("mito-server-http-test-{}", Uuid::new_v4()));
        std::fs::create_dir_all(&tmp_dir).expect("test temporary directory should be created");
        AppState {
            jobs: Arc::new(RwLock::new(HashMap::new())),
            cancel_flags: Arc::new(RwLock::new(HashMap::new())),
            analysis_slots: Arc::new(Semaphore::new(2)),
            job_slots: Arc::new(Semaphore::new(4)),
            worker_threads: 1,
            min_mapping_quality: 20,
            min_base_quality: 10,
            excluded_snp_flags: 0xF00,
            numt_threshold: 0.30,
            max_evidence_observations: 5_000_000,
            max_phase_links: 1_000_000,
            evidence_page_size: 4096,
            engine_version: Arc::from("0.5.0-dev"),
            schema_version: Arc::from("0.5"),
            error_schema_version: Arc::from("1.0"),
            htslib_enabled: true,
            api_key: api_key.map(Arc::<str>::from),
            max_upload_bytes: 1024 * 1024,
            tmp_dir,
        }
    }

    async fn response_json(response: axum::response::Response) -> serde_json::Value {
        let bytes = to_bytes(response.into_body(), 1024 * 1024)
            .await
            .expect("response body should be readable");
        serde_json::from_slice(&bytes).expect("response body should be JSON")
    }

    fn get_request(uri: &str, token: Option<&str>) -> Request<Body> {
        let mut builder = Request::builder().method("GET").uri(uri);
        if let Some(token) = token {
            builder = builder.header(header::AUTHORIZATION, format!("Bearer {token}"));
        }
        builder
            .body(Body::empty())
            .expect("test request should be valid")
    }

    fn post_request(uri: &str, token: Option<&str>) -> Request<Body> {
        let mut builder = Request::builder().method("POST").uri(uri);
        if let Some(token) = token {
            builder = builder.header(header::AUTHORIZATION, format!("Bearer {token}"));
        }
        builder
            .body(Body::empty())
            .expect("test request should be valid")
    }

    fn multipart_request(file_name: &str, content: &str) -> Request<Body> {
        let boundary = "mito-architect-test-boundary";
        let body = format!(
            "--{boundary}\r\nContent-Disposition: form-data; name=\"file\"; filename=\"{file_name}\"\r\nContent-Type: application/octet-stream\r\n\r\n{content}\r\n--{boundary}--\r\n"
        );
        Request::builder()
            .method("POST")
            .uri("/upload")
            .header(header::AUTHORIZATION, "Bearer laboratory-secret")
            .header(
                header::CONTENT_TYPE,
                format!("multipart/form-data; boundary={boundary}"),
            )
            .body(Body::from(body))
            .expect("multipart test request should be valid")
    }

    fn multipart_evidence_request(file_first: bool) -> Request<Body> {
        let boundary = "mito-architect-evidence-boundary";
        let file = format!(
            "--{boundary}\r\nContent-Disposition: form-data; name=\"file\"; filename=\"phase.sam\"\r\nContent-Type: application/octet-stream\r\n\r\n@SQ\tSN:MT\tLN:16569\nphase\t0\tMT\t1\t60\t12M\t*\t0\t0\tGAACGCAGGTCT\tIIIIIIIIIIII\tMI:Z:M1\tRX:Z:AAA\tDX:Z:duplex\n\r\n"
        );
        let options = format!(
            "--{boundary}\r\nContent-Disposition: form-data; name=\"evidence_graph\"\r\n\r\ntrue\r\n--{boundary}\r\nContent-Disposition: form-data; name=\"evidence_page_size\"\r\n\r\n1\r\n--{boundary}\r\nContent-Disposition: form-data; name=\"molecule_id_tag\"\r\n\r\nMI\r\n--{boundary}\r\nContent-Disposition: form-data; name=\"umi_tag\"\r\n\r\nRX\r\n--{boundary}\r\nContent-Disposition: form-data; name=\"duplex_tag\"\r\n\r\nDX\r\n"
        );
        let body = if file_first {
            format!("{file}{options}--{boundary}--\r\n")
        } else {
            format!("{options}{file}--{boundary}--\r\n")
        };
        Request::builder()
            .method("POST")
            .uri("/upload")
            .header(header::AUTHORIZATION, "Bearer laboratory-secret")
            .header(
                header::CONTENT_TYPE,
                format!("multipart/form-data; boundary={boundary}"),
            )
            .body(Body::from(body))
            .expect("evidence multipart request should be valid")
    }

    #[test]
    fn result_contract_requires_runtime_metadata() {
        let valid = json!({
            "metadata": {
                "schema_version": "0.5",
                "sv_event_schema_version": "1.0",
                "complex_sv_event_schema_version": "1.0",
                "clinical_annotation_schema_version": "1.0",
                "engine_version": "0.5.0-dev",
                "reference_length": 16569,
                "resources": []
            },
            "filter_stats": {
                "passed_reads": 2,
                "numt_assessment": { "mode": "competitive_alignment" }
            },
            "reads": [],
            "variants": [],
            "svs": [],
            "complex_events": [],
            "clusters": []
        });
        validate_result_contract(&valid).expect("valid result should satisfy contract");

        let missing_schema = json!({
            "metadata": {
                "engine_version": "0.5.0-dev",
                "reference_length": 16569
            },
            "filter_stats": {
                "passed_reads": 2,
                "numt_assessment": { "mode": "competitive_alignment" }
            },
            "reads": [],
            "variants": [],
            "svs": [],
            "clusters": []
        });
        let error =
            validate_result_contract(&missing_schema).expect_err("schema_version must be required");
        assert!(error.contains("/metadata/schema_version"));
    }

    #[test]
    fn upload_names_are_sanitized_and_extension_checked() {
        assert_eq!(
            sanitize_file_name("../patient sample.sam"),
            ".._patient_sample.sam"
        );
        validate_upload_name("sample.fastq").expect("FASTQ should be accepted");
        validate_upload_name("sample.bam").expect("BAM should be accepted");

        let error = validate_upload_name("sample.exe").expect_err("unknown inputs rejected");
        assert_eq!(error.status, StatusCode::BAD_REQUEST);
    }

    #[test]
    fn api_key_comparison_rejects_length_and_content_mismatches() {
        assert!(constant_time_eq(b"laboratory-secret", b"laboratory-secret"));
        assert!(!constant_time_eq(
            b"laboratory-secret",
            b"laboratory-secreu"
        ));
        assert!(!constant_time_eq(
            b"laboratory-secret",
            b"laboratory-secret-extra"
        ));
    }

    #[test]
    fn api_errors_have_a_stable_machine_readable_envelope() {
        let error = ApiError::new(StatusCode::PAYLOAD_TOO_LARGE, "too large");
        assert_eq!(error.body.schema_version, "1.0");
        assert_eq!(error.body.code, "MITO-API-E1005");
        assert_eq!(error.body.message, "too large");
        assert!(!error.body.retryable);
    }

    #[test]
    fn engine_open_failures_do_not_expose_server_paths() {
        let failure = public_engine_failure(
            mito_ffi::MitoError::new(
                "MITO-E1101",
                "could not open input file: /tmp/private/sample.bam",
            ),
            uuid::Uuid::nil(),
        );
        assert_eq!(failure.code, "MITO-E1101");
        assert_eq!(failure.message, "input could not be opened");
        assert!(!failure.message.contains("/tmp/private"));
    }

    #[tokio::test]
    async fn router_preserves_http_contract_across_auth_and_job_states() {
        let state = test_state(Some("laboratory-secret"));
        let tmp_dir = state.tmp_dir.clone();
        let app = build_app(state.clone(), CorsLayer::new());

        let response = app
            .clone()
            .oneshot(get_request("/healthz", None))
            .await
            .expect("health request should complete");
        assert_eq!(response.status(), StatusCode::OK);
        assert_eq!(response_json(response).await, json!({ "status": "ok" }));

        let response = app
            .clone()
            .oneshot(get_request("/readyz", None))
            .await
            .expect("readiness request should complete");
        assert_eq!(response.status(), StatusCode::OK);
        assert_eq!(
            response_json(response).await,
            json!({
                "status": "ready",
                "engine_version": "0.5.0-dev",
                "schema_version": "0.5",
                "error_schema_version": "1.0",
                "htslib_enabled": true,
                "available_analysis_slots": 2,
                "available_job_slots": 4
            })
        );

        let unknown_id = Uuid::nil();
        let response = app
            .clone()
            .oneshot(get_request(&format!("/status/{unknown_id}"), None))
            .await
            .expect("unauthorized request should complete");
        assert_eq!(response.status(), StatusCode::UNAUTHORIZED);
        assert_eq!(
            response_json(response).await,
            json!({
                "error": {
                    "schema_version": "1.0",
                    "code": "MITO-API-E1002",
                    "message": "invalid or missing API key",
                    "retryable": false
                }
            })
        );

        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/status/{unknown_id}"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("not-found request should complete");
        assert_eq!(response.status(), StatusCode::NOT_FOUND);
        assert_eq!(
            response_json(response).await,
            json!({
                "error": {
                    "schema_version": "1.0",
                    "code": "MITO-API-E1003",
                    "message": "job not found",
                    "retryable": false
                }
            })
        );

        let response = app
            .clone()
            .oneshot(get_request("/status/not-a-uuid", Some("laboratory-secret")))
            .await
            .expect("malformed job ID request should complete");
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);
        assert_eq!(
            response_json(response).await,
            json!({
                "error": {
                    "schema_version": "1.0",
                    "code": "MITO-API-E1001",
                    "message": "job ID must be a UUID",
                    "retryable": false
                }
            })
        );

        let response = app
            .clone()
            .oneshot(post_request("/upload", Some("laboratory-secret")))
            .await
            .expect("malformed multipart request should complete");
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);
        assert_eq!(
            response_json(response).await,
            json!({
                "error": {
                    "schema_version": "1.0",
                    "code": "MITO-API-E1001",
                    "message": "invalid multipart upload",
                    "retryable": false
                }
            })
        );

        let response = app
            .clone()
            .oneshot(multipart_request("payload.exe", "x"))
            .await
            .expect("invalid extension upload should complete");
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);
        assert_eq!(
            response_json(response).await,
            json!({
                "error": {
                    "schema_version": "1.0",
                    "code": "MITO-API-E1001",
                    "message": "expected .fastq, .fq, .sam, .bam, or .cram input",
                    "retryable": false
                }
            })
        );

        let job_id = Uuid::new_v4();
        state.jobs.write().await.insert(
            job_id,
            JobRecord {
                status: JobStatus::Queued,
                progress: 0,
                input_path: tmp_dir.join("input.sam"),
                result: None,
                result_summary: None,
                evidence_pages: Vec::new(),
                evidence_search_index: None,
                html_report: None,
                error: None,
                cancel_requested: false,
                created_at: 1,
            },
        );
        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/status/{job_id}"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("status request should complete");
        assert_eq!(response.status(), StatusCode::OK);
        assert_eq!(
            response_json(response).await,
            json!({
                "job_id": job_id,
                "status": "queued",
                "progress": 0,
                "error": null
            })
        );

        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/result/{job_id}"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("pending result request should complete");
        assert_eq!(response.status(), StatusCode::ACCEPTED);
        assert_eq!(
            response_json(response).await,
            json!({
                "error": {
                    "schema_version": "1.0",
                    "code": "MITO-API-E1007",
                    "message": "analysis result is not ready yet",
                    "retryable": true
                }
            })
        );

        {
            let mut jobs = state.jobs.write().await;
            let job = jobs.get_mut(&job_id).expect("job should exist");
            job.status = JobStatus::Done;
            job.progress = 100;
            job.result = Some(Bytes::from_static(
                br#"{"metadata":{"schema_version":"0.5"}}"#,
            ));
            job.result_summary = job.result.clone();
            job.html_report = Some(Bytes::from_static(b"<!doctype html><title>result</title>"));
        }
        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/result/{job_id}"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("completed result request should complete");
        assert_eq!(response.status(), StatusCode::OK);
        assert_eq!(
            response
                .headers()
                .get(header::CONTENT_TYPE)
                .and_then(|value| value.to_str().ok()),
            Some("application/json; charset=utf-8")
        );
        assert_eq!(
            response_json(response).await,
            json!({ "metadata": { "schema_version": "0.5" } })
        );

        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/download/{job_id}"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("download request should complete");
        assert_eq!(response.status(), StatusCode::OK);
        assert_eq!(
            response
                .headers()
                .get(header::CONTENT_TYPE)
                .and_then(|value| value.to_str().ok()),
            Some("text/html; charset=utf-8")
        );
        let content_disposition = response
            .headers()
            .get(header::CONTENT_DISPOSITION)
            .and_then(|value| value.to_str().ok())
            .expect("download should include content disposition");
        assert_eq!(
            content_disposition,
            format!("attachment; filename=\"mito-{job_id}.html\"")
        );
        let body = to_bytes(response.into_body(), 1024)
            .await
            .expect("download body should be readable");
        assert_eq!(
            body,
            Bytes::from_static(b"<!doctype html><title>result</title>")
        );

        let response = app
            .oneshot(post_request(
                &format!("/cancel/{job_id}"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("completed cancellation request should complete");
        assert_eq!(response.status(), StatusCode::CONFLICT);
        assert_eq!(
            response_json(response).await,
            json!({
                "error": {
                    "schema_version": "1.0",
                    "code": "MITO-API-E1004",
                    "message": "completed jobs cannot be cancelled",
                    "retryable": false
                }
            })
        );

        std::fs::remove_dir_all(tmp_dir).expect("test temporary directory should be removed");
    }

    #[tokio::test]
    async fn evidence_graph_selection_is_independent_of_multipart_field_order() {
        for file_first in [false, true] {
            let state = test_state(Some("laboratory-secret"));
            let tmp_dir = state.tmp_dir.clone();
            let app = build_app(state.clone(), CorsLayer::new());
            let response = app
                .clone()
                .oneshot(multipart_evidence_request(file_first))
                .await
                .expect("evidence upload should complete");
            assert_eq!(response.status(), StatusCode::OK);
            let uploaded = response_json(response).await;
            assert_eq!(
                uploaded
                    .get("result_schema")
                    .and_then(serde_json::Value::as_str),
                Some("0.6")
            );
            let job_id = uploaded
                .get("job_id")
                .and_then(serde_json::Value::as_str)
                .expect("upload should return a job ID");

            let mut completed = false;
            for _ in 0..2_000 {
                let response = app
                    .clone()
                    .oneshot(get_request(
                        &format!("/status/{job_id}"),
                        Some("laboratory-secret"),
                    ))
                    .await
                    .expect("status request should complete");
                let status = response_json(response).await;
                match status.get("status").and_then(serde_json::Value::as_str) {
                    Some("done") => {
                        completed = true;
                        break;
                    }
                    Some("error") | Some("cancelled") => {
                        panic!("schema 0.6 upload failed: {status}");
                    }
                    _ => tokio::time::sleep(std::time::Duration::from_millis(10)).await,
                }
            }
            assert!(completed, "schema 0.6 upload did not complete in time");
            let response = app
                .clone()
                .oneshot(get_request(
                    &format!("/result/{job_id}"),
                    Some("laboratory-secret"),
                ))
                .await
                .expect("result request should complete");
            assert_eq!(response.status(), StatusCode::OK);
            let result = response_json(response).await;
            assert_eq!(
                result
                    .pointer("/metadata/schema_version")
                    .and_then(serde_json::Value::as_str),
                Some("0.6")
            );
            assert_eq!(
                result
                    .pointer("/molecules/0/identity_policy")
                    .and_then(serde_json::Value::as_str),
                Some("sam_tag:MI")
            );
            let response = app
                .clone()
                .oneshot(get_request(
                    &format!("/result/{job_id}/summary"),
                    Some("laboratory-secret"),
                ))
                .await
                .expect("result summary request should complete");
            assert_eq!(response.status(), StatusCode::OK);
            let summary = response_json(response).await;
            assert_eq!(
                summary
                    .pointer("/evidence_encoding/observation_storage")
                    .and_then(Value::as_str),
                Some("remote_http_pages")
            );
            let expected_page_endpoint = format!("/result/{job_id}/evidence/{{page_index}}");
            assert_eq!(
                summary
                    .pointer("/evidence_encoding/observation_page_endpoint")
                    .and_then(Value::as_str),
                Some(expected_page_endpoint.as_str())
            );
            let expected_search_endpoint = format!("/result/{job_id}/evidence");
            assert_eq!(
                summary
                    .pointer("/evidence_encoding/observation_search_endpoint")
                    .and_then(Value::as_str),
                Some(expected_search_endpoint.as_str())
            );
            assert_eq!(
                summary
                    .pointer("/observation_pages")
                    .and_then(Value::as_array)
                    .map(Vec::len),
                Some(0)
            );
            let expected_search_count = summary
                .pointer("/evidence_encoding/observation_count")
                .and_then(Value::as_u64)
                .expect("summary should declare the global observation count");
            let response = app
                .clone()
                .oneshot(get_request(
                    &format!("/result/{job_id}/evidence/0"),
                    Some("laboratory-secret"),
                ))
                .await
                .expect("evidence page request should complete");
            assert_eq!(response.status(), StatusCode::OK);
            let page = response_json(response).await;
            assert_eq!(page.get("index").and_then(Value::as_u64), Some(0));
            assert_eq!(page.get("count").and_then(Value::as_u64), Some(1));
            let expected_state = page
                .pointer("/columns/state/0")
                .and_then(Value::as_str)
                .expect("evidence page should contain one stored state")
                .to_string();
            let response = app
                .clone()
                .oneshot(get_request(
                    &format!("/result/{job_id}/evidence?molecule_id=MI:M1&limit=1"),
                    Some("laboratory-secret"),
                ))
                .await
                .expect("global evidence search should complete");
            assert_eq!(response.status(), StatusCode::OK);
            let search = response_json(response).await;
            assert_eq!(
                search.get("schema_version").and_then(Value::as_str),
                Some("1.0")
            );
            assert_eq!(
                search.get("total_matches").and_then(Value::as_u64),
                Some(expected_search_count)
            );
            assert_eq!(
                search
                    .pointer("/rows/0/molecule_id")
                    .and_then(Value::as_str),
                Some("MI:M1")
            );
            assert_eq!(
                search.pointer("/rows/0/page_index").and_then(Value::as_u64),
                Some(0)
            );
            assert_eq!(
                search.pointer("/rows/0/row_index").and_then(Value::as_u64),
                Some(0)
            );
            assert_eq!(
                search.pointer("/rows/0/state").and_then(Value::as_str),
                Some(expected_state.as_str())
            );
            assert_eq!(search.get("next_cursor").and_then(Value::as_u64), Some(1));
            let response = app
                .clone()
                .oneshot(get_request(
                    &format!("/result/{job_id}/evidence?molecule_id=MI:M1&cursor=1&limit=1"),
                    Some("laboratory-secret"),
                ))
                .await
                .expect("second evidence-search cursor should complete");
            assert_eq!(response.status(), StatusCode::OK);
            let second_search_page = response_json(response).await;
            assert_eq!(
                second_search_page
                    .pointer("/rows/0/page_index")
                    .and_then(Value::as_u64),
                Some(1)
            );
            assert!(second_search_page
                .get("next_cursor")
                .is_some_and(Value::is_null));
            let response = app
                .clone()
                .oneshot(get_request(
                    &format!("/result/{job_id}/evidence?event_id=missing"),
                    Some("laboratory-secret"),
                ))
                .await
                .expect("empty evidence search should complete");
            assert_eq!(response.status(), StatusCode::OK);
            assert_eq!(
                response_json(response)
                    .await
                    .get("total_matches")
                    .and_then(Value::as_u64),
                Some(0)
            );
            let response = app
                .clone()
                .oneshot(get_request(
                    &format!("/result/{job_id}/evidence?state=NOT_CALLABLE"),
                    Some("laboratory-secret"),
                ))
                .await
                .expect("invalid sparse-state search should complete");
            assert_eq!(response.status(), StatusCode::BAD_REQUEST);
            let response = app
                .clone()
                .oneshot(get_request(
                    &format!("/result/{job_id}/evidence/999"),
                    Some("laboratory-secret"),
                ))
                .await
                .expect("missing evidence page request should complete");
            assert_eq!(response.status(), StatusCode::NOT_FOUND);
            std::fs::remove_dir_all(tmp_dir).expect("test temporary directory should be removed");
        }
    }
}
