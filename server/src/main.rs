use axum::body::{Body, Bytes};
use axum::extract::{DefaultBodyLimit, Multipart, Path as AxumPath, Request, State};
use axum::http::{header, HeaderValue, Method, StatusCode};
use axum::middleware::{self, Next};
use axum::response::{IntoResponse, Response};
use axum::routing::{get, post};
use axum::{Json, Router};
use mito_ffi::{AnalyzeOptions, MitoEngine};
use serde::Serialize;
use serde_json::Value;
use std::collections::HashMap;
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
    html_report: Option<Bytes>,
    error: Option<FailureBody>,
    cancel_requested: bool,
    created_at: u64,
}

#[derive(Debug, Serialize)]
struct UploadResponse {
    job_id: Uuid,
    status: JobStatus,
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

    let protected = Router::new()
        .route("/upload", post(upload))
        .route("/status/:job_id", get(status))
        .route("/result/:job_id", get(result))
        .route("/download/:job_id", get(download))
        .route("/cancel/:job_id", post(cancel))
        .route_layer(middleware::from_fn_with_state(
            state.clone(),
            require_api_key,
        ));
    let app = Router::new()
        .route("/healthz", get(health))
        .route("/readyz", get(ready))
        .merge(protected)
        .layer(DefaultBodyLimit::max(max_upload_bytes))
        .layer(cors)
        .with_state(state);

    let listener = tokio::net::TcpListener::bind(addr).await?;
    println!("mito-server listening on http://{addr}");
    axum::serve(listener, app).await?;
    Ok(())
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
    mut multipart: Multipart,
) -> Result<Json<UploadResponse>, ApiError> {
    let job_permit = state
        .job_slots
        .clone()
        .try_acquire_owned()
        .map_err(|_| ApiError::new(StatusCode::TOO_MANY_REQUESTS, "analysis queue is full"))?;
    while let Some(mut field) = multipart
        .next_field()
        .await
        .map_err(|error| ApiError::new(StatusCode::BAD_REQUEST, error.to_string()))?
    {
        if field.name() != Some("file") {
            continue;
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

        let record = JobRecord {
            status: JobStatus::Queued,
            progress: 0,
            input_path: input_path.clone(),
            result: None,
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
        tokio::spawn(async move {
            run_job(job_state, job_id, input_path, cancel_flag, job_permit).await;
        });

        return Ok(Json(UploadResponse {
            job_id,
            status: JobStatus::Queued,
        }));
    }

    Err(ApiError::new(
        StatusCode::BAD_REQUEST,
        "multipart upload must include a file field named 'file'",
    ))
}

async fn status(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<Uuid>,
) -> Result<Json<StatusResponse>, ApiError> {
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
    AxumPath(job_id): AxumPath<Uuid>,
) -> Result<Json<StatusResponse>, ApiError> {
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
    AxumPath(job_id): AxumPath<Uuid>,
) -> Result<Response, ApiError> {
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

async fn download(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<Uuid>,
) -> Result<Response, ApiError> {
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
    _job_permit: OwnedSemaphorePermit,
) {
    let Ok(_permit) = state.analysis_slots.clone().acquire_owned().await else {
        update_job(
            &state,
            job_id,
            JobStatus::Error,
            100,
            None,
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
    update_job(&state, job_id, JobStatus::Processing, 12, None, None, None).await;

    let input_for_worker = input_path.clone();
    let threads = state.worker_threads;
    let min_mapping_quality = state.min_mapping_quality;
    let min_base_quality = state.min_base_quality;
    let excluded_snp_flags = state.excluded_snp_flags;
    let numt_threshold = state.numt_threshold;
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
                    update_job(
                        &state,
                        job_id,
                        JobStatus::Done,
                        100,
                        Some(Bytes::from(json_text)),
                        Some(Bytes::from(html)),
                        None,
                    )
                    .await;
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
            update_job(
                &state,
                job_id,
                JobStatus::Error,
                100,
                None,
                None,
                Some(failure),
            )
            .await;
        }
        Err(error) => {
            eprintln!("job {job_id} analysis task failed: {error}");
            update_job(
                &state,
                job_id,
                JobStatus::Error,
                100,
                None,
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
    result: Option<Bytes>,
    html_report: Option<Bytes>,
    error: Option<FailureBody>,
) {
    if let Some(job) = state.jobs.write().await.get_mut(&job_id) {
        if job.status == JobStatus::Cancelled && status != JobStatus::Cancelled {
            return;
        }
        job.status = status;
        job.progress = progress.min(100);
        if result.is_some() {
            job.result = result;
        }
        if html_report.is_some() {
            job.html_report = html_report;
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
    require_string(value, "/metadata/engine_version")?;
    require_u64(value, "/metadata/reference_length")?;
    require_array(value, "/reads")?;
    require_array(value, "/variants")?;
    require_array(value, "/svs")?;
    require_array(value, "/clusters")?;
    require_array(value, "/metadata/resources")?;
    require_object(value, "/filter_stats/numt_assessment")?;
    require_u64(value, "/filter_stats/passed_reads")?;
    Ok(())
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
        constant_time_eq, public_engine_failure, sanitize_file_name, validate_result_contract,
        validate_upload_name, ApiError,
    };
    use axum::http::StatusCode;
    use serde_json::json;

    #[test]
    fn result_contract_requires_runtime_metadata() {
        let valid = json!({
            "metadata": {
                "schema_version": "0.4",
                "engine_version": "0.4.0",
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
            "clusters": []
        });
        validate_result_contract(&valid).expect("valid result should satisfy contract");

        let missing_schema = json!({
            "metadata": {
                "engine_version": "0.4.0",
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
}
