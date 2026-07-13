use std::ffi::{CStr, CString};
use std::fmt;
use std::os::raw::{c_char, c_void};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

unsafe extern "C" {
    fn mito_engine_new() -> *mut c_void;
    fn mito_engine_delete(engine: *mut c_void);
    fn mito_engine_has_htslib() -> bool;
    fn mito_engine_version() -> *const c_char;
    fn mito_engine_schema_version() -> *const c_char;
    fn mito_engine_error_schema_version() -> *const c_char;
    fn mito_engine_analyze_with_config_v2(
        engine: *mut c_void,
        input_path: *const c_char,
        ref_path: *const c_char,
        filter_numt: bool,
        threads: usize,
        min_mapping_quality: u8,
        min_base_quality: u8,
        excluded_snp_flags: u16,
        numt_threshold: f64,
        allow_development_tags: bool,
        should_cancel: Option<unsafe extern "C" fn(*mut c_void) -> bool>,
        cancel_user_data: *mut c_void,
    ) -> *const c_char;
    fn mito_engine_get_last_error() -> *const c_char;
    fn mito_engine_get_last_error_code() -> *const c_char;
    fn mito_engine_free_string(value: *const c_char);
}

pub struct MitoEngine {
    raw: *mut c_void,
}

// SAFETY: AnalysisEngine owns no thread-affine state. The Rust wrapper is not Sync, so one engine
// cannot be called concurrently through safe Rust, and destruction occurs on the receiving thread.
unsafe impl Send for MitoEngine {}

impl MitoEngine {
    pub fn capabilities() -> EngineCapabilities {
        EngineCapabilities {
            engine_version: static_c_string(unsafe { mito_engine_version() }),
            schema_version: static_c_string(unsafe { mito_engine_schema_version() }),
            error_schema_version: static_c_string(unsafe { mito_engine_error_schema_version() }),
            htslib: unsafe { mito_engine_has_htslib() },
        }
    }

    pub fn new() -> Result<Self, MitoError> {
        let raw = unsafe { mito_engine_new() };
        if raw.is_null() {
            return Err(last_error());
        }
        Ok(Self { raw })
    }

    pub fn analyze(
        &self,
        input_path: &Path,
        reference_path: Option<&Path>,
    ) -> Result<String, MitoError> {
        self.analyze_with_options(input_path, reference_path, AnalyzeOptions::default())
    }

    pub fn analyze_with_options(
        &self,
        input_path: &Path,
        reference_path: Option<&Path>,
        options: AnalyzeOptions,
    ) -> Result<String, MitoError> {
        self.analyze_with_cancellation(input_path, reference_path, options, None)
    }

    pub fn analyze_with_cancel_flag(
        &self,
        input_path: &Path,
        reference_path: Option<&Path>,
        options: AnalyzeOptions,
        cancel_flag: &AtomicBool,
    ) -> Result<String, MitoError> {
        self.analyze_with_cancellation(input_path, reference_path, options, Some(cancel_flag))
    }

    fn analyze_with_cancellation(
        &self,
        input_path: &Path,
        reference_path: Option<&Path>,
        options: AnalyzeOptions,
        cancel_flag: Option<&AtomicBool>,
    ) -> Result<String, MitoError> {
        let input = path_to_cstring(input_path, "input")?;
        let reference = match reference_path {
            Some(path) => Some(path_to_cstring(path, "reference")?),
            None => None,
        };

        let reference_ptr = reference
            .as_ref()
            .map_or(std::ptr::null(), |value| value.as_ptr());
        let (callback, user_data) = cancel_flag.map_or((None, std::ptr::null_mut()), |flag| {
            (
                Some(cancel_callback as unsafe extern "C" fn(*mut c_void) -> bool),
                flag as *const AtomicBool as *mut c_void,
            )
        });
        let result = unsafe {
            mito_engine_analyze_with_config_v2(
                self.raw,
                input.as_ptr(),
                reference_ptr,
                options.filter_numt,
                options.threads.max(1),
                options.min_mapping_quality,
                options.min_base_quality,
                options.excluded_snp_flags,
                options.numt_threshold,
                options.allow_development_tags,
                callback,
                user_data,
            )
        };

        if result.is_null() {
            return Err(last_error());
        }

        let json = unsafe { CStr::from_ptr(result) }
            .to_string_lossy()
            .into_owned();
        unsafe { mito_engine_free_string(result) };
        Ok(json)
    }
}

fn static_c_string(value: *const c_char) -> String {
    if value.is_null() {
        return "unknown".to_owned();
    }
    unsafe { CStr::from_ptr(value) }
        .to_string_lossy()
        .into_owned()
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EngineCapabilities {
    pub engine_version: String,
    pub schema_version: String,
    pub error_schema_version: String,
    pub htslib: bool,
}

fn path_to_cstring(path: &Path, label: &str) -> Result<CString, MitoError> {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        CString::new(path.as_os_str().as_bytes()).map_err(|_| {
            MitoError::new(
                "MITO-E1001",
                format!("{label} path contains an interior NUL byte"),
            )
        })
    }
    #[cfg(not(unix))]
    {
        let value = path.to_str().ok_or_else(|| {
            MitoError::new("MITO-E1001", format!("{label} path is not valid UTF-8"))
        })?;
        CString::new(value).map_err(|_| {
            MitoError::new(
                "MITO-E1001",
                format!("{label} path contains an interior NUL byte"),
            )
        })
    }
}

unsafe extern "C" fn cancel_callback(user_data: *mut c_void) -> bool {
    if user_data.is_null() {
        return false;
    }
    let flag = &*(user_data as *const AtomicBool);
    flag.load(Ordering::Relaxed)
}

impl Drop for MitoEngine {
    fn drop(&mut self) {
        if !self.raw.is_null() {
            unsafe { mito_engine_delete(self.raw) };
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct AnalyzeOptions {
    pub filter_numt: bool,
    pub threads: usize,
    pub min_mapping_quality: u8,
    pub min_base_quality: u8,
    pub excluded_snp_flags: u16,
    pub numt_threshold: f64,
    pub allow_development_tags: bool,
}

impl Default for AnalyzeOptions {
    fn default() -> Self {
        Self {
            filter_numt: true,
            threads: 1,
            min_mapping_quality: 20,
            min_base_quality: 10,
            excluded_snp_flags: 0xF00,
            numt_threshold: 0.30,
            allow_development_tags: false,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MitoError {
    pub code: String,
    pub message: String,
}

impl MitoError {
    pub fn new(code: impl Into<String>, message: impl Into<String>) -> Self {
        Self {
            code: code.into(),
            message: message.into(),
        }
    }
}

impl fmt::Display for MitoError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "[{}] {}", self.code, self.message)
    }
}

impl std::error::Error for MitoError {}

fn last_error() -> MitoError {
    let code = static_c_string(unsafe { mito_engine_get_last_error_code() });
    let message = static_c_string(unsafe { mito_engine_get_last_error() });
    MitoError::new(
        if code == "unknown" || code.is_empty() {
            "MITO-E9001".to_owned()
        } else {
            code
        },
        if message == "unknown" || message.is_empty() {
            "unknown mito engine error".to_owned()
        } else {
            message
        },
    )
}

#[cfg(test)]
mod tests {
    use super::{AnalyzeOptions, MitoEngine};
    use std::fs;
    use std::io::Write;
    use std::sync::atomic::{AtomicBool, Ordering};

    #[test]
    fn cancelled_analysis_returns_error() {
        let input_path = std::env::temp_dir().join(format!(
            "mito_ffi_cancel_{}_{}.fastq",
            std::process::id(),
            unique_suffix()
        ));
        let mut input = fs::File::create(&input_path).expect("temp input");
        writeln!(input, "@read-1").expect("write header");
        writeln!(input, "GATCACAGGT").expect("write sequence");
        writeln!(input, "+").expect("write plus");
        writeln!(input, "IIIIIIIIII").expect("write quality");
        drop(input);

        let cancel = AtomicBool::new(true);
        let engine = MitoEngine::new().expect("engine");
        let error = engine
            .analyze_with_cancel_flag(
                &input_path,
                None,
                AnalyzeOptions {
                    filter_numt: true,
                    threads: 1,
                    ..AnalyzeOptions::default()
                },
                &cancel,
            )
            .expect_err("analysis should be cancelled");
        assert_eq!(error.code, "MITO-E1501");
        assert_eq!(error.message, "analysis cancelled");
        assert!(cancel.load(Ordering::Relaxed));
        let _ = fs::remove_file(input_path);
    }

    #[test]
    fn external_failures_preserve_stable_error_codes() {
        let valid_input = temp_path("valid", "fastq");
        fs::write(&valid_input, "@read-1\nGATCACAGGT\n+\nIIIIIIIIII\n").expect("valid input");
        let malformed_input = temp_path("malformed", "fastq");
        fs::write(&malformed_input, "@read-1\nGATC\nnot-plus\nIIII\n").expect("malformed input");
        let fasta_input = temp_path("reads", "fasta");
        fs::write(&fasta_input, ">read-1\nGATC\n").expect("FASTA input");
        let invalid_reference = temp_path("invalid-reference", "fasta");
        fs::write(&invalid_reference, "GATC\n").expect("invalid reference");
        let empty_input = temp_path("empty", "fastq");
        fs::write(&empty_input, "").expect("empty input");
        let missing_input = temp_path("missing", "fastq");
        let missing_reference = temp_path("missing-reference", "fasta");

        let engine = MitoEngine::new().expect("engine");
        let cases = [
            (
                &malformed_input,
                None,
                AnalyzeOptions::default(),
                "MITO-E1103",
            ),
            (&fasta_input, None, AnalyzeOptions::default(), "MITO-E1102"),
            (&empty_input, None, AnalyzeOptions::default(), "MITO-E1104"),
            (
                &missing_input,
                None,
                AnalyzeOptions::default(),
                "MITO-E1101",
            ),
            (
                &valid_input,
                Some(invalid_reference.as_path()),
                AnalyzeOptions::default(),
                "MITO-E1202",
            ),
            (
                &valid_input,
                Some(missing_reference.as_path()),
                AnalyzeOptions::default(),
                "MITO-E1201",
            ),
            (
                &valid_input,
                None,
                AnalyzeOptions {
                    numt_threshold: f64::NAN,
                    ..AnalyzeOptions::default()
                },
                "MITO-E1001",
            ),
        ];
        for (input, reference, options, expected_code) in cases {
            let error = engine
                .analyze_with_options(input, reference, options)
                .expect_err("negative case must fail");
            assert_eq!(error.code, expected_code, "unexpected error: {error}");
        }

        for path in [
            valid_input,
            malformed_input,
            fasta_input,
            invalid_reference,
            empty_input,
        ] {
            let _ = fs::remove_file(path);
        }
    }

    fn temp_path(label: &str, extension: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "mito_ffi_{label}_{}_{}.{}",
            std::process::id(),
            unique_suffix(),
            extension
        ))
    }

    fn unique_suffix() -> u128 {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|duration| duration.as_nanos())
            .unwrap_or_default()
    }
}
