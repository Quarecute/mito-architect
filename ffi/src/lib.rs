use anyhow::{anyhow, Context, Result};
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_void};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

unsafe extern "C" {
    fn mito_engine_new() -> *mut c_void;
    fn mito_engine_delete(engine: *mut c_void);
    fn mito_engine_analyze_with_options(
        engine: *mut c_void,
        input_path: *const c_char,
        ref_path: *const c_char,
        filter_numt: bool,
        threads: usize,
    ) -> *const c_char;
    fn mito_engine_analyze_with_cancel(
        engine: *mut c_void,
        input_path: *const c_char,
        ref_path: *const c_char,
        filter_numt: bool,
        threads: usize,
        should_cancel: Option<unsafe extern "C" fn(*mut c_void) -> bool>,
        cancel_user_data: *mut c_void,
    ) -> *const c_char;
    fn mito_engine_get_last_error() -> *const c_char;
    fn mito_engine_free_string(value: *const c_char);
}

pub struct MitoEngine {
    raw: *mut c_void,
}

// SAFETY: AnalysisEngine owns no thread-affine state. The Rust wrapper is not Sync, so one engine
// cannot be called concurrently through safe Rust, and destruction occurs on the receiving thread.
unsafe impl Send for MitoEngine {}

impl MitoEngine {
    pub fn new() -> Result<Self> {
        let raw = unsafe { mito_engine_new() };
        if raw.is_null() {
            return Err(anyhow!(last_error()));
        }
        Ok(Self { raw })
    }

    pub fn analyze(&self, input_path: &Path, reference_path: Option<&Path>) -> Result<String> {
        self.analyze_with_options(input_path, reference_path, AnalyzeOptions::default())
    }

    pub fn analyze_with_options(
        &self,
        input_path: &Path,
        reference_path: Option<&Path>,
        options: AnalyzeOptions,
    ) -> Result<String> {
        self.analyze_with_cancellation(input_path, reference_path, options, None)
    }

    pub fn analyze_with_cancel_flag(
        &self,
        input_path: &Path,
        reference_path: Option<&Path>,
        options: AnalyzeOptions,
        cancel_flag: &AtomicBool,
    ) -> Result<String> {
        self.analyze_with_cancellation(input_path, reference_path, options, Some(cancel_flag))
    }

    fn analyze_with_cancellation(
        &self,
        input_path: &Path,
        reference_path: Option<&Path>,
        options: AnalyzeOptions,
        cancel_flag: Option<&AtomicBool>,
    ) -> Result<String> {
        let input = path_to_cstring(input_path, "input")?;
        let reference = match reference_path {
            Some(path) => Some(path_to_cstring(path, "reference")?),
            None => None,
        };

        let reference_ptr = reference
            .as_ref()
            .map_or(std::ptr::null(), |value| value.as_ptr());
        let result = if let Some(flag) = cancel_flag {
            unsafe {
                mito_engine_analyze_with_cancel(
                    self.raw,
                    input.as_ptr(),
                    reference_ptr,
                    options.filter_numt,
                    options.threads.max(1),
                    Some(cancel_callback),
                    flag as *const AtomicBool as *mut c_void,
                )
            }
        } else {
            unsafe {
                mito_engine_analyze_with_options(
                    self.raw,
                    input.as_ptr(),
                    reference_ptr,
                    options.filter_numt,
                    options.threads.max(1),
                )
            }
        };

        if result.is_null() {
            return Err(anyhow!(last_error()));
        }

        let json = unsafe { CStr::from_ptr(result) }
            .to_string_lossy()
            .into_owned();
        unsafe { mito_engine_free_string(result) };
        Ok(json)
    }
}

fn path_to_cstring(path: &Path, label: &str) -> Result<CString> {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        CString::new(path.as_os_str().as_bytes())
            .with_context(|| format!("{label} path contains an interior NUL byte"))
    }
    #[cfg(not(unix))]
    {
        let value = path
            .to_str()
            .with_context(|| format!("{label} path is not valid UTF-8"))?;
        CString::new(value).with_context(|| format!("{label} path contains an interior NUL byte"))
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
}

impl Default for AnalyzeOptions {
    fn default() -> Self {
        Self {
            filter_numt: true,
            threads: 1,
        }
    }
}

fn last_error() -> String {
    let ptr = unsafe { mito_engine_get_last_error() };
    if ptr.is_null() {
        return "unknown mito engine error".to_string();
    }
    let message = unsafe { CStr::from_ptr(ptr) }.to_string_lossy();
    if message.is_empty() {
        "unknown mito engine error".to_string()
    } else {
        message.into_owned()
    }
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
                },
                &cancel,
            )
            .expect_err("analysis should be cancelled");
        assert!(error.to_string().contains("analysis cancelled"));
        assert!(cancel.load(Ordering::Relaxed));
        let _ = fs::remove_file(input_path);
    }

    fn unique_suffix() -> u128 {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|duration| duration.as_nanos())
            .unwrap_or_default()
    }
}
