use std::path::PathBuf;
use std::process::Command;
use std::{env, fs};

fn main() {
    println!("cargo:rerun-if-changed=../CMakeLists.txt");
    println!("cargo:rerun-if-changed=../core");
    println!("cargo:rerun-if-changed=include");
    println!("cargo:rerun-if-changed=src/mito_c_api.cpp");

    let profile = env::var("PROFILE").unwrap_or_else(|_| "debug".to_string());
    let mut config = cmake::Config::new("..");
    config.define("MITO_BUILD_TESTS", "OFF");
    config.profile(if profile == "release" {
        "Release"
    } else {
        "Debug"
    });
    println!("cargo:rerun-if-env-changed=MITO_HDBSCAN_CPP_ROOT");
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_PATH");
    if let Ok(output) = Command::new("pkg-config")
        .args(["--variable=pcfiledir", "htslib"])
        .output()
    {
        if output.status.success() {
            let directory = String::from_utf8_lossy(&output.stdout).trim().to_owned();
            if !directory.is_empty() {
                println!("cargo:rerun-if-changed={directory}/htslib.pc");
            }
        }
    }
    let hdbscan_root = env::var("MITO_HDBSCAN_CPP_ROOT").unwrap_or_default();
    config.define("MITO_HDBSCAN_CPP_ROOT", hdbscan_root);

    let dst = config.build();
    let lib_dir = PathBuf::from(&dst).join("lib");
    println!("cargo:rustc-link-search=native={}", lib_dir.display());
    println!("cargo:rustc-link-lib=static=mito_ffi");
    println!("cargo:rustc-link-lib=static=mito_core");
    if lib_dir.join("libmito_hdbscan_cpp.a").exists() {
        println!("cargo:rustc-link-lib=static=mito_hdbscan_cpp");
    }
    if !emit_cmake_cache_link_args(
        &PathBuf::from(&dst).join("build").join("CMakeCache.txt"),
        "HTSLIB",
    ) {
        emit_pkg_config_link_args("htslib");
    }

    if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=dylib=c++");
    } else if cfg!(target_env = "msvc") {
        println!("cargo:rustc-link-lib=dylib=msvcp140");
    } else {
        println!("cargo:rustc-link-lib=dylib=stdc++");
    }
}

fn emit_cmake_cache_link_args(cache_path: &PathBuf, prefix: &str) -> bool {
    let Ok(cache) = fs::read_to_string(cache_path) else {
        return false;
    };
    if cache_value(&cache, &format!("{prefix}_FOUND")).as_deref() != Some("1") {
        return false;
    }
    let libdir = cache_value(&cache, &format!("{prefix}_LIBDIR"));
    if let Some(libdir) = &libdir {
        if !libdir.is_empty() {
            println!("cargo:rustc-link-search=native={libdir}");
        }
    }
    if let Some(libraries) = cache_value(&cache, &format!("{prefix}_LIBRARIES")) {
        for lib in libraries.split(';').filter(|lib| !lib.is_empty()) {
            if lib == "hts"
                && libdir
                    .as_ref()
                    .is_some_and(|dir| PathBuf::from(dir).join("libhts.a").exists())
            {
                println!("cargo:rustc-link-lib=static=hts");
                for dep in ["deflate", "lzma", "bz2", "z"] {
                    println!("cargo:rustc-link-lib=dylib={dep}");
                }
            } else {
                println!("cargo:rustc-link-lib=dylib={lib}");
            }
        }
    }
    true
}

fn cache_value(cache: &str, key: &str) -> Option<String> {
    let prefix = format!("{key}:");
    cache.lines().find_map(|line| {
        line.strip_prefix(&prefix)
            .and_then(|rest| rest.split_once('='))
            .map(|(_, value)| value.to_string())
    })
}

fn emit_pkg_config_link_args(package: &str) {
    let search_paths = Command::new("pkg-config")
        .args(["--libs-only-L", package])
        .output();
    if let Ok(output) = search_paths {
        if output.status.success() {
            for token in String::from_utf8_lossy(&output.stdout).split_whitespace() {
                if let Some(path) = token.strip_prefix("-L") {
                    println!("cargo:rustc-link-search=native={path}");
                }
            }
        }
    }

    let libs = Command::new("pkg-config")
        .args(["--libs-only-l", package])
        .output();
    if let Ok(output) = libs {
        if output.status.success() {
            for token in String::from_utf8_lossy(&output.stdout).split_whitespace() {
                if let Some(lib) = token.strip_prefix("-l") {
                    println!("cargo:rustc-link-lib=dylib={lib}");
                }
            }
        }
    }
}
