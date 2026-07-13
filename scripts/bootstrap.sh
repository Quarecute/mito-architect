#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE="check"
RUN_VERIFY=0
ASSUME_YES=0
INCLUDE_DOCKER=0

usage() {
  cat <<'EOF'
Usage: bash scripts/bootstrap.sh [options]

Install or audit the native Mito-Architect development environment.

Options:
  --check             Audit dependencies without changing the host (default).
  --install           Install system and locked project dependencies.
  --verify            Run the canonical verification gate after the audit/install.
  --include-docker    Audit/install Docker in addition to the native toolchain.
  --yes               Pass the package manager's non-interactive confirmation flag.
  -h, --help          Show this help.

Supported installers: Debian/Ubuntu (apt-get) and Arch Linux (pacman plus
yay/paru for bioinformatics packages when they are not in configured repos).
Other Linux distributions can always use --check and follow the reported gaps.
EOF
}

while (($#)); do
  case "$1" in
    --check) MODE="check" ;;
    --install) MODE="install" ;;
    --verify) RUN_VERIFY=1 ;;
    --include-docker) INCLUDE_DOCKER=1 ;;
    --yes) ASSUME_YES=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "error: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

log() {
  printf '[bootstrap] %s\n' "$*"
}

run_privileged() {
  if [[ "${EUID:-$(id -u)}" -eq 0 ]]; then
    "$@"
  elif command -v sudo >/dev/null 2>&1; then
    sudo "$@"
  else
    echo "error: system package installation needs root or sudo" >&2
    return 1
  fi
}

install_arch() {
  local pacman_args=(-S --needed)
  local aur_args=(-S --needed)
  if [[ "$ASSUME_YES" -eq 1 ]]; then
    pacman_args+=(--noconfirm)
    aur_args+=(--noconfirm)
  fi

  local official=(base-devel cmake ninja pkgconf nodejs npm curl gzip jq ca-certificates)
  if ! command -v cargo >/dev/null 2>&1 || ! command -v rustc >/dev/null 2>&1; then
    official+=(rust)
  fi
  if [[ "$INCLUDE_DOCKER" -eq 1 ]]; then
    official+=(docker docker-compose)
  fi
  run_privileged pacman "${pacman_args[@]}" "${official[@]}"

  local bio=(htslib samtools minimap2 bcftools sra-tools)
  local helper=""
  if command -v yay >/dev/null 2>&1; then
    helper="yay"
  elif command -v paru >/dev/null 2>&1; then
    helper="paru"
  fi
  if [[ -n "$helper" ]]; then
    "$helper" "${aur_args[@]}" "${bio[@]}"
  else
    echo "error: Arch bioinformatics packages require yay/paru or equivalent manual installation: ${bio[*]}" >&2
    return 1
  fi
}

install_debian() {
  local apt_args=(install --no-install-recommends)
  if [[ "$ASSUME_YES" -eq 1 ]]; then
    apt_args+=(-y)
  fi
  local packages=(
    build-essential ca-certificates cmake ninja-build pkg-config libhts-dev
    samtools minimap2 bcftools sra-toolkit nodejs npm rustc cargo curl gzip jq
  )
  if [[ "$INCLUDE_DOCKER" -eq 1 ]]; then
    packages+=(docker.io docker-compose-plugin)
  fi
  run_privileged apt-get update
  run_privileged apt-get "${apt_args[@]}" "${packages[@]}"
}

install_system_dependencies() {
  if [[ ! -r /etc/os-release ]]; then
    echo "error: cannot identify the operating system; use --check" >&2
    return 1
  fi
  # shellcheck disable=SC1091
  source /etc/os-release
  case "${ID:-}:${ID_LIKE:-}" in
    arch:*|*:arch*) install_arch ;;
    debian:*|ubuntu:*|*:debian*) install_debian ;;
    *)
      echo "error: automatic installation is not supported for ${PRETTY_NAME:-this OS}; use --check" >&2
      return 1
      ;;
  esac
}

version_ge() {
  [[ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n1)" == "$2" ]]
}

declare -a ERRORS=()

require_command() {
  local command_name="$1"
  local purpose="$2"
  if command -v "$command_name" >/dev/null 2>&1; then
    log "ok: $command_name ($(command -v "$command_name"))"
  else
    ERRORS+=("missing command '$command_name' ($purpose)")
  fi
}

require_minimum_version() {
  local command_name="$1"
  local actual="$2"
  local minimum="$3"
  if version_ge "$actual" "$minimum"; then
    log "ok: $command_name $actual (minimum $minimum)"
  else
    ERRORS+=("$command_name $actual is older than required $minimum")
  fi
}

audit_dependencies() {
  ERRORS=()
  local commands=(
    "c++:C++20 compiler"
    "cmake:C++ configuration"
    "pkg-config:native library discovery"
    "cargo:Rust dependency/build driver"
    "rustc:Rust compiler"
    "node:frontend runtime"
    "npm:locked frontend dependency driver"
    "samtools:BAM/CRAM workflow"
    "minimap2:competitive long-read mapping"
    "bcftools:VCF interoperability"
    "prefetch:bounded public fixture retrieval"
    "fasterq-dump:public FASTQ conversion"
    "curl:resource retrieval"
    "gzip:compressed clinical resources"
    "jq:deterministic verification support"
    "sha256sum:resource integrity verification"
  )
  if [[ "$INCLUDE_DOCKER" -eq 1 ]]; then
    commands+=("docker:container build/deployment")
  fi

  local item
  for item in "${commands[@]}"; do
    require_command "${item%%:*}" "${item#*:}"
  done

  if command -v cmake >/dev/null 2>&1; then
    require_minimum_version "cmake" "$(cmake --version | awk 'NR == 1 {print $3}')" "3.24"
  fi
  if command -v rustc >/dev/null 2>&1; then
    require_minimum_version "rustc" "$(rustc --version | awk '{print $2}')" "1.80"
  fi
  if command -v node >/dev/null 2>&1; then
    require_minimum_version "node" "$(node --version | sed 's/^v//')" "18.0"
  fi

  if pkg-config --exists htslib 2>/dev/null; then
    log "ok: htslib $(pkg-config --modversion htslib)"
  else
    ERRORS+=("htslib development package is missing from pkg-config")
  fi

  if [[ "${#ERRORS[@]}" -ne 0 ]]; then
    printf 'dependency audit failed:\n' >&2
    printf '  - %s\n' "${ERRORS[@]}" >&2
    printf 'Install a newer Rust or Node toolchain if the distribution package is below the reported minimum.\n' >&2
    return 1
  fi
  log "native dependency audit passed"
}

install_project_dependencies() {
  cd "$ROOT_DIR"
  log "resolving locked Rust dependencies"
  cargo fetch --locked
  log "installing locked JavaScript dependencies"
  npm ci

  # A prior build without htslib can otherwise remain cached after libhts is
  # installed. Cleaning only the bridge is deliberate and bounded.
  cargo clean -p mito-ffi
  log "building the workspace with the detected native capabilities"
  cargo build --workspace --locked
  npm run build
}

cd "$ROOT_DIR"
if [[ "$MODE" == "install" ]]; then
  log "installing supported host packages"
  install_system_dependencies
fi

audit_dependencies

if [[ "$MODE" == "install" ]]; then
  install_project_dependencies
  cargo run -p mito-cli --locked -- doctor
fi

if [[ "$RUN_VERIFY" -eq 1 ]]; then
  bash scripts/verify.sh
fi

log "complete"
