#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${MITO_SANITIZER_BUILD_DIR:-$ROOT_DIR/build-sanitize}"

cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DMITO_BUILD_TESTS=ON \
  -DMITO_ENABLE_SANITIZERS=ON
cmake --build "$BUILD_DIR" -j"${MITO_BUILD_JOBS:-2}"
ASAN_OPTIONS="detect_leaks=${MITO_ASAN_DETECT_LEAKS:-0}:strict_string_checks=1" \
UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=1" \
  ctest --test-dir "$BUILD_DIR" --output-on-failure
