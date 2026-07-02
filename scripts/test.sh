#!/usr/bin/env bash
set -euo pipefail

BUILD_TYPE="${1:-Debug}"

case "$(uname -s)" in
  MINGW*|MSYS*|CYGWIN*)
    echo "MSVC builds on Windows must run from PowerShell (sourcing vcvars64.bat from Git Bash hangs)." >&2
    echo "Run instead: .\\scripts\\test.ps1 $BUILD_TYPE" >&2
    exit 1
    ;;
  *)
    cmake -S . -B build-tests -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build build-tests --target vmc_tests
    ctest --test-dir build-tests --output-on-failure
    ;;
esac
