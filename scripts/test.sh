#!/usr/bin/env bash
set -euo pipefail

BUILD_TYPE="Debug"
FP32=0

for arg in "$@"; do
  case "$arg" in
    --fp32)
      FP32=1
      ;;
    *)
      BUILD_TYPE="$arg"
      ;;
  esac
done

if [ "$FP32" -eq 1 ]; then
  FP64_FLAG="OFF"
else
  FP64_FLAG="ON"
fi

case "$(uname -s)" in
  MINGW*|MSYS*|CYGWIN*)
    echo "MSVC builds on Windows must run from PowerShell (sourcing vcvars64.bat from Git Bash hangs)." >&2
    echo "Run instead: .\\scripts\\test.ps1 $([ "$FP32" -eq 1 ] && echo "-Fp32 ")-BuildType $BUILD_TYPE" >&2
    exit 1
    ;;
  *)
    cmake -S . -B build-tests -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
          -DFP_64=$FP64_FLAG -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build build-tests --target vmc_tests
    ctest --test-dir build-tests --output-on-failure
    ;;
esac
