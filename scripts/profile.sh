#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR=""
CUDA=0
FP32=0

for arg in "$@"; do
  case "$arg" in
    --cuda)
      CUDA=1
      ;;
    --fp32)
      FP32=1
      ;;
    *)
      BUILD_DIR="$arg"
      ;;
  esac
done

if [ "$CUDA" -eq 1 ]; then
  CUDA_FLAG="ON"
  DEFAULT_DIR="build-prof-cuda"
else
  CUDA_FLAG="OFF"
  DEFAULT_DIR="build-prof"
fi
BUILD_DIR="${BUILD_DIR:-$DEFAULT_DIR}"

if [ "$FP32" -eq 1 ]; then
  FP64_FLAG="OFF"
else
  FP64_FLAG="ON"
fi

case "$(uname -s)" in
  MINGW*|MSYS*|CYGWIN*)
    echo "CUDA/MSVC builds on Windows must run from PowerShell (sourcing vcvars64.bat from Git Bash hangs)." >&2
    echo "Run instead: .\\scripts\\profile.ps1 $([ "$CUDA" -eq 1 ] && echo "-Cuda ")$([ "$FP32" -eq 1 ] && echo "-Fp32 ")-BuildDir $BUILD_DIR" >&2
    exit 1
    ;;
  *)
    cmake -S . -B "$BUILD_DIR" -DVMC_ENABLE_CUDA=$CUDA_FLAG -DFP_64=$FP64_FLAG -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo
    cmake --build "$BUILD_DIR" --target vmc
    ;;
esac

echo "Profiler build ready: ./$BUILD_DIR/vmc"
