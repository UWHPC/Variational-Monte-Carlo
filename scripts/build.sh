#!/usr/bin/env bash
set -euo pipefail

BUILD_TYPE="Release"
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
      BUILD_TYPE="$arg"
      ;;
  esac
done

if [ "$CUDA" -eq 1 ]; then
  BUILD_DIR="build-cuda"
  CUDA_FLAG="ON"
else
  BUILD_DIR="build"
  CUDA_FLAG="OFF"
fi

if [ "$FP32" -eq 1 ]; then
  FP64_FLAG="OFF"
else
  FP64_FLAG="ON"
fi

case "$(uname -s)" in
  MINGW*|MSYS*|CYGWIN*)
    echo "CUDA/MSVC builds on Windows must run from PowerShell (sourcing vcvars64.bat from Git Bash hangs)." >&2
    echo "Run instead: .\\scripts\\build.ps1 $([ "$CUDA" -eq 1 ] && echo "-Cuda ")$([ "$FP32" -eq 1 ] && echo "-Fp32 ")-BuildType $BUILD_TYPE" >&2
    exit 1
    ;;
  *)
    cmake -S . -B "$BUILD_DIR" -DVMC_ENABLE_CUDA=$CUDA_FLAG -DFP_64=$FP64_FLAG -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
    cmake --build "$BUILD_DIR" --target vmc
    ;;
esac
