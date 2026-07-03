#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR=""
CUDA=0
FP32=0
FAST_MATH=0

for arg in "$@"; do
  case "$arg" in
    --cuda)
      CUDA=1
      ;;
    --fp32)
      FP32=1
      ;;
    --vmc-fast-math)
      FAST_MATH=1
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

if [ "$FP32" -eq 1 ] || [ "$FAST_MATH" -eq 1 ]; then
  FP64_FLAG="OFF"
else
  FP64_FLAG="ON"
fi

if [ "$FAST_MATH" -eq 1 ]; then
  FAST_MATH_FLAG="ON"
else
  FAST_MATH_FLAG="OFF"
fi

case "$(uname -s)" in
  MINGW*|MSYS*|CYGWIN*)
    echo "CUDA/MSVC builds on Windows must run from PowerShell (sourcing vcvars64.bat from Git Bash hangs)." >&2
    echo "Run instead: .\\scripts\\profile.ps1 $([ "$CUDA" -eq 1 ] && echo "-Cuda ")$([ "$FP32" -eq 1 ] && echo "-Fp32 ")$([ "$FAST_MATH" -eq 1 ] && echo "-Vmc-fast-math ")-BuildDir $BUILD_DIR" >&2
    exit 1
    ;;
  *)
    cmake -S . -B "$BUILD_DIR" -DVMC_ENABLE_CUDA=$CUDA_FLAG -DFP_64=$FP64_FLAG -DVMC_FAST_MATH=$FAST_MATH_FLAG -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo
    cmake --build "$BUILD_DIR" --target vmc
    ;;
esac

echo "Profiler build ready: ./$BUILD_DIR/vmc"
