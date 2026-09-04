#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CXX_COMPILER="/usr/bin/g++-15"
CUDA_COMPILER="/usr/local/cuda-13.3/bin/nvcc"

CUDA=0
FP64_FLAG="ON"
FAST_MATH_FLAG="OFF"

for arg in "$@"; do
  case "$arg" in
    --cuda)
      CUDA=1
      ;;
    --fp32)
      FP64_FLAG="OFF"
      ;;
    --ffast-math)
      FAST_MATH_FLAG="ON"
      ;;
    *)
      echo "unknown argument: $arg" >&2
      exit 2
      ;;
  esac
done

if [ ! -x "$CXX_COMPILER" ]; then
  echo "missing C++ compiler: $CXX_COMPILER" >&2
  exit 1
fi

if [ "$CUDA" -eq 1 ]; then
  if [ ! -x "$CUDA_COMPILER" ]; then
    echo "missing CUDA compiler: $CUDA_COMPILER" >&2
    exit 1
  fi

  BUILD_DIR="build-debug-cuda"
  CUDA_FLAG="ON"
else
  BUILD_DIR="build-debug"
  CUDA_FLAG="OFF"
fi

CMAKE_ARGS=(
  -S "$REPO_ROOT"
  -B "$REPO_ROOT/$BUILD_DIR"
  -G Ninja
  -DCMAKE_CXX_COMPILER="$CXX_COMPILER"
  -DVMC_ENABLE_CUDA="$CUDA_FLAG"
  -DFP_64="$FP64_FLAG"
  -DVMC_FAST_MATH="$FAST_MATH_FLAG"
  -DBUILD_TESTING=OFF
  -DCMAKE_BUILD_TYPE=Debug
)

if [ "$CUDA" -eq 1 ]; then
  CMAKE_ARGS+=(
    -DCMAKE_CUDA_COMPILER="$CUDA_COMPILER"
    -DCMAKE_CUDA_HOST_COMPILER="$CXX_COMPILER"
  )
fi

cmake "${CMAKE_ARGS[@]}"
cmake --build "$REPO_ROOT/$BUILD_DIR" --target vmc

cd "$REPO_ROOT"
exec "./$BUILD_DIR/vmc"
