#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CXX_COMPILER="/usr/bin/g++-15"
CUDA_COMPILER="/usr/local/cuda-13.3/bin/nvcc"

if [ ! -x "$CXX_COMPILER" ]; then
  echo "missing clangd host compiler: $CXX_COMPILER" >&2
  exit 1
fi

if [ ! -x "$CUDA_COMPILER" ]; then
  echo "missing clangd CUDA compiler: $CUDA_COMPILER" >&2
  exit 1
fi

cmake -S "$REPO_ROOT" -B "$REPO_ROOT/build-cpp" -G Ninja \
  -DVMC_ENABLE_CUDA=OFF \
  -DFP_64=ON \
  -DVMC_FAST_MATH=OFF \
  -DBUILD_TESTING=ON \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER="$CXX_COMPILER"

cmake -S "$REPO_ROOT" -B "$REPO_ROOT/build-cuda" -G Ninja \
  -DVMC_ENABLE_CUDA=ON \
  -DFP_64=ON \
  -DVMC_FAST_MATH=OFF \
  -DBUILD_TESTING=ON \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER="$CXX_COMPILER" \
  -DCMAKE_CUDA_COMPILER="$CUDA_COMPILER" \
  -DCMAKE_CUDA_HOST_COMPILER="$CXX_COMPILER" \
  -DCMAKE_CUDA_ARCHITECTURES=86

echo "clangd databases configured in build-cpp and build-cuda"
