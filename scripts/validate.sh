#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cxx_compiler="/usr/bin/g++-15"
cuda_compiler="/usr/local/cuda-13.3/bin/nvcc"

run_cpu=0
run_cuda=0
backend_selected=0
fp64_flag="ON"
fast_math_flag="OFF"

for argument in "$@"; do
  case "$argument" in
    --cpu)
      run_cpu=1
      backend_selected=1
      ;;
    --cuda)
      run_cuda=1
      backend_selected=1
      ;;
    --fp32)
      fp64_flag="OFF"
      ;;
    --ffast-math)
      fast_math_flag="ON"
      ;;
    *)
      echo "unknown argument: $argument" >&2
      exit 2
      ;;
  esac
done

if [ "$backend_selected" -eq 0 ]; then
  run_cpu=1
  run_cuda=1
fi

if [ ! -x "$cxx_compiler" ]; then
  echo "missing C++ compiler: $cxx_compiler" >&2
  exit 1
fi

cuda_toolchain_available() {
  if [ ! -x "$cuda_compiler" ]; then
    cuda_unavailable_reason="missing CUDA compiler: $cuda_compiler"
    return 1
  fi

  local cuda_root
  cuda_root="$(dirname "$(dirname "$cuda_compiler")")"
  if [ ! -e "$cuda_root/targets/x86_64-linux/lib/libcusolver.so" ] \
     && [ ! -e "$cuda_root/lib64/libcusolver.so" ]; then
    cuda_unavailable_reason="CUDA toolkit is missing cuSOLVER libraries under $cuda_root"
    return 1
  fi

  return 0
}

run_pipeline() {
  local backend="$1"
  local build_dir="$repo_root/build-validation-$backend"
  local cuda_flag="OFF"
  local cmake_arguments=(
    -S "$repo_root"
    -B "$build_dir"
    -G Ninja
    -DCMAKE_CXX_COMPILER="$cxx_compiler"
    -DFP_64="$fp64_flag"
    -DVMC_FAST_MATH="$fast_math_flag"
    -DBUILD_TESTING=OFF
    -DVMC_BUILD_VALIDATION=ON
    -DCMAKE_BUILD_TYPE=Release
  )

  local cached_dependency_root="$repo_root/build-test-cpu/_deps"
  if [ -d "$cached_dependency_root/xpu-src" ]; then
    cmake_arguments+=(
      -DFETCHCONTENT_SOURCE_DIR_XPU="$cached_dependency_root/xpu-src"
    )
  fi
  if [ -d "$cached_dependency_root/catch2-src" ]; then
    cmake_arguments+=(
      -DFETCHCONTENT_SOURCE_DIR_CATCH2="$cached_dependency_root/catch2-src"
    )
  fi

  if [ -d "$build_dir/_deps/xpu-src" ] && [ -d "$build_dir/_deps/catch2-src" ]; then
    cmake_arguments+=(-DFETCHCONTENT_UPDATES_DISCONNECTED=ON)
  fi

  if [ "$backend" = "cuda" ]; then
    cuda_flag="ON"
    cmake_arguments+=(
      -DCMAKE_CUDA_COMPILER="$cuda_compiler"
      -DCMAKE_CUDA_HOST_COMPILER="$cxx_compiler"
    )
  fi
  cmake_arguments+=("-DVMC_ENABLE_CUDA=$cuda_flag")

  echo "==> configuring $backend validation"
  cmake "${cmake_arguments[@]}"

  echo "==> compiling $backend validation"
  cmake --build "$build_dir" --target vmc_validations

  if [ "$backend" = "cuda" ]; then
    if ! command -v nvidia-smi >/dev/null 2>&1 || ! nvidia-smi -L >/dev/null 2>&1; then
      echo "==> CUDA validation compiled; no NVIDIA GPU detected, skipping execution"
      return
    fi
  fi

  echo "==> running $backend validation"
  ctest --test-dir "$build_dir" --output-on-failure
}

if [ "$run_cpu" -eq 1 ]; then
  run_pipeline cpu
fi

if [ "$run_cuda" -eq 1 ]; then
  if ! cuda_toolchain_available; then
    if [ "$backend_selected" -eq 1 ]; then
      echo "$cuda_unavailable_reason" >&2
      exit 1
    fi
    echo "==> $cuda_unavailable_reason; skipping CUDA validation"
  else
    run_pipeline cuda
  fi
fi
