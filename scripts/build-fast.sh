#!/usr/bin/env bash
set -euo pipefail

FP32=0
FAST_MATH=0
for arg in "$@"; do
  case "$arg" in
    --fp32)
      FP32=1
      ;;
    --vmc-fast-math)
      FAST_MATH=1
      ;;
  esac
done

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
    echo "MSVC builds on Windows must run from PowerShell (sourcing vcvars64.bat from Git Bash hangs)." >&2
    echo "Run instead: .\\scripts\\build.ps1 $([ "$FP32" -eq 1 ] && echo "-Fp32 ")$([ "$FAST_MATH" -eq 1 ] && echo "-Vmc-fast-math ")" >&2
    exit 1
    ;;
  *)
    cmake -S . -B build -DFP_64=$FP64_FLAG -DVMC_FAST_MATH=$FAST_MATH_FLAG -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release
    cmake --build build --target vmc
    ;;
esac
