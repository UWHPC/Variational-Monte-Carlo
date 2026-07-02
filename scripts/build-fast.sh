#!/usr/bin/env bash
set -euo pipefail

FP32=0
for arg in "$@"; do
  case "$arg" in
    --fp32)
      FP32=1
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
    echo "Run instead: .\\scripts\\build.ps1 $([ "$FP32" -eq 1 ] && echo "-Fp32 ")" >&2
    exit 1
    ;;
  *)
    cmake -S . -B build -DFP_64=$FP64_FLAG -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release
    cmake --build build --target vmc
    ;;
esac
