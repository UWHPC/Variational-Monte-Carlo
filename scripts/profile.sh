#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="${1:-build-prof}"

cmake -S . -B "$BUILD_DIR" -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build "$BUILD_DIR" --target vmc
echo "Profiler build ready: ./$BUILD_DIR/vmc"
