#!/usr/bin/env bash
set -euo pipefail

BUILD_TYPE="${1:-Release}"

cmake -S . -B build -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
cmake --build build --target vmc