#!/usr/bin/env bash
set -euo pipefail

BUILD_TYPE="${1:-RelWithDebInfo}"
BUILD_DIR="${2:-build-prof}"

cmake -S . -B "$BUILD_DIR" -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE="$BUILD_TYPE" -DPROFILE_MODE=ON
cmake --build "$BUILD_DIR" --target vmc
"./$BUILD_DIR/vmc"
