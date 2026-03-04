#!/usr/bin/env bash
set -euo pipefail

cmake -S . -B build -DBUILD_TESTING=OFF
cmake --build build --target vmc
