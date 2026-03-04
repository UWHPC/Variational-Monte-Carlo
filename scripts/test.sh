#!/usr/bin/env bash
set -euo pipefail

cmake -S . -B build-tests -DBUILD_TESTING=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-tests --target vmc_tests
ctest --test-dir build-tests --output-on-failure
