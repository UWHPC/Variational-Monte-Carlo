#!/usr/bin/env bash
set -euo pipefail

cmake --build build
cmake -S . -B build
