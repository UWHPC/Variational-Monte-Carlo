# Variational-Monte-Carlo

## Prerequisites

- Linux or WSL2
- CMake 3.25+
- Ninja
- GCC 15
- Git
- CUDA 13.3 for CUDA builds

## Run

Build and run an optimized CPU executable:

`./scripts/run.sh`

Available options:

- `--cuda` enables CUDA.
- `--fp32` uses single precision.
- `-ffast-math` enables unsafe floating-point optimizations.

Options can be combined, for example:

`./scripts/run.sh --cuda --fp32 -ffast-math`

## Debug

Build and run with debug instrumentation:

`./scripts/debug.sh`

`debug.sh` accepts the same options as `run.sh`.

## Run Tests

Build and run the test suite:

`./scripts/test.sh`

`test.sh` accepts the same options as `run.sh`.

## LSP / clangd

This repository uses separate CPU and CUDA compilation databases so clangd can
parse both regular C++ and CUDA translation units deterministically.

Generate them with:

`./scripts/configure-clangd.sh`
