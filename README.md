# Variational-Monte-Carlo

## Prerequisites

- Linux or WSL2
- CMake 3.25+
- Ninja
- GCC 15
- OpenBLAS and LAPACKE development files
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

By default, clangd uses the CPU compilation database so IntelliSense works with
or without CUDA. When CUDA 13.3 is installed, the setup script also generates a
separate CUDA compilation database in `build-cuda`.

When VS Code prompts for recommended extensions, install **clangd**. The
workspace disables the Microsoft C/C++ language server to prevent two
IntelliSense engines from competing.

Generate them with:

`./scripts/configure-clangd.sh`

The CPU database is generated even when CUDA is not installed.

After generating the database for the first time, run **clangd: Restart language
server** from the VS Code command palette (or reload the window).
