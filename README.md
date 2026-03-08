# Variational-Monte-Carlo

## Prerequisites

- CMake 3.20+
- A C++23 compiler
  - Windows: Visual Studio 2022 with "Desktop development with C++"
  - Linux: GCC 13+ or Clang 16+
  - macOS: AppleClang that supports C++23 features you use, or Homebrew LLVM
- Git
- clang-format
- clang-tidy

## Build

The easiest way to build is to execute: \
`./scripts/build.sh` or `./scripts/build.ps1`

Script parameters:
- `./scripts/build.sh [BUILD_TYPE]` (default: `Release`)
- `./scripts/build.ps1 [-BuildType <BUILD_TYPE>]` (default: `Release`)

Alternatively, a manual build involves: \
`cmake -S . -B build -DBUILD_TESTING=OFF` \
`cmake --build build --target vmc`

## Profiling Mode

You can enable a lightweight built-in profiler that:
- disables per-step JSON output (so I/O does not dominate samples),
- prints a timing breakdown for core simulation phases.

Build and run:
`./scripts/profile.sh` or `./scripts/profile.ps1`

Script parameters:
- `./scripts/profile.sh [BUILD_TYPE] [BUILD_DIR]` (defaults: `RelWithDebInfo`, `build-prof`)
- `./scripts/profile.ps1 [-BuildType <BUILD_TYPE>] [-BuildDir <BUILD_DIR>]`
  (defaults: `RelWithDebInfo`, `build-prof`)

Alternatively, manually:
`cmake -S . -B build-prof -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo -DPROFILE_MODE=ON`
`cmake --build build-prof --target vmc`
`./build-prof/vmc`

You will see `[profile]` lines for:
- warmup and measurement Metropolis time,
- `evaluate_log_psi`,
- derivative evaluation,
- local energy evaluation,
- blocking-statistics overhead,
- output writing overhead.

## Run Tests

The easiest way to run tests is to execute: \
`./scripts/test.sh` or `./scripts/test.ps1`

Script parameters:
- `./scripts/test.sh [BUILD_TYPE]` (default: `Debug`)
- `./scripts/test.ps1 [-BuildType <BUILD_TYPE>]` (default: `Debug`)

Alternatively, manually testing involves: \
`cmake -S . -B build-tests -DBUILD_TESTING=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON` \
`cmake --build build-tests --target vmc_tests` \
`ctest --test-dir build-tests`

## LSP / clangd

This repo includes a `.clangd` file that points clangd to the
`build-tests` compilation database so test files can resolve Catch2 headers.

Generate it once with:
`cmake -S . -B build-tests -DBUILD_TESTING=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON`
