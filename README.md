# Variational-Monte-Carlo

## Prerequisites

- CMake 3.20+
- A C++23 compiler
  - Windows: Visual Studio 2022 with "Desktop development with C++"
  - Linux: GCC 13+ or Clang 16+
  - macOS: AppleClang that supports C++23 features you use, or Homebrew LLVM
- Git
- clang-format

## Build

The easiest way to build is to execute: \
`./scripts/build.sh` or `./scripts/build.ps1`

Alternatively, a manual build involves: \
`cmake -S . -B build -DBUILD_TESTING=OFF` \
`cmake --build build --target vmc`

## Run Tests

The easiest way to run tests is to execute: \
`./scripts/test.sh` or `./scripts/test.ps1`

Alternatively, manually testing involves: \
`cmake -S . -B build-tests -DBUILD_TESTING=ON` \
`cmake --build build-tests --target vmc_tests` \
`ctest --test-dir build-tests`
