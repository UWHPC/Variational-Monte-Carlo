# Variational-Monte-Carlo

## Prerequisites

- CMake 3.20+
- A C++23 compiler
  - Windows: Visual Studio 2022 with "Desktop development with C++"
  - Linux: GCC 13+ or Clang 16+
  - macOS: AppleClang that supports C++23 features you use, or Homebrew LLVM
- Git

## Build

`cmake -S . -B build` \
`cmake --build build`

## Run Tests

`ctest --test-dir build`
