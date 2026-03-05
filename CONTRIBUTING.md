# Developer workflow

## Prerequisites

- CMake 3.20+
- A C++ compiler (Clang, GCC, or MSVC)
- LLVM tools: clang-format and clang-tidy
- Python 3.9+

## One command to run everything

From the repo root:

python tools/check.py all

## Format code

python tools/check.py format

## Run only formatting check

python tools/check.py format-check

## Run tests

python tools/check.py test

## Run clang-tidy

python tools/check.py tidy-check
