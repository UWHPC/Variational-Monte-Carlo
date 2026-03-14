#pragma once

#include <algorithm>
#include <cstdlib>
#include <memory>

// Restrict pointers:
#if defined(__GNUC__) || defined(__clang__)
#define RESTRICT __restrict__
#elif defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

// Sincos support for all compilers
#if defined(__APPLE__)
#define PORTABLE_SINCOS(theta, sin, cos) __sincos((theta), &(sin), &(cos))
#elif defined(__GNUC__) || defined(__clang__)
#define PORTABLE_SINCOS(theta, sin, cos) sincos((theta), &(sin), &(cos))
#else
#define PORTABLE_SINCOS(theta, sin, cos) \
    do { \
        (sin) = std::sin(theta); \
        (cos) = std::cos(theta); \
    } while (0)
#endif