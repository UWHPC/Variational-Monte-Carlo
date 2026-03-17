#pragma once

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <math.h>

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
inline void PORTABLE_SINCOS(double theta, double* s, double* c) {
    __sincos(theta, s, c);
}
#elif defined(__GNUC__) || defined(__clang__)
#include <math.h>
inline void PORTABLE_SINCOS(double theta, double* s, double* c) {
    sincos(theta, s, c);
}
#else
inline void PORTABLE_SINCOS(double theta, double* s, double* c) {
    *s = std::sin(theta);
    *c = std::cos(theta);
}
#endif

// Hint for the compiler that pointers are aligned
#if defined(__GNUC__) || defined(__clang__)
    #define ASSUME_ALIGNED(ptr, align) \
        (ptr) = static_cast<decltype(ptr)>(__builtin_assume_aligned((ptr), (align)))
#elif defined(_MSC_VER)
    #define ASSUME_ALIGNED(ptr, align) \
        __assume((reinterpret_cast<uintptr_t>(ptr) % (align)) == 0)
#else
    #define ASSUME_ALIGNED(ptr, align) ((void)0)
#endif
