#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <concepts>
#include <memory>
#include <type_traits>

template <typename T>
concept arithmetic = std::is_arithmetic_v<T>;

#ifdef FP_64
using real_t = double;
#else
using real_t = float;
#endif

using complex_t = std::complex<real_t>;

[[nodiscard]] constexpr real_t operator""_r(long double value) {
  return static_cast<real_t>(value);
}
[[nodiscard]] constexpr real_t operator""_r(unsigned long long value) {
  return static_cast<real_t>(value);
}

// Restrict pointers:
#if defined(__GNUC__) || defined(__clang__)
#define RESTRICT __restrict__
#elif defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#if defined(__AVX512F__)
constexpr std::size_t SIMD_BYTES{64};
#elif defined(__AVX2__) || defined(__AVX__)
constexpr std::size_t SIMD_BYTES{32};
#elif defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64) || defined(__ARM_NEON) ||           \
    defined(__aarch64__)
constexpr std::size_t SIMD_BYTES{16};
#else
constexpr std::size_t SIMD_BYTES{alignof(std::max_align_t)};
#endif

// Sincos support for all compilers
#if defined(__APPLE__)
inline void PORTABLE_SINCOS(double theta, double* s, double* c) { __sincos(theta, s, c); }
inline void PORTABLE_SINCOS(float theta, float* s, float* c) { __sincosf(theta, s, c); }
#elif defined(__GNUC__) || defined(__clang__)
#include <math.h>
inline void PORTABLE_SINCOS(double theta, double* s, double* c) { sincos(theta, s, c); }
inline void PORTABLE_SINCOS(float theta, float* s, float* c) { sincosf(theta, s, c); }
#else
inline void PORTABLE_SINCOS(double theta, double* s, double* c) {
  *s = std::sin(theta);
  *c = std::cos(theta);
}
inline void PORTABLE_SINCOS(float theta, float* s, float* c) {
  *s = std::sin(theta);
  *c = std::cos(theta);
}
#endif

// Hint for the compiler that pointers are aligned
#if defined(__GNUC__) || defined(__clang__)
#define ASSUME_ALIGNED(ptr, align)                                                                 \
  (ptr) = static_cast<decltype(ptr)>(__builtin_assume_aligned((ptr), (align)))
#elif defined(_MSC_VER)
#define ASSUME_ALIGNED(ptr, align) __assume((reinterpret_cast<uintptr_t>(ptr) % (align)) == 0)
#else
#define ASSUME_ALIGNED(ptr, align) (static_cast<void>(0))
#endif
