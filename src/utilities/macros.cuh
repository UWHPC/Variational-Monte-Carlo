#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <math.h>
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

#if defined(__CUDACC__)
#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#define CUDA_CHECK(call) \
  do { \
    cudaError_t cuda_check_result_{(call)}; \
    if (cuda_check_result_ != cudaSuccess) { \
      std::fprintf( \
        stderr, \
        "CUDA error at %s:%d: %s\n", \
        __FILE__, \
        __LINE__, \
        cudaGetErrorString(cuda_check_result_) \
      ); \
      std::abort(); \
    } \
  } while (0)

#define CUSOLVER_CHECK(call) \
  do { \
    cusolverStatus_t cusolver_check_result_{(call)}; \
    if (cusolver_check_result_ != CUSOLVER_STATUS_SUCCESS) { \
      std::fprintf( \
        stderr, \
        "cuSOLVER error at %s:%d: %d\n", \
        __FILE__, \
        __LINE__, \
        static_cast<int>(cusolver_check_result_) \
      ); \
      std::abort(); \
    } \
  } while (0)
#endif

// Cuda Callable
#if defined(__CUDACC__)
  #define CUDA_CALLABLE __host__ __device__
#else
  #define CUDA_CALLABLE
#endif

// Restrict pointers:
#if defined(__GNUC__) || defined(__clang__)
  #define RESTRICT __restrict__
#elif defined(_MSC_VER)
  #define RESTRICT __restrict
#else
  #define RESTRICT
#endif

// SIMD Bytes
#if defined(__AVX512F__)
  constexpr std::size_t SIMD_BYTES{64};
#elif defined(__AVX2__) || defined(__AVX__)
  constexpr std::size_t SIMD_BYTES{32};
#elif defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64) || defined(__ARM_NEON) || defined(__aarch64__)
  constexpr std::size_t SIMD_BYTES{16};
#else
  constexpr std::size_t SIMD_BYTES{alignof(std::max_align_t)};
#endif

// Assume aligned
#if defined(__GNUC__) || defined(__clang__)
#define ASSUME_ALIGNED(ptr, align) \
  (ptr) = static_cast<decltype(ptr)>(__builtin_assume_aligned((ptr), (align)))
#elif defined(_MSC_VER)
  #define ASSUME_ALIGNED(ptr, align) __assume((reinterpret_cast<uintptr_t>(ptr) % (align)) == 0)
#else
  #define ASSUME_ALIGNED(ptr, align) (static_cast<void>(0))
#endif

[[nodiscard]] CUDA_CALLABLE constexpr real_t operator""_r(long double value) {
  return static_cast<real_t>(value);
}
[[nodiscard]] CUDA_CALLABLE constexpr real_t operator""_r(unsigned long long value) {
  return static_cast<real_t>(value);
}