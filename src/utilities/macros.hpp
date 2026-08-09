#pragma once

#include <xpu/config.hpp>

#include <cstdint>
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

#ifdef XPU_CUDA
#include <cusolverDn.h>

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

[[nodiscard]] CUDA_CALLABLE
constexpr real_t operator""_r(long double value) {
  return static_cast<real_t>(value);
}
[[nodiscard]] CUDA_CALLABLE 
constexpr real_t operator""_r(unsigned long long value) {
  return static_cast<real_t>(value);
}