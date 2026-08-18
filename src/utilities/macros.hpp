#pragma once

#include <xpu/config.hpp>

#ifdef FP_64
  using real_t = double;
  inline constexpr real_t EPSILON{1e-12};
#else
  using real_t = float;
  inline constexpr real_t EPSILON{1e-5};
#endif

[[nodiscard]] CUDA_CALLABLE
constexpr real_t operator""_r(long double value) {
  return static_cast<real_t>(value);
}
[[nodiscard]] CUDA_CALLABLE 
constexpr real_t operator""_r(unsigned long long value) {
  return static_cast<real_t>(value);
}