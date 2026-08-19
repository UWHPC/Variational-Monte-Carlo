#pragma once

#include <xpu/config.hpp>
#include <concepts>

#ifdef FP_64
  using fp_t = double;
  inline constexpr fp_t EPSILON{1e-12};
#else
  using fp_t = float;
  inline constexpr fp_t EPSILON{1e-5};
#endif

[[nodiscard]] CUDA_CALLABLE
constexpr fp_t operator""_fp(long double value) {
  return static_cast<fp_t>(value);
}
[[nodiscard]] CUDA_CALLABLE 
constexpr fp_t operator""_fp(unsigned long long value) {
  return static_cast<fp_t>(value);
}

template <typename C, typename T> [[nodiscard]] CUDA_CALLABLE
inline constexpr C scast(T arg) {
  return static_cast<C>(arg);
}
