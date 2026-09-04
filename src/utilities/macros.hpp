#pragma once

#include <xpu/config.hpp>
#include <concepts>

#if defined(FP_64)

using fp_t = double;

#else

using fp_t = float;

#endif

template <typename C, typename T>
[[nodiscard]] CUDA_CALLABLE
constexpr auto scast(const T& arg) noexcept -> C {
  return static_cast<C>(arg);
}

#if defined(FP_64)

inline constexpr auto epsilon{scast<fp_t>(1e-12)};

#else

inline constexpr auto epsilon{scast<fp_t>(1e-6)};

#endif

[[nodiscard]] CUDA_CALLABLE
constexpr auto operator""_fp(long double value) noexcept -> fp_t {
  return scast<fp_t>(value);
}

[[nodiscard]] CUDA_CALLABLE
constexpr auto operator""_fp(unsigned long long value) noexcept -> fp_t {
  return scast<fp_t>(value);
}
