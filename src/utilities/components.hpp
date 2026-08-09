#pragma once

#include <xpu/xpu.hpp>
#include <cstddef>

enum class Axis { X, Y, Z, NUM };
enum class Derivatives { GRAD_X, GRAD_Y, GRAD_Z, LAP, NUM };

template <typename T> CUDA_CALLABLE
inline constexpr std::size_t idx(T idx) {
  return static_cast<std::size_t>(idx);
}

template <typename A, typename B> CUDA_CALLABLE
inline constexpr std::size_t enum_index(A a, B b) {
  return idx(a) + idx(b);
}
