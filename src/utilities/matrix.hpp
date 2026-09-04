#pragma once

#include "macros.hpp"
#include <xpu/math.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

inline int lower_upper_decomp(
  fp_t* lower_upper,
  int* pivot,
  std::size_t N,
  std::size_t stride
) {
  int swap_count{};

  for (std::size_t row = 0; row < N; ++row) {
    pivot[row] = scast<int>(row);
  }

  for (std::size_t col = 0; col < N; ++col) {
    std::size_t pivot_row{col};
    fp_t max_abs{xpu::abs(lower_upper[col * stride + col])};

    for (std::size_t row = col + 1; row < N; ++row) {
      const fp_t value{xpu::abs(lower_upper[row * stride + col])};
      if (value > max_abs) {
        max_abs = value;
        pivot_row = row;
      }
    }

    constexpr fp_t PIVOT_TOLERANCE{1e-12_fp};
    if (max_abs < PIVOT_TOLERANCE) {
      lower_upper[col * stride + col] = 0.0_fp;
      continue;
    }

    if (pivot_row != col) {
      for (std::size_t col2 = 0; col2 < N; ++col2) {
        std::swap(lower_upper[col * stride + col2], lower_upper[pivot_row * stride + col2]);
      }
      std::swap(pivot[col], pivot[pivot_row]);
      ++swap_count;
    }

    const fp_t pivot_value{lower_upper[col * stride + col]};
    for (std::size_t row = col + 1; row < N; ++row) {
      lower_upper[row * stride + col] /= pivot_value;
      const fp_t multiplier{lower_upper[row * stride + col]};
      for (std::size_t col2 = col + 1; col2 < N; ++col2) {
        lower_upper[row * stride + col2] -= multiplier * lower_upper[col * stride + col2];
      }
    }
  }

  return swap_count;
}

inline void solve_lower_upper(
  const fp_t* lower_upper,
  const int* pivot,
  const fp_t* b,
  fp_t* x,
  std::size_t N,
  std::size_t stride
) {
  for (std::size_t row = 0; row < N; ++row) {
    const std::size_t permuted_row{scast<std::size_t>(pivot[row])};
    x[row] = b[permuted_row];
  }

  for (std::size_t row = 0; row < N; ++row) {
    fp_t sum{x[row]};
    for (std::size_t col = 0; col < row; ++col) {
      sum -= lower_upper[row * stride + col] * x[col];
    }
    x[row] = sum;
  }

  for (std::size_t rev = 0; rev < N; ++rev) {
    const std::size_t row{N - 1 - rev};
    fp_t sum{x[row]};
    for (std::size_t col = row + 1; col < N; ++col) {
      sum -= lower_upper[row * stride + col] * x[col];
    }

    const fp_t diag{lower_upper[row * stride + row]};
    x[row] = (xpu::abs(diag) > std::numeric_limits<fp_t>::min()) ? (sum / diag) : 0.0_fp;
  }
}

CUDA_CALLABLE [[nodiscard]]
inline bool is_canonical(int n_x, int n_y, int n_z) {
  if (n_x > 0) {
    return true;
  }
  if (n_x < 0) {
    return false;
  }
  if (n_y > 0) {
    return true;
  }
  if (n_y < 0) {
    return false;
  }
  if (n_z > 0) {
    return true;
  }
  if (n_z < 0) {
    return false;
  }

  return true;
}