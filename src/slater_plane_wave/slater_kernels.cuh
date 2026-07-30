#pragma once

#include "../utilities/macros.cuh"
#include "../utilities/math.cuh"
#include <cstdint>

CUDA_CALLABLE inline
void build_row_cell(
  std::size_t idx, std::size_t row_offset,
  const real_t* RESTRICT sin_cache, const real_t* RESTRICT cos_cache,
  const std::size_t k_idx, const std::uint8_t orb_t_idx,
  real_t* RESTRICT row
) {
  const real_t type{static_cast<real_t>(orb_t_idx)};

  const std::size_t cache_idx{row_offset + k_idx};

  const real_t cos_term{cos_cache[cache_idx]};
  const real_t sin_term{sin_cache[cache_idx]};

  row[idx] = cos_term + type * (sin_term - cos_term);
}

CUDA_CALLABLE inline
real_t determinant_ratio_term(
  std::size_t idx, std::size_t row_offset,
  const real_t* RESTRICT new_row, const real_t* RESTRICT inv_det
) {
  return new_row[idx] * inv_det[row_offset + idx];
}

CUDA_CALLABLE inline
void trig_cell(
  std::size_t idx, std::size_t row_offset,
  real_t p_x, real_t p_y, real_t p_z,
  const real_t* RESTRICT k_x, const real_t* RESTRICT k_y, const real_t* RESTRICT k_z,
  real_t* RESTRICT sin_cache, real_t* RESTRICT cos_cache
) {
  const real_t dot{
    k_x[idx] * p_x +
    k_y[idx] * p_y +
    k_z[idx] * p_z
  };

  const std::size_t cache_idx{row_offset + idx};

  vmc::sincos(dot, &sin_cache[cache_idx], &cos_cache[cache_idx]);
}

CUDA_CALLABLE inline
real_t log_abs_det_term(
  std::size_t idx, std::size_t mat_S,
  const real_t* RESTRICT lower_upper
) {
  const real_t U_diag{lower_upper[idx * mat_S + idx]};

  return vmc::log(vmc::abs(U_diag));
}

struct [[nodiscard]] SlaterDerivativeTerms {
  real_t grad_x;
  real_t grad_y;
  real_t grad_z;
  real_t laplacian;
};

CUDA_CALLABLE inline
SlaterDerivativeTerms slater_derivative_terms(
  std::size_t idx, std::size_t row_offset, std::size_t mat_offset,
  const real_t* RESTRICT k_x, const real_t* RESTRICT k_y, const real_t* RESTRICT k_z,
  const std::size_t k_idx, const std::uint8_t orb_t_idx,
  const real_t* RESTRICT inv_det,
  const real_t* RESTRICT sin_cache, const real_t* RESTRICT cos_cache
) {
  const real_t kx{k_x[k_idx]};
  const real_t ky{k_y[k_idx]};
  const real_t kz{k_z[k_idx]};

  const real_t k_sq{kx * kx + ky * ky + kz * kz};

  const real_t type{static_cast<real_t>(orb_t_idx)};

  const std::size_t cache_idx{row_offset + k_idx};

  const real_t cos_term{cos_cache[cache_idx]};
  const real_t sin_term{sin_cache[cache_idx]};

  const real_t grad_factor{-sin_term + type * (sin_term + cos_term)};
  const real_t lap_factor{-cos_term + type * (cos_term - sin_term)};

  const real_t weight{inv_det[mat_offset + idx]};

  return {
    weight * kx * grad_factor,
    weight * ky * grad_factor,
    weight * kz * grad_factor,
    weight * k_sq * lap_factor
  };
}

CUDA_CALLABLE inline
void inv_det_scale_cell(
  std::size_t idx, std::size_t row_offset,
  const real_t* RESTRICT inv_d_col, real_t inv_ratio,
  real_t* RESTRICT inv_det
) {
  inv_det[row_offset + idx] = inv_d_col[idx] * inv_ratio;
}

CUDA_CALLABLE inline
void inv_det_update_cell(
  std::size_t idx, std::size_t row_offset,
  const real_t* RESTRICT inv_d_col, real_t factor,
  real_t* RESTRICT inv_det
) {
  inv_det[row_offset + idx] -= inv_d_col[idx] * factor;
}
