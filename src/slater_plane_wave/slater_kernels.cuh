#pragma once

#include "../utilities/macros.cuh"
#include "../utilities/math.cuh"
#include <cstdint>

CUDA_CALLABLE inline
void build_row_cell(
  std::size_t idx, std::size_t particle, std::size_t trig_S,
  const real_t* RESTRICT sin_cache, const real_t* RESTRICT cos_cache,
  const std::size_t k_idx, const std::uint8_t orb_t_idx,
  real_t* RESTRICT row
) {
  const real_t type{static_cast<real_t>(orb_t_idx)};

  const real_t cos_term{cos_cache[particle * trig_S + k_idx]};
  const real_t sin_term{sin_cache[particle * trig_S + k_idx]};

  row[idx] = cos_term + type * (sin_term - cos_term);
}

CUDA_CALLABLE inline
real_t inv_det_dot_term(
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
