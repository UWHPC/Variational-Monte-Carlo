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

/*
TODO: is this really necessary...
Note to self: what I'm thinking is: theres two ideas.
Physics shouldn't be maintained in two places but with the way this function is called,
we end up maintaining the indexing in two places while the kernel just has one multiplication
seems weird and not needed, probably should rethink where the abstraction lives.
*/
CUDA_CALLABLE inline
real_t determinant_ratio_term(
  const real_t new_row_idx, const real_t inv_det_idx
) {
  return new_row_idx * inv_det_idx;
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
