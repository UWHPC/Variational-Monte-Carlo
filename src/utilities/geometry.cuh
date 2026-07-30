#pragma once

#include "macros.cuh"
#include "math.cuh"

CUDA_CALLABLE inline
real_t minimum_image(real_t displ, real_t L) {
  const real_t neg_L{-1.0_r * L};
  const real_t half_L{0.5_r * L};
  const real_t neg_half_L{-1.0_r * half_L};

  return displ + L * (displ <= neg_half_L) + neg_L * (displ > half_L);
}

struct [[nodiscard]] Displacement {
  real_t x;
  real_t y;
  real_t z;
  real_t dist;
};

CUDA_CALLABLE inline
Displacement minimum_image_displacement(
  real_t from_x, real_t from_y, real_t from_z,
  real_t to_x, real_t to_y, real_t to_z,
  real_t L
) {
  const real_t displ_x{minimum_image(from_x - to_x, L)};
  const real_t displ_y{minimum_image(from_y - to_y, L)};
  const real_t displ_z{minimum_image(from_z - to_z, L)};

  const real_t dist_sq{
    displ_x * displ_x +
    displ_y * displ_y +
    displ_z * displ_z
  };

  return {displ_x, displ_y, displ_z, vmc::sqrt(dist_sq)};
}
