#pragma once

#include "../utilities/macros.cuh"
#include "../utilities/math.cuh"
#include "../utilities/geometry.cuh"

CUDA_CALLABLE inline
real_t reciprocal_energy_term(
  std::size_t idx,
  const real_t* RESTRICT g_weights,
  const real_t* RESTRICT sum_real, const real_t* RESTRICT sum_imag
) {
  const real_t re{sum_real[idx]};
  const real_t im{sum_imag[idx]};

  return g_weights[idx] * (re * re + im * im);
}

CUDA_CALLABLE inline
real_t real_energy_term(
  std::size_t idx,
  real_t self_x, real_t self_y, real_t self_z,
  real_t L, real_t alpha,
  const real_t* RESTRICT pos_x, const real_t* RESTRICT pos_y, const real_t* RESTRICT pos_z
) {
  const Displacement d{
    minimum_image_displacement(self_x, self_y, self_z, pos_x[idx], pos_y[idx], pos_z[idx], L)
  };

  // Protect against 1.0 / 0.0 generating NaN
  const real_t inv_r{(d.dist < 1e-12_r) ? 1.0_r : 1.0_r / d.dist};

  return vmc::erfc(alpha * d.dist) * inv_r;
}

CUDA_CALLABLE inline
real_t real_energy_delta_term(
  std::size_t idx, std::size_t moved,
  real_t old_x, real_t old_y, real_t old_z,
  real_t new_x, real_t new_y, real_t new_z,
  real_t L, real_t alpha,
  const real_t* RESTRICT pos_x, const real_t* RESTRICT pos_y, const real_t* RESTRICT pos_z
) {
  // Branchless mask to safely skip the moved particle
  const real_t valid_mask{(idx == moved) ? 0.0_r : 1.0_r};

  const real_t erfc_old{
    real_energy_term(idx, old_x, old_y, old_z, L, alpha, pos_x, pos_y, pos_z)
  };
  const real_t erfc_new{
    real_energy_term(idx, new_x, new_y, new_z, L, alpha, pos_x, pos_y, pos_z)
  };

  return valid_mask * (erfc_new - erfc_old);
}

struct [[nodiscard]] StructureFactorTerm {
  real_t cos_term;
  real_t sin_term;
};

CUDA_CALLABLE inline
StructureFactorTerm structure_factor_term(
  std::size_t idx,
  real_t g_x, real_t g_y, real_t g_z,
  const real_t* RESTRICT pos_x, const real_t* RESTRICT pos_y, const real_t* RESTRICT pos_z
) {
  const real_t G_dot_r{
    g_x * pos_x[idx] +
    g_y * pos_y[idx] +
    g_z * pos_z[idx]
  };

  real_t sin_temp{};
  real_t cos_temp{};

  vmc::sincos(G_dot_r, &sin_temp, &cos_temp);

  return {cos_temp, sin_temp};
}

CUDA_CALLABLE inline
StructureFactorTerm structure_factor_delta(
  std::size_t idx,
  const real_t* RESTRICT g_x, const real_t* RESTRICT g_y, const real_t* RESTRICT g_z,
  real_t old_x, real_t old_y, real_t old_z,
  real_t new_x, real_t new_y, real_t new_z
) {
  const real_t old_dot{
    g_x[idx] * old_x +
    g_y[idx] * old_y +
    g_z[idx] * old_z
  };
  const real_t new_dot{
    g_x[idx] * new_x +
    g_y[idx] * new_y +
    g_z[idx] * new_z
  };

  real_t old_sin{}, old_cos{};
  real_t new_sin{}, new_cos{};

  vmc::sincos(old_dot, &old_sin, &old_cos);
  vmc::sincos(new_dot, &new_sin, &new_cos);

  return {new_cos - old_cos, new_sin - old_sin};
}

CUDA_CALLABLE inline
real_t kinetic_energy_term(
  std::size_t idx,
  const real_t* RESTRICT grad_x, const real_t* RESTRICT grad_y, const real_t* RESTRICT grad_z,
  const real_t* RESTRICT lap
) {
  // ||Grad(logPsi)||^2 for particle idx
  const real_t grad_sq{
    grad_x[idx] * grad_x[idx] +
    grad_y[idx] * grad_y[idx] +
    grad_z[idx] * grad_z[idx]
  };

  // Lapl(logPsi) + ||Grad(logPsi)||^2
  return lap[idx] + grad_sq;
}
