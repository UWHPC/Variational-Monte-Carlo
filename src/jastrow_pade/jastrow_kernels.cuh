#pragma once

#include "../utilities/macros.cuh"
#include "../utilities/math.cuh"
#include "../utilities/geometry.cuh"

CUDA_CALLABLE inline
real_t jastrow_value_term(
  std::size_t idx, std::size_t self,
  real_t self_x, real_t self_y, real_t self_z,
  real_t L, real_t a, real_t b,
  const real_t* RESTRICT p_x, const real_t* RESTRICT p_y, const real_t* RESTRICT p_z
) {
  const real_t mask{(idx == self) ? 0.0_r : 1.0_r};

  const Displacement d{
    minimum_image_displacement(self_x, self_y, self_z, p_x[idx], p_y[idx], p_z[idx], L)
  };

  const real_t inv_denom{1.0_r / (1.0_r + b * d.dist)};

  return mask * a * d.dist * inv_denom;
}

CUDA_CALLABLE inline
real_t jastrow_delta_term(
  std::size_t idx, std::size_t moved,
  real_t old_x, real_t old_y, real_t old_z,
  real_t new_x, real_t new_y, real_t new_z,
  real_t L, real_t a, real_t b,
  const real_t* RESTRICT p_x, const real_t* RESTRICT p_y, const real_t* RESTRICT p_z
) {
  const real_t valid_mask{(idx == moved) ? 0.0_r : 1.0_r};

  const Displacement d_old{
    minimum_image_displacement(old_x, old_y, old_z, p_x[idx], p_y[idx], p_z[idx], L)
  };
  const Displacement d_new{
    minimum_image_displacement(new_x, new_y, new_z, p_x[idx], p_y[idx], p_z[idx], L)
  };

  const real_t denom_old{1.0_r / (1.0_r + b * d_old.dist)};
  const real_t denom_new{1.0_r / (1.0_r + b * d_new.dist)};

  return valid_mask * a * (d_new.dist * denom_new - d_old.dist * denom_old);
}

struct [[nodiscard]] JastrowPairFactors {
  real_t grad_factor;
  real_t lap_pair;
};

CUDA_CALLABLE inline
JastrowPairFactors jastrow_pair_factors(
  real_t dist, real_t mask,
  real_t a, real_t b, real_t neg2ab
) {
  const bool degenerate{dist < 1e-12_r};

  const real_t inv_dist{degenerate ? 1.0_r : 1.0_r / dist};
  const real_t active{degenerate ? 0.0_r : mask};

  const real_t inv_denom{1.0_r / (1.0_r + b * dist)};
  const real_t inv_denom_sq{inv_denom * inv_denom};

  const real_t first_deriv{a * inv_denom_sq};
  const real_t second_deriv{neg2ab * inv_denom_sq * inv_denom};

  return {
    active * first_deriv * inv_dist,
    active * (second_deriv + 2.0_r * first_deriv * inv_dist)
  };
}

struct [[nodiscard]] JastrowDerivativeTerms {
  real_t grad_x;
  real_t grad_y;
  real_t grad_z;
  real_t laplacian;
};

CUDA_CALLABLE inline
JastrowDerivativeTerms jastrow_derivative_terms(
  std::size_t idx, std::size_t self,
  real_t self_x, real_t self_y, real_t self_z,
  real_t L, real_t a, real_t b, real_t neg2ab,
  const real_t* RESTRICT p_x, const real_t* RESTRICT p_y, const real_t* RESTRICT p_z
) {
  const real_t mask{(idx == self) ? 0.0_r : 1.0_r};

  const Displacement d{
    minimum_image_displacement(self_x, self_y, self_z, p_x[idx], p_y[idx], p_z[idx], L)
  };

  const JastrowPairFactors f{jastrow_pair_factors(d.dist, mask, a, b, neg2ab)};

  return {
    f.grad_factor * d.x,
    f.grad_factor * d.y,
    f.grad_factor * d.z,
    f.lap_pair
  };
}

// Applies the in-place update for particle `idx` and returns the contribution
// the moved particle accumulates.
CUDA_CALLABLE inline
JastrowDerivativeTerms jastrow_move_cell(
  std::size_t idx, std::size_t moved,
  real_t old_x, real_t old_y, real_t old_z,
  real_t new_x, real_t new_y, real_t new_z,
  real_t L, real_t a, real_t b, real_t neg2ab,
  const real_t* RESTRICT p_x, const real_t* RESTRICT p_y, const real_t* RESTRICT p_z,
  real_t* RESTRICT grad_x, real_t* RESTRICT grad_y, real_t* RESTRICT grad_z,
  real_t* RESTRICT laplacian
) {
  const real_t mask{(idx == moved) ? 0.0_r : 1.0_r};

  const Displacement d_old{
    minimum_image_displacement(old_x, old_y, old_z, p_x[idx], p_y[idx], p_z[idx], L)
  };
  const Displacement d_new{
    minimum_image_displacement(new_x, new_y, new_z, p_x[idx], p_y[idx], p_z[idx], L)
  };

  const JastrowPairFactors f_old{jastrow_pair_factors(d_old.dist, mask, a, b, neg2ab)};
  const JastrowPairFactors f_new{jastrow_pair_factors(d_new.dist, mask, a, b, neg2ab)};

  const JastrowDerivativeTerms old_terms{
    f_old.grad_factor * d_old.x,
    f_old.grad_factor * d_old.y,
    f_old.grad_factor * d_old.z,
    f_old.lap_pair
  };
  const JastrowDerivativeTerms new_terms{
    f_new.grad_factor * d_new.x,
    f_new.grad_factor * d_new.y,
    f_new.grad_factor * d_new.z,
    f_new.lap_pair
  };

  grad_x[idx] += old_terms.grad_x - new_terms.grad_x;
  grad_y[idx] += old_terms.grad_y - new_terms.grad_y;
  grad_z[idx] += old_terms.grad_z - new_terms.grad_z;

  laplacian[idx] += new_terms.laplacian - old_terms.laplacian;

  return new_terms;
}
