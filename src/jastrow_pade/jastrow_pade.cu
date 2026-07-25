#include "jastrow_pade.cuh"

#include <cmath>
#include <cstddef>

#ifdef VMC_CUDA_BACKEND
namespace {

__global__
void cudaComputeDerivatives(
  std::size_t moved, const std::size_t N,
  real_t old_x, real_t old_y, real_t old_z,
  real_t new_x, real_t new_y, real_t new_z,
  const real_t L, const real_t half_L,
  const real_t a_local, const real_t b_local, const real_t neg2ab,
  const real_t* RESTRICT p_x, const real_t* RESTRICT p_y, const real_t* RESTRICT p_z,
  real_t* RESTRICT grad_x, real_t* RESTRICT grad_y, real_t* RESTRICT grad_z,
  real_t* RESTRICT laplacian
) {
  const auto [i]{vmc::cudaThreadIdx<1>()};
  if (i >= N || i == moved) return; 

  real_t displ_old_x{old_x - p_x[i]};
  real_t displ_old_y{old_y - p_y[i]};
  real_t displ_old_z{old_z - p_z[i]};

  displ_old_x += L * (displ_old_x <= -half_L) + -L * (displ_old_x > half_L);
  displ_old_y += L * (displ_old_y <= -half_L) + -L * (displ_old_y > half_L);
  displ_old_z += L * (displ_old_z <= -half_L) + -L * (displ_old_z > half_L);

  const real_t dist_old{
    vmc::sqrt(
      displ_old_x * displ_old_x +
      displ_old_y * displ_old_y +
      displ_old_z * displ_old_z
    )
  };

  real_t grad_factor_old{0.0_r};
  real_t lap_pair_old{0.0_r};

  if (dist_old >= 1e-12_r) {
    const real_t inv_dist_old{1.0_r / dist_old};
    const real_t inv_denom_old{1.0_r / (1.0_r + b_local * dist_old)};
    const real_t inv_denom_sq_old{inv_denom_old * inv_denom_old};

    const real_t first_deriv_old{a_local * inv_denom_sq_old};
    const real_t second_deriv_old{neg2ab * inv_denom_sq_old * inv_denom_old};

    grad_factor_old = a_local * first_deriv_old * inv_dist_old;
    lap_pair_old = second_deriv_old + 2.0_r * first_deriv_old * inv_dist_old;
  }

  real_t displ_new_x{new_x - p_x[i]};
  real_t displ_new_y{new_y - p_y[i]};
  real_t displ_new_z{new_z - p_z[i]};

  displ_new_x += L * (displ_new_x <= -half_L) + -L * (displ_new_x > half_L);
  displ_new_y += L * (displ_new_y <= -half_L) + -L * (displ_new_y > half_L);
  displ_new_z += L * (displ_new_z <= -half_L) + -L * (displ_new_z > half_L);

  const real_t dist_new{
    vmc::sqrt(
      displ_new_x * displ_new_x +
      displ_new_y * displ_new_y +
      displ_new_z * displ_new_z
    )
  };

  real_t grad_factor_new{0.0_r};
  real_t lap_pair_new{0.0_r};

  if (dist_old >= 1e-12_r) {
    const real_t inv_dist_new{1.0_r / dist_new};
    const real_t inv_denom_new{1.0_r / (1.0_r + b_local * dist_new)};
    const real_t inv_denom_sq_new{inv_denom_new * inv_denom_new};

    const real_t first_deriv_new{a_local * inv_denom_sq_new};
    const real_t second_deriv_new{neg2ab * inv_denom_sq_new * inv_denom_new};

    grad_factor_new = a_local * first_deriv_new * inv_dist_new;
    lap_pair_new = second_deriv_new + 2.0_r * first_deriv_new * inv_dist_new;
  }

  atomicAdd(&grad_x[moved], grad_factor_new * displ_new_x);
  atomicAdd(&grad_y[moved], grad_factor_new * displ_new_y);
  atomicAdd(&grad_z[moved], grad_factor_new * displ_new_z);
  atomicAdd(&laplacian[moved], lap_pair_new);

  grad_x[i] += grad_factor_old * displ_old_x - grad_factor_new * displ_new_x;
  grad_y[i] += grad_factor_old * displ_old_y - grad_factor_new * displ_new_y;
  grad_z[i] += grad_factor_old * displ_old_z - grad_factor_new * displ_new_z;

  laplacian[i] += lap_pair_new - lap_pair_old;
}

}
#endif

void JastrowPade::update_derivatives_for_move(
  const Particles& particles,
  std::size_t moved,
  real_t old_x,
  real_t old_y,
  real_t old_z,
  real_t* RESTRICT grad_x,
  real_t* RESTRICT grad_y,
  real_t* RESTRICT grad_z,
  real_t* RESTRICT laplacian
) const noexcept {
#ifdef VMC_CUDA_BACKEND
  grad_x[moved] = 0.0_r;
  grad_y[moved] = 0.0_r;
  grad_z[moved] = 0.0_r;
  laplacian[moved] = 0.0_r;

  const std::size_t num_particles{particles.size()};
  const real_t L{box_length_};
  const real_t half_L{0.5_r * L};

  const auto p{particles.pos()};

  const real_t a_local{a()};
  const real_t b_local{b()};
  const real_t neg2ab{-2.0_r * a_local * b_local};

  const real_t new_x{p.x_[moved]};
  const real_t new_y{p.y_[moved]};
  const real_t new_z{p.z_[moved]};

  dim3 computeDerivativesThreads(256);
  dim3 computeDerivativesBlocks(vmc::cudaNumBlocks(num_particles, computeDerivativesThreads.x));

  cudaComputeDerivatives<<<computeDerivativesBlocks, computeDerivativesThreads>>>(
    moved, num_particles,
    old_x, old_y, old_z, 
    new_x, new_y, new_z, 
    L, half_L, 
    a_local, b_local, neg2ab,
    p.x_, p.y_, p.z_, 
    grad_x, grad_y, grad_z, 
    laplacian
  );

  #else
  grad_x[moved] = 0.0_r;
  grad_y[moved] = 0.0_r;
  grad_z[moved] = 0.0_r;
  laplacian[moved] = 0.0_r;

  const std::size_t num_particles{particles.size()};
  const real_t L{box_length_};
  const real_t neg_L{-1.0_r * L};

  const real_t half_L{0.5_r * L};
  const real_t neg_half_L{-1.0_r * half_L};

  const auto p{particles.pos().align()};

  ASSUME_ALIGNED(grad_x, SIMD_BYTES);
  ASSUME_ALIGNED(grad_y, SIMD_BYTES);
  ASSUME_ALIGNED(grad_z, SIMD_BYTES);
  ASSUME_ALIGNED(laplacian, SIMD_BYTES);

  const real_t a_local{a()};
  const real_t b_local{b()};
  const real_t neg2ab{-2.0_r * a_local * b_local};

  const real_t new_x{p.x_[moved]};
  const real_t new_y{p.y_[moved]};
  const real_t new_z{p.z_[moved]};

  real_t m_grad_x{}, m_grad_y{}, m_grad_z{}, m_lap{};

  // Not vectorized: loop contains control flow
  #pragma omp simd reduction(+ : m_grad_x, m_grad_y, m_grad_z, m_lap)
  for (std::size_t j = 0; j < num_particles; ++j) {
    const bool is_moved{j == moved};

    real_t displ_old_x{old_x - p.x_[j]};
    real_t displ_old_y{old_y - p.y_[j]};
    real_t displ_old_z{old_z - p.z_[j]};

    displ_old_x += L * (displ_old_x <= neg_half_L) + neg_L * (displ_old_x > half_L);
    displ_old_y += L * (displ_old_y <= neg_half_L) + neg_L * (displ_old_y > half_L);
    displ_old_z += L * (displ_old_z <= neg_half_L) + neg_L * (displ_old_z > half_L);

    const real_t dist_old{
      vmc::sqrt(
        displ_old_x * displ_old_x +
        displ_old_y * displ_old_y +
        displ_old_z * displ_old_z
      )
    };

    const real_t inv_dist_old{(dist_old < 1e-12_r) ? 1.0_r : 1.0_r / dist_old};
    const real_t mask_old{(is_moved || dist_old < 1e-12_r) ? 0.0_r : 1.0_r};

    const real_t inv_denom_old{1.0_r / (1.0_r + b_local * dist_old)};
    const real_t inv_denom_sq_old{inv_denom_old * inv_denom_old};

    const real_t first_deriv_old{a_local * inv_denom_sq_old};
    const real_t second_deriv_old{neg2ab * inv_denom_sq_old * inv_denom_old};

    const real_t grad_factor_old{mask_old * first_deriv_old * inv_dist_old};
    const real_t lap_pair_old{mask_old * (second_deriv_old + 2.0_r * first_deriv_old * inv_dist_old)};

    real_t displ_new_x{new_x - p.x_[j]};
    real_t displ_new_y{new_y - p.y_[j]};
    real_t displ_new_z{new_z - p.z_[j]};

    displ_new_x += L * (displ_new_x <= neg_half_L) + neg_L * (displ_new_x > half_L);
    displ_new_y += L * (displ_new_y <= neg_half_L) + neg_L * (displ_new_y > half_L);
    displ_new_z += L * (displ_new_z <= neg_half_L) + neg_L * (displ_new_z > half_L);

    const real_t dist_new{
      vmc::sqrt(
        displ_new_x * displ_new_x +
        displ_new_y * displ_new_y +
        displ_new_z * displ_new_z
      )
    };

    const real_t inv_dist_new{(dist_new < 1e-12_r) ? 1.0_r : 1.0_r / dist_new};
    const real_t mask_new{(is_moved || dist_new < 1e-12_r) ? 0.0_r : 1.0_r};

    const real_t inv_denom_new{1.0_r / (1.0_r + b_local * dist_new)};
    const real_t inv_denom_sq_new{inv_denom_new * inv_denom_new};

    const real_t first_deriv_new{a_local * inv_denom_sq_new};
    const real_t second_deriv_new{neg2ab * inv_denom_sq_new * inv_denom_new};

    const real_t grad_factor_new{mask_new * first_deriv_new * inv_dist_new};
    const real_t lap_pair_new{mask_new * (second_deriv_new + 2.0_r * first_deriv_new * inv_dist_new)};

    m_grad_x += grad_factor_new * displ_new_x;
    m_grad_y += grad_factor_new * displ_new_y;
    m_grad_z += grad_factor_new * displ_new_z;
    m_lap += lap_pair_new;

    grad_x[j] += grad_factor_old * displ_old_x - grad_factor_new * displ_new_x;
    grad_y[j] += grad_factor_old * displ_old_y - grad_factor_new * displ_new_y;
    grad_z[j] += grad_factor_old * displ_old_z - grad_factor_new * displ_new_z;

    laplacian[j] += lap_pair_new - lap_pair_old;
  }

  grad_x[moved] = m_grad_x;
  grad_y[moved] = m_grad_y;
  grad_z[moved] = m_grad_z;
  laplacian[moved] = m_lap;
#endif
}
