#include "jastrow_pade.cuh"

#ifdef VMC_CUDA_BACKEND

namespace {

__global__
void cudaAddDerivatives(
  std::size_t num_particles,
  real_t L, real_t a, real_t b,
  const real_t* RESTRICT p_x, const real_t* RESTRICT p_y, const real_t* RESTRICT p_z,
  real_t* RESTRICT grad_x, real_t* RESTRICT grad_y, real_t* RESTRICT grad_z,
  real_t* RESTRICT laplacian
) {
  const auto [i, j]{vmc::cudaThreadIdx<2>()};
  if (i >= num_particles || j >= num_particles) { return; }
  if (i == j) { return; }

  auto displ_x{p_x[i] - p_x[j]};
  auto displ_y{p_y[i] - p_y[j]};
  auto displ_z{p_z[i] - p_z[j]};

  const auto neg_L{-1.0_r * L};
  const auto half_L{0.5_r * L};
  const auto neg_half_L{-1.0_r * half_L};


  displ_x += L * (displ_x <= neg_half_L) + neg_L * (displ_x > half_L);
  displ_y += L * (displ_y <= neg_half_L) + neg_L * (displ_y > half_L);
  displ_z += L * (displ_z <= neg_half_L) + neg_L * (displ_z > half_L);

  const auto dist_sq{
    displ_x * displ_x +
    displ_y * displ_y +
    displ_z * displ_z
  };

  const auto dist{vmc::sqrt(dist_sq)};

  if (dist < 1e-12_r) { return; }
  const real_t inv_dist{1.0_r / dist};

  const auto denom{1.0_r / (1.0_r + b * dist)};
  const auto denom_sq{denom * denom};
  const auto denom_cb{denom_sq * denom};

  const auto first_deriv{a * denom_sq};
  const auto second_deriv{(-2.0_r * a * b) * denom_cb};

  const auto grad_factor{first_deriv * inv_dist};
  const auto laplacian_pair{second_deriv + 2.0_r * first_deriv * inv_dist};

  atomicAdd(&grad_x[i], grad_factor * displ_x);
  atomicAdd(&grad_y[i], grad_factor * displ_y);
  atomicAdd(&grad_z[i], grad_factor * displ_z);
  atomicAdd(&laplacian[i], laplacian_pair);
}

}
#endif

void JastrowPade::add_derivatives(
  const Particles& particles,
  real_t* RESTRICT grad_x,
  real_t* RESTRICT grad_y,
  real_t* RESTRICT grad_z,
  real_t* RESTRICT laplacian
) const noexcept {
#ifdef VMC_CUDA_BACKEND
  const std::size_t num_particles{particles.size()};

  dim3 addDerivativesThreads(16, 16);
  dim3 addDerivativesBlocks(
    vmc::cudaNumBlocks(num_particles, addDerivativesThreads.x),
    vmc::cudaNumBlocks(num_particles, addDerivativesThreads.y)
  );

  cudaAddDerivatives<<<addDerivativesBlocks, addDerivativesThreads>>>(
    num_particles,
    box_length_, this->a(), this->b(),
    particles.pos().x_, particles.pos().y_, particles.pos().z_,
    grad_x, grad_y, grad_z, laplacian
  );
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

#else
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
  const real_t neg_two_a_b{-2.0_r * a_local * b_local};

  for (std::size_t i = 0; i < num_particles; ++i) {
    real_t d_grad_x{}, d_grad_y{}, d_grad_z{}, d_lap{};

    // Not vectorized: loop contains control flow
    #pragma omp simd reduction(+ : d_grad_x, d_grad_y, d_grad_z, d_lap)
    for (std::size_t j = 0; j < num_particles; ++j) {
      const real_t valid_idx{i == j ? 0.0_r : 1.0_r};

      real_t displ_x{p.x_[i] - p.x_[j]};
      real_t displ_y{p.y_[i] - p.y_[j]};
      real_t displ_z{p.z_[i] - p.z_[j]};

      displ_x += L * (displ_x <= neg_half_L) + neg_L * (displ_x > half_L);
      displ_y += L * (displ_y <= neg_half_L) + neg_L * (displ_y > half_L);
      displ_z += L * (displ_z <= neg_half_L) + neg_L * (displ_z > half_L);

      const real_t dist_sq{
        displ_x * displ_x +
        displ_y * displ_y +
        displ_z * displ_z
      };
      const real_t dist{vmc::sqrt(dist_sq)};

      const bool degenerate{dist < 1e-12_r};
      const real_t inv_dist{degenerate ? 1.0_r : 1.0_r / dist};
      const real_t mask{degenerate ? 0.0_r : valid_idx};

      const real_t denom{1.0_r / (1.0_r + b_local * dist)};
      const real_t denom_sq{denom * denom};
      const real_t denom_cb{denom_sq * denom};

      const real_t first_deriv{a_local * denom_sq};
      const real_t second_deriv{neg_two_a_b * denom_cb};

      const real_t grad_factor{mask * first_deriv * inv_dist};
      const real_t laplacian_pair{mask * (second_deriv + 2.0_r * first_deriv * inv_dist)};

      d_grad_x += grad_factor * displ_x;
      d_grad_y += grad_factor * displ_y;
      d_grad_z += grad_factor * displ_z;

      d_lap += laplacian_pair;
    }

    grad_x[i] += d_grad_x;
    grad_y[i] += d_grad_y;
    grad_z[i] += d_grad_z;
    laplacian[i] += d_lap;
  }
#endif
}