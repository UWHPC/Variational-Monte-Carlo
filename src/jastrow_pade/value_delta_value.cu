#include "jastrow_pade.cuh"

#ifdef VMC_CUDA_BACKEND
#include <cstddef>
namespace {

__global__ void cudaDeltaValue(
  std::size_t num_particles, std::size_t moved,
  real_t old_x, real_t old_y, real_t old_z,
  real_t new_x, real_t new_y, real_t new_z,
  real_t L, real_t a, real_t b,
  const real_t* RESTRICT p_x, const real_t* RESTRICT p_y, const real_t* RESTRICT p_z,
  real_t* RESTRICT delta
) {
  const std::size_t j{blockIdx.x * blockDim.x + threadIdx.x};
  if (j >= num_particles) { return; }

  const real_t neg_L{-1.0_r * L};
  const real_t half_L{0.5_r * L};
  const real_t neg_half_L{-1.0_r * half_L};

  const real_t valid_mask{(j == moved) ? 0.0_r : 1.0_r};

  real_t displ_old_x{old_x - p_x[j]};
  real_t displ_old_y{old_y - p_y[j]};
  real_t displ_old_z{old_z - p_z[j]};

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
  const real_t denom_old{1.0_r / (1.0_r + b * dist_old)};

  real_t displ_new_x{new_x - p_x[j]};
  real_t displ_new_y{new_y - p_y[j]};
  real_t displ_new_z{new_z - p_z[j]};

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
  const real_t denom_new{1.0_r / (1.0_r + b * dist_new)};

  const real_t term{valid_mask * a * (dist_new * denom_new - dist_old * denom_old)};

  atomicAdd(delta, term);
}

}

#else
#include <cmath>
#include <cstddef>
#endif

real_t JastrowPade::delta_value(
  const Particles& particles,
  std::size_t moved,
  real_t old_x,
  real_t old_y,
  real_t old_z
) const noexcept {
#ifdef VMC_CUDA_BACKEND
  const std::size_t num_particles{particles.size()};

  AlignedSoA<real_t> delta{1, 1};

  dim3 deltaValueThreads(256);
  dim3 deltaValueBlocks(
    vmc::cudaNumBlocks(num_particles, deltaValueThreads.x)
  );

  cudaDeltaValue<<<deltaValueBlocks, deltaValueThreads>>>(
    num_particles, moved,
    old_x, old_y, old_z,
    particles.pos().x_[moved], particles.pos().y_[moved], particles.pos().z_[moved],
    box_length_, a(), b(),
    particles.pos().x_, particles.pos().y_, particles.pos().z_,
    delta[0]
  );

  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
  return *delta[0];
#else
  const std::size_t num_particles{particles.size()};
  const real_t L{box_length_};
  const real_t neg_L{-1.0_r * L};
  const real_t half_L{0.5_r * L};
  const real_t neg_half_L{-1.0_r * half_L};

  const auto p{particles.pos().align()};

  const real_t a_local{a()};
  const real_t b_local{b()};

  const real_t new_x{p.x_[moved]};
  const real_t new_y{p.y_[moved]};
  const real_t new_z{p.z_[moved]};

  real_t delta{};

  // Not vectorized: loop contains control flow
  #pragma omp simd reduction(+ : delta)
  for (std::size_t j = 0; j < num_particles; ++j) {
    const real_t valid_mask{(j == moved) ? 0.0_r : 1.0_r};

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
    const real_t denom_old{1.0_r / (1.0_r + b_local * dist_old)};

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
    const real_t denom_new{1.0_r / (1.0_r + b_local * dist_new)};

    delta += valid_mask * a_local * (dist_new * denom_new - dist_old * denom_old);
  }

  return delta;
#endif
}
