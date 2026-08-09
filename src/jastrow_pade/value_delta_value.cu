#include <xpu/xpu.hpp>
#include "jastrow_pade.hpp"

#ifdef XPU_CUDA
#include <cstddef>
namespace {

__global__
void cudaValue(
  std::size_t num_particles,
  real_t L, real_t a, real_t b,
  const real_t* RESTRICT p_x, const real_t* RESTRICT p_y, const real_t* RESTRICT p_z,
  real_t* RESTRICT jastrow_pade
) {
  const std::size_t i{blockIdx.x * blockDim.x + threadIdx.x};
  if (i >= num_particles) { return; }

  const std::size_t j{blockIdx.y * blockDim.y + threadIdx.y};
  if (j >= num_particles) { return; }

  const real_t mask{i == j ? 0.0_r : 1.0_r};

  const real_t neg_L{-1.0_r * L};
  const real_t half_L{0.5_r * L};
  const real_t neg_half_L{-1.0_r * half_L};

  real_t displ_x{p_x[i] - p_x[j]};
  real_t displ_y{p_y[i] - p_y[j]};
  real_t displ_z{p_z[i] - p_z[j]};

  displ_x += L * (displ_x <= neg_half_L) + neg_L * (displ_x > half_L);
  displ_y += L * (displ_y <= neg_half_L) + neg_L * (displ_y > half_L);
  displ_z += L * (displ_z <= neg_half_L) + neg_L * (displ_z > half_L);

  const real_t dist_sq{
    displ_x * displ_x +
    displ_y * displ_y +
    displ_z * displ_z
  };
  const real_t dist{xpu::sqrt(dist_sq)};

  const real_t denom{1.0_r + b * dist};
  const real_t inv_denom{1.0_r / denom};

  atomicAdd(jastrow_pade, mask * a * dist * inv_denom);
}

__global__
void cudaDeltaValue(
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
    xpu::sqrt(
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
    xpu::sqrt(
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

real_t JastrowPade::value(const Particles& particles) const noexcept {
#ifdef XPU_CUDA
  const std::size_t num_particles{particles.size()};

  xpu::buffer<real_t> jastrow_pade{1};

  dim3 valueThreads(16, 16);
  dim3 valueBlocks(
    xpu::block_per_dim(num_particles, valueThreads.x),
    xpu::block_per_dim(num_particles, valueThreads.y)
  );

  cudaValue<<<valueBlocks, valueThreads>>>(
    num_particles,
    box_length_, this->a(), this->b(),
    particles.pos().x_, particles.pos().y_, particles.pos().z_,
    jastrow_pade.data()
  );

  xpu::cuda_check(cudaGetLastError());
  xpu::cuda_check(cudaDeviceSynchronize());
  real_t value_host{};
  xpu::cuda_check(cudaMemcpy(&value_host, jastrow_pade.data(), sizeof(real_t), cudaMemcpyDeviceToHost));
  return 0.5_r * value_host;
#else
  const std::size_t num_particles{particles.size()};
  const real_t L{box_length_};
  const real_t neg_L{-1.0_r * L};
  const real_t half_L{0.5_r * L};
  const real_t neg_half_L{-1.0_r * half_L};

  const auto p{particles.pos()};

  const real_t a_local{a()};
  const real_t b_local{b()};

  real_t jastrow_pade{};
  for (std::size_t i = 0; i < num_particles; ++i) {
    real_t local_jastrow{};

    // Not vectorized: loop contains control flow
    #pragma omp simd reduction(+ : local_jastrow)
    for (std::size_t j = 0; j < num_particles; ++j) {
      const real_t mask{i == j ? 0.0_r : 1.0_r};

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
      const real_t dist{xpu::sqrt(dist_sq)};

      const real_t denom{1.0_r + b_local * dist};
      const real_t inv_denom{1.0_r / denom};

      local_jastrow += mask * a_local * dist * inv_denom;
    }
    jastrow_pade += local_jastrow;
  }
  return 0.5_r * jastrow_pade;
#endif
}

real_t JastrowPade::delta_value(
  const Particles& particles,
  std::size_t moved,
  real_t old_x,
  real_t old_y,
  real_t old_z
) const noexcept {
#ifdef XPU_CUDA
  const std::size_t num_particles{particles.size()};

  xpu::buffer<real_t> delta{1};

  dim3 deltaValueThreads(256);
  dim3 deltaValueBlocks(
    xpu::block_per_dim(num_particles, deltaValueThreads.x)
  );

  const auto moved_pos{particles.pos()};
  real_t new_x{}, new_y{}, new_z{};
  xpu::cuda_check(cudaMemcpy(&new_x, moved_pos.x_ + moved, sizeof(real_t), cudaMemcpyDeviceToHost));
  xpu::cuda_check(cudaMemcpy(&new_y, moved_pos.y_ + moved, sizeof(real_t), cudaMemcpyDeviceToHost));
  xpu::cuda_check(cudaMemcpy(&new_z, moved_pos.z_ + moved, sizeof(real_t), cudaMemcpyDeviceToHost));

  cudaDeltaValue<<<deltaValueBlocks, deltaValueThreads>>>(
    num_particles, moved,
    old_x, old_y, old_z,
    new_x, new_y, new_z,
    box_length_, a(), b(),
    particles.pos().x_, particles.pos().y_, particles.pos().z_,
    delta.data()
  );

  xpu::cuda_check(cudaGetLastError());
  xpu::cuda_check(cudaDeviceSynchronize());
  real_t delta_host{};
  xpu::cuda_check(cudaMemcpy(&delta_host, delta.data(), sizeof(real_t), cudaMemcpyDeviceToHost));
  return delta_host;
#else
  const std::size_t num_particles{particles.size()};
  const real_t L{box_length_};
  const real_t neg_L{-1.0_r * L};
  const real_t half_L{0.5_r * L};
  const real_t neg_half_L{-1.0_r * half_L};

  const auto p{particles.pos()};

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
      xpu::sqrt(
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
      xpu::sqrt(
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
