#include "jastrow_pade.cuh"
#include "jastrow_kernels.cuh"

#ifdef VMC_CUDA_BACKEND
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

  atomicAdd(
    jastrow_pade,
    jastrow_value_term(j, i, p_x[i], p_y[i], p_z[i], L, a, b, p_x, p_y, p_z)
  );
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

  atomicAdd(
    delta,
    jastrow_delta_term(
      j, moved,
      old_x, old_y, old_z,
      new_x, new_y, new_z,
      L, a, b,
      p_x, p_y, p_z
    )
  );
}

}

#else
#include <cmath>
#include <cstddef>
#endif

real_t JastrowPade::value(const Particles& particles) const noexcept {
#ifdef VMC_CUDA_BACKEND
  const std::size_t num_particles{particles.size()};

  AlignedSoA<real_t> jastrow_pade{1, 1};

  dim3 valueThreads(16, 16);
  dim3 valueBlocks(
    vmc::cudaNumBlocks(num_particles, valueThreads.x),
    vmc::cudaNumBlocks(num_particles, valueThreads.y)
  );

  cudaValue<<<valueBlocks, valueThreads>>>(
    num_particles,
    box_length_, this->a(), this->b(),
    particles.pos().x_, particles.pos().y_, particles.pos().z_,
    jastrow_pade[0]
  );

  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
  return 0.5_r * *jastrow_pade[0];
#else
  const std::size_t num_particles{particles.size()};
  const real_t L{box_length_};

  const auto p{particles.pos().align()};

  const real_t a_local{a()};
  const real_t b_local{b()};

  real_t jastrow_pade{};
  for (std::size_t i = 0; i < num_particles; ++i) {
    const real_t self_x{p.x_[i]};
    const real_t self_y{p.y_[i]};
    const real_t self_z{p.z_[i]};

    real_t local_jastrow{};

    // Not vectorized: loop contains control flow
    #pragma omp simd reduction(+ : local_jastrow)
    for (std::size_t j = 0; j < num_particles; ++j) {
      local_jastrow += jastrow_value_term(
        j, i,
        self_x, self_y, self_z,
        L, a_local, b_local,
        p.x_, p.y_, p.z_
      );
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
    delta += jastrow_delta_term(
      j, moved,
      old_x, old_y, old_z,
      new_x, new_y, new_z,
      L, a_local, b_local,
      p.x_, p.y_, p.z_
    );
  }

  return delta;
#endif
}
