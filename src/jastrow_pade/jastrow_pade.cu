#include "jastrow_pade.cuh"
#include "jastrow_kernels.cuh"

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
  if (i >= N || i == moved) { return; }

  const JastrowDerivativeTerms moved_terms{
    jastrow_move_cell(
      i, moved,
      old_x, old_y, old_z,
      new_x, new_y, new_z,
      L, a_local, b_local, neg2ab,
      p_x, p_y, p_z,
      grad_x, grad_y, grad_z, laplacian
    )
  };

  atomicAdd(&grad_x[moved], moved_terms.grad_x);
  atomicAdd(&grad_y[moved], moved_terms.grad_y);
  atomicAdd(&grad_z[moved], moved_terms.grad_z);
  atomicAdd(&laplacian[moved], moved_terms.laplacian);
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
    const JastrowDerivativeTerms moved_terms{
      jastrow_move_cell(
        j, moved,
        old_x, old_y, old_z,
        new_x, new_y, new_z,
        L, a_local, b_local, neg2ab,
        p.x_, p.y_, p.z_,
        grad_x, grad_y, grad_z, laplacian
      )
    };

    m_grad_x += moved_terms.grad_x;
    m_grad_y += moved_terms.grad_y;
    m_grad_z += moved_terms.grad_z;
    m_lap += moved_terms.laplacian;
  }

  grad_x[moved] = m_grad_x;
  grad_y[moved] = m_grad_y;
  grad_z[moved] = m_grad_z;
  laplacian[moved] = m_lap;
#endif
}
