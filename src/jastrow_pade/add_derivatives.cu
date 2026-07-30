#include "jastrow_pade.cuh"
#include "jastrow_kernels.cuh"

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

  const JastrowDerivativeTerms terms{
    jastrow_derivative_terms(
      j, i,
      p_x[i], p_y[i], p_z[i],
      L, a, b, -2.0_r * a * b,
      p_x, p_y, p_z
    )
  };

  atomicAdd(&grad_x[i], terms.grad_x);
  atomicAdd(&grad_y[i], terms.grad_y);
  atomicAdd(&grad_z[i], terms.grad_z);
  atomicAdd(&laplacian[i], terms.laplacian);
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

  const auto p{particles.pos().align()};

  ASSUME_ALIGNED(grad_x, SIMD_BYTES);
  ASSUME_ALIGNED(grad_y, SIMD_BYTES);
  ASSUME_ALIGNED(grad_z, SIMD_BYTES);
  ASSUME_ALIGNED(laplacian, SIMD_BYTES);

  const real_t a_local{a()};
  const real_t b_local{b()};
  const real_t neg_two_a_b{-2.0_r * a_local * b_local};

  for (std::size_t i = 0; i < num_particles; ++i) {
    const real_t self_x{p.x_[i]};
    const real_t self_y{p.y_[i]};
    const real_t self_z{p.z_[i]};

    real_t d_grad_x{}, d_grad_y{}, d_grad_z{}, d_lap{};

    // Not vectorized: loop contains control flow
    #pragma omp simd reduction(+ : d_grad_x, d_grad_y, d_grad_z, d_lap)
    for (std::size_t j = 0; j < num_particles; ++j) {
      const JastrowDerivativeTerms terms{
        jastrow_derivative_terms(
          j, i,
          self_x, self_y, self_z,
          L, a_local, b_local, neg_two_a_b,
          p.x_, p.y_, p.z_
        )
      };

      d_grad_x += terms.grad_x;
      d_grad_y += terms.grad_y;
      d_grad_z += terms.grad_z;
      d_lap += terms.laplacian;
    }

    grad_x[i] += d_grad_x;
    grad_y[i] += d_grad_y;
    grad_z[i] += d_grad_z;
    laplacian[i] += d_lap;
  }
#endif
}