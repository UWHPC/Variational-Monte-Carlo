#include "slater_plane_wave.cuh"
#include "slater_kernels.cuh"

#ifdef VMC_CUDA_BACKEND

namespace {

__global__
void cudaAccumulateDerivatives(
  std::size_t size_orb,
  const real_t* RESTRICT gx,
  const real_t* RESTRICT gy,
  const real_t* RESTRICT gz,
  real_t* RESTRICT g_mag,
  real_t* RESTRICT lap
) {
  const auto [i]{vmc::cudaThreadIdx<1>()};
  if (i >= size_orb) { return; }

  g_mag[i] = (
    gx[i] * gx[i] +
    gy[i] * gy[i] +
    gz[i] * gz[i]
  );

  lap[i] -= g_mag[i];
}

__global__
void cudaAddDerivatives(
  std::size_t size_orb, std::size_t mat_S, std::size_t trig_S,
  const real_t* RESTRICT kx, const real_t* RESTRICT ky, const real_t* RESTRICT kz,
  const std::size_t* RESTRICT k_index, const uint8_t* RESTRICT orb_t,
  const real_t* RESTRICT inv_det, const real_t* RESTRICT sin_cache, const real_t* RESTRICT cos_cache,
  real_t* RESTRICT grad_x, real_t* RESTRICT grad_y, real_t* RESTRICT grad_z, real_t* RESTRICT lap
) {
  const auto [i, j]{vmc::cudaThreadIdx<2>()};
  if (i >= size_orb || j >= size_orb) { return; }

  const SlaterDerivativeTerms terms{
    slater_derivative_terms(
      j, i * trig_S, i * mat_S,
      kx, ky, kz,
      k_index[j], orb_t[j],
      inv_det, sin_cache, cos_cache
    )
  };

  atomicAdd(&grad_x[i], terms.grad_x);
  atomicAdd(&grad_y[i], terms.grad_y);
  atomicAdd(&grad_z[i], terms.grad_z);
  atomicAdd(&lap[i], terms.laplacian);
}

} // namespace

#endif

void SlaterPlaneWave::add_derivatives(
  real_t* RESTRICT grad_x,
  real_t* RESTRICT grad_y,
  real_t* RESTRICT grad_z,
  real_t* RESTRICT laplacian
) const noexcept {
#ifdef VMC_CUDA_BACKEND
  AlignedSoA<real_t> grad_mag{this->num_orbitals(), 1};

  dim3 addDerivativesThreads(16, 16);
  dim3 addDerivativesBlocks(
    vmc::cudaNumBlocks(this->num_orbitals(), addDerivativesThreads.x),
    vmc::cudaNumBlocks(this->num_orbitals(), addDerivativesThreads.y)
  );
  cudaAddDerivatives<<<addDerivativesBlocks, addDerivativesThreads>>>(
    this->num_orbitals(), this->matrix_row_stride(), this->trig_row_stride(),
    this->k_vector().x_, this->k_vector().y_, this->k_vector().z_,
    this->orbital_k_index(), this->orbital_type(),
    this->inv_determinant(), this->sin_cache(), this->cos_cache(),
    grad_x, grad_y, grad_z, laplacian
  );

  dim3 accumulateDerivativesThreads(256);
  dim3 accumulateDerivativesBlocks(
    vmc::cudaNumBlocks(this->num_orbitals(), accumulateDerivativesThreads.x)
  );
  cudaAccumulateDerivatives<<<accumulateDerivativesBlocks, accumulateDerivativesThreads>>>(
    this->num_orbitals(),
    grad_x, grad_y, grad_z, grad_mag[0],
    laplacian
  );

  CUDA_CHECK(cudaDeviceSynchronize());
#else
  const std::size_t N{this->num_orbitals()};
  const std::size_t S{this->matrix_row_stride()};

  const auto kv{this->k_vector().align()};

  const auto& k_index{this->orbital_k_index()};
  const auto& o_type{this->orbital_type()};

  const real_t* RESTRICT inv_det{this->inv_determinant()}; ASSUME_ALIGNED(inv_det, SIMD_BYTES);
  const real_t* RESTRICT cos_cache{this->cos_cache()}; ASSUME_ALIGNED(cos_cache, SIMD_BYTES);
  const real_t* RESTRICT sin_cache{this->sin_cache()}; ASSUME_ALIGNED(sin_cache, SIMD_BYTES);

  ASSUME_ALIGNED(grad_x, SIMD_BYTES);
  ASSUME_ALIGNED(grad_y, SIMD_BYTES);
  ASSUME_ALIGNED(grad_z, SIMD_BYTES);
  ASSUME_ALIGNED(laplacian, SIMD_BYTES);

  for (std::size_t particle = 0; particle < N; ++particle) {
    real_t d_log_det_dx{}, d_log_det_dy{}, d_log_det_dz{};
    real_t laplace_det_term{};

    const std::size_t offset{particle * this->trig_row_stride()};
    const std::size_t p_offset{particle * S};

    // Added reduction clauses so the compiler can safely vectorize the accumulators
    // Not vectorized: loop-carried data dependency
    #pragma omp simd reduction(+ : d_log_det_dx, d_log_det_dy, d_log_det_dz, laplace_det_term)
    for (std::size_t orbital = 0; orbital < N; ++orbital) {
      const SlaterDerivativeTerms terms{
        slater_derivative_terms(
          orbital, offset, p_offset,
          kv.x_, kv.y_, kv.z_,
          k_index[orbital], o_type[orbital],
          inv_det, sin_cache, cos_cache
        )
      };

      d_log_det_dx += terms.grad_x;
      d_log_det_dy += terms.grad_y;
      d_log_det_dz += terms.grad_z;
      laplace_det_term += terms.laplacian;
    }

    const real_t grad_sq{
      d_log_det_dx * d_log_det_dx +
      d_log_det_dy * d_log_det_dy +
      d_log_det_dz * d_log_det_dz
    };

    grad_x[particle] += d_log_det_dx;
    grad_y[particle] += d_log_det_dy;
    grad_z[particle] += d_log_det_dz;

    laplacian[particle] += (laplace_det_term - grad_sq);
  }
#endif
}
