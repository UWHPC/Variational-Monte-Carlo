#include "slater_plane_wave.cuh"

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
  const std::size_t i{blockIdx.x * blockDim.x + threadIdx.x};
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
  const real_t* RESTRICT inv_det, const real_t* RESTRICT s_cache, const real_t* RESTRICT c_cache,
  real_t* RESTRICT grad_x, real_t* RESTRICT grad_y, real_t* RESTRICT grad_z, real_t* RESTRICT lap
) {
  const std::size_t i{blockIdx.x * blockDim.x + threadIdx.x};
  if (i >= size_orb) { return; }

  const std::size_t j{blockIdx.y * blockDim.y + threadIdx.y};
  if (j >= size_orb) { return; }

  const std::size_t trig_offset{i * trig_S};
  const std::size_t mat_offset{i * mat_S};
  const std::size_t inv_offset{mat_offset + j};

  const std::size_t k{k_index[j]};

  const real_t k_mag{
    kx[k] * kx[k] +
    ky[k] * ky[k] +
    kz[k] * kz[k]
  };

  const real_t trig_t{static_cast<real_t>(orb_t[j])};
  const std::size_t trig_idx{trig_offset + k};

  const real_t sin_k{s_cache[trig_idx]};
  const real_t cos_k{c_cache[trig_idx]};

  const real_t weight{inv_det[inv_offset]};
  const real_t grad_factor{
    weight * (
      -sin_k + trig_t * (sin_k + cos_k)
    )
  };
  const real_t lap_factor{
    weight * (
      -cos_k + trig_t * (cos_k - sin_k)
    )
  };

  atomicAdd(&grad_x[i], grad_factor * kx[k]);
  atomicAdd(&grad_y[i], grad_factor * ky[k]);
  atomicAdd(&grad_z[i], grad_factor * kz[k]);
  atomicAdd(&lap[i], lap_factor * k_mag);
}

} // namespace

void SlaterPlaneWave::add_derivatives(
  real_t* RESTRICT grad_x,
  real_t* RESTRICT grad_y,
  real_t* RESTRICT grad_z,
  real_t* RESTRICT laplacian
) const noexcept {
  const auto kv{this->k_vector().align()};
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
}

#else

void SlaterPlaneWave::add_derivatives(
  real_t* RESTRICT grad_x,
  real_t* RESTRICT grad_y,
  real_t* RESTRICT grad_z,
  real_t* RESTRICT laplacian
) const noexcept {
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
      const std::size_t k_idx{k_index[orbital]};

      const real_t k_x_orbital{kv.x_[k_idx]};
      const real_t k_y_orbital{kv.y_[k_idx]};
      const real_t k_z_orbital{kv.z_[k_idx]};

      const real_t k_sq{
        k_x_orbital * k_x_orbital +
        k_y_orbital * k_y_orbital +
        k_z_orbital * k_z_orbital
      };

      const real_t type{static_cast<real_t>(o_type[orbital])};

      const real_t cos_term{cos_cache[offset + k_idx]};
      const real_t sin_term{sin_cache[offset + k_idx]};

      const real_t grad_factor{-sin_term + type * (sin_term + cos_term)};
      const real_t lap_factor{-cos_term + type * (cos_term - sin_term)};

      const real_t weight{inv_det[p_offset + orbital]};

      d_log_det_dx += weight * k_x_orbital * grad_factor;
      d_log_det_dy += weight * k_y_orbital * grad_factor;
      d_log_det_dz += weight * k_z_orbital * grad_factor;
      laplace_det_term += weight * k_sq * lap_factor;
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
}

#endif
