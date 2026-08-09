#include <xpu/xpu.hpp>
#include "slater_plane_wave.hpp"
#include "particles/particles.hpp"
#include <xpu/soa.hpp>

#include <cmath>
#include <cstddef>
#include <cstring>

SlaterPlaneWave::SlaterPlaneWave(const Particles& particles, real_t box_lengthL)
  : num_orbitals_{particles.count()}
  , matrix_row_stride_{xpu::handle_pad<real_t>(particles.count())}
  , matrix_size_{matrix_row_stride_ * particles.count()}
  , box_length_{box_lengthL}
  , orbital_k_index_(particles.count())
  , orbital_type_(particles.count())
  , int_vec_{particles.count()}
  , fp_vec_{particles.count()}
  , trig_cache_{0}
  , trig_scratch_{0}
  , matrices_{matrix_row_stride_ * particles.count()}
{
  this->initialize(particles);
#ifdef XPU_CUDA
  this->init_cuda_scratch();
#endif
};

#ifdef XPU_CUDA
namespace {

__global__
void kComputeSK(
  std::size_t N,
  std::size_t S,
  std::size_t particle,
  const real_t* RESTRICT new_row,
  const real_t* RESTRICT inv_det,
  real_t* RESTRICT s_k_array
) {
  const std::size_t m{blockIdx.x * blockDim.x + threadIdx.x};
  const std::size_t k{blockIdx.y * blockDim.y + threadIdx.y};

  if (m >= N || k >= N) return;
  if (k == particle) return;

  const real_t product{new_row[m] * inv_det[k * S + m]};
  atomicAdd(&s_k_array[k], product);
}

__global__
void kUpdateInverse(
  std::size_t N,
  std::size_t S,
  std::size_t particle,
  const real_t* RESTRICT inv_d_col,
  const real_t* RESTRICT s_k_array,
  real_t inv_ratio,
  real_t* RESTRICT inv_det
) {
  const std::size_t j{blockIdx.x * blockDim.x + threadIdx.x};
  const std::size_t k{blockIdx.y * blockDim.y + threadIdx.y};

  if (j >= N || k >= N) return;

  const std::size_t idx{k * S + j};

  if (k == particle) {
    inv_det[idx] = inv_d_col[j] * inv_ratio;
  } else {
    const real_t factor{s_k_array[k] * inv_ratio};
    inv_det[idx] -= inv_d_col[j] * factor;
  }
}

}

#endif

void SlaterPlaneWave::accept_move(
  std::size_t particle,
  const real_t* new_row,
  real_t ratio
) noexcept {
  #ifdef XPU_CUDA
  const std::size_t N{num_orbitals()};
  const std::size_t S{matrix_row_stride()};
  const real_t inv_ratio{1.0_r / ratio};

  if (!std::isfinite(inv_ratio)) return;

  real_t* RESTRICT inv_det{inv_determinant()};
  real_t* RESTRICT det_matrix{determinant()};
  real_t* RESTRICT inv_d_col{this->inv_d_col()};
  real_t* RESTRICT s_k_array{this->solution()};

  const std::size_t p_offset{particle * S};

  xpu::cu_check(cudaMemcpyAsync(inv_d_col, &inv_det[p_offset], N * sizeof(real_t), cudaMemcpyDeviceToDevice));

  xpu::cu_check(cudaMemsetAsync(s_k_array, 0, N * sizeof(real_t)));

  dim3 threads(16, 16);
  dim3 blocks(xpu::block_per_dim(N, threads.x), xpu::block_per_dim(N, threads.y));

  kComputeSK<<<blocks, threads>>>(
    N, S, particle, new_row, inv_det, s_k_array
  );
  xpu::cu_check(cudaGetLastError());

  kUpdateInverse<<<blocks, threads>>>(
    N, S, particle, inv_d_col, s_k_array, inv_ratio, inv_det
  );
  xpu::cu_check(cudaGetLastError());
  xpu::cu_check(cudaMemcpyAsync(&det_matrix[p_offset], new_row, N * sizeof(real_t), cudaMemcpyDeviceToDevice));
  xpu::cu_check(cudaDeviceSynchronize());
  return;
  #else
  const std::size_t N{num_orbitals()};
  const std::size_t S{matrix_row_stride()};
  const real_t inv_ratio{1.0_r / ratio};

  if (!std::isfinite(inv_ratio)) {
    return;
  }

  real_t* RESTRICT inv_det{inv_determinant()};
  real_t* RESTRICT det_matrix{determinant()};
  real_t* RESTRICT inv_d_col{this->inv_d_col()};
  const std::size_t p_offset{particle * S}; 

  // Cache particle row column j for inv_D before changing
  std::memcpy(inv_d_col, &inv_det[p_offset], N * sizeof(real_t));

  // Follows Sherman-Morrison update (branchless)
  for (std::size_t k = 0; k < N; ++k) {
    if (k == particle)
      continue;

    const std::size_t k_offset{k * S};
    real_t s_k{};

    #pragma omp simd reduction(+ : s_k)
    for (std::size_t m = 0; m < N; ++m) {
      s_k += new_row[m] * inv_det[k_offset + m];
    }

    const real_t factor{s_k * inv_ratio};

    #pragma omp simd
    for (std::size_t j = 0; j < N; ++j) {
      inv_det[k_offset + j] -= inv_d_col[j] * factor;
    }
  }

  #pragma omp simd
  for (std::size_t j = 0; j < N; ++j) {
    inv_det[p_offset + j] = inv_d_col[j] * inv_ratio;
  }

  std::memcpy(&det_matrix[p_offset], new_row, N * sizeof(real_t));
  #endif
}
