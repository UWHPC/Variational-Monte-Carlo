#include "slater_plane_wave.cuh"
#include "slater_kernels.cuh"
#include "../utilities/matrix.hpp"
#include "particles/particles.cuh"
#include "utilities/aligned_soa.cuh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <vector>

// Real plane-wave basis from complex exponentials
// The complex basis e^{ik·r} = cos(k·r) + i·sin(k·r) gives two
// linearly independent real orbitals per k-vector:
//   φ_cos(r) = cos(k·r)
//   φ_sin(r) = sin(k·r)
//
// Since cos(k·r) = cos(-k·r) and sin(k·r) = -sin(-k·r),
// the ±k pair maps to the same {cos, sin} pair. We deduplicate
// by keeping only the canonical +n representative (first nonzero
// component positive).
//
// k = 0 is special: sin(0) = 0 for all r, so it contributes
// only one orbital (cos(0) = 1).
//
// Orbital count: N = 1 + 2 × (number of nonzero canonical k-vectors)
// Closed shells: N = 1, 7, 19, 27, 33, 57, ...
SlaterPlaneWave::SlaterPlaneWave(const Particles& particles, real_t box_lengthL)
: num_orbitals_{particles.size()}
, matrix_row_stride_{AlignedSoA<real_t>::round_up(particles.size())}
, matrix_size_{matrix_row_stride_ * particles.size()}
, box_length_{box_lengthL}
, orbital_k_index_(particles.size(), NUM_ORB_K)
, orbital_type_(particles.size(), NUM_ORB_TYPE)
, int_vec_{particles.size(), NUM_INT_VECTORS}
, fp_vec_{particles.size(), NUM_DOUBLE_VECTORS}
, trig_cache_{}
, matrices_{matrix_row_stride_ * particles.size(), NUM_MATRIX}
{
  this->initialize(particles);
#ifdef VMC_CUDA_BACKEND
  this->init_cuda_scratch();
#endif
};

#ifdef VMC_CUDA_BACKEND
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

  atomicAdd(&s_k_array[k], inv_det_dot_term(m, k * S, new_row, inv_det));
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
    inv_det_scale_cell(j, k * S, inv_d_col, inv_ratio, inv_det);
  } else {
    inv_det_update_cell(j, k * S, inv_d_col, s_k_array[k] * inv_ratio, inv_det);
  }
}

}

#endif

void SlaterPlaneWave::accept_move(
  std::size_t particle,
  const real_t* new_row,
  real_t ratio
) noexcept {
  #ifdef VMC_CUDA_BACKEND
  const std::size_t N{num_orbitals()};
  const std::size_t S{matrix_row_stride()};
  const real_t inv_ratio{1.0_r / ratio};

  if (!std::isfinite(inv_ratio)) return;

  real_t* RESTRICT inv_det{inv_determinant()};
  real_t* RESTRICT det_matrix{determinant()};
  real_t* RESTRICT inv_d_col{this->inv_d_col()};
  real_t* RESTRICT s_k_array{this->solution()};

  const std::size_t p_offset{particle * S};

  CUDA_CHECK(cudaMemcpyAsync(inv_d_col, &inv_det[p_offset], N * sizeof(real_t), cudaMemcpyDeviceToDevice));

  CUDA_CHECK(cudaMemsetAsync(s_k_array, 0, N * sizeof(real_t)));

  dim3 threads(16, 16);
  dim3 blocks(vmc::cudaNumBlocks(N, threads.x), vmc::cudaNumBlocks(N, threads.y));

  kComputeSK<<<blocks, threads>>>(
    N, S, particle, new_row, inv_det, s_k_array
  );
  CUDA_CHECK(cudaGetLastError());

  kUpdateInverse<<<blocks, threads>>>(
    N, S, particle, inv_d_col, s_k_array, inv_ratio, inv_det
  );
  CUDA_CHECK(cudaGetLastError());

  CUDA_CHECK(cudaMemcpyAsync(&det_matrix[p_offset], new_row, N * sizeof(real_t), cudaMemcpyDeviceToDevice));

  CUDA_CHECK(cudaDeviceSynchronize());

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
  ASSUME_ALIGNED(inv_det, SIMD_BYTES);
  ASSUME_ALIGNED(det_matrix, SIMD_BYTES);
  ASSUME_ALIGNED(inv_d_col, SIMD_BYTES);

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
      s_k += inv_det_dot_term(m, k_offset, new_row, inv_det);
    }

    const real_t factor{s_k * inv_ratio};

    #pragma omp simd
    for (std::size_t j = 0; j < N; ++j) {
      inv_det_update_cell(j, k_offset, inv_d_col, factor, inv_det);
    }
  }

  #pragma omp simd
  for (std::size_t j = 0; j < N; ++j) {
    inv_det_scale_cell(j, p_offset, inv_d_col, inv_ratio, inv_det);
  }

  std::memcpy(&det_matrix[p_offset], new_row, N * sizeof(real_t));
  #endif
}
