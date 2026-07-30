#include "slater_plane_wave.cuh"
#include "slater_kernels.cuh"

#ifdef VMC_CUDA_BACKEND
#include <cstddef>
namespace {

__global__ void cudaBuildRow(
  std::size_t N, std::size_t trig_S, std::size_t particle,
  const real_t* RESTRICT sin_cache, const real_t* RESTRICT cos_cache,
  const std::size_t* RESTRICT k_idx, const std::uint8_t* RESTRICT orb_t,
  real_t* RESTRICT row
) {
  const auto [i]{vmc::cudaThreadIdx<1>()};
  if (i >= N) { return; }

  build_row_cell(i, particle, trig_S, sin_cache, cos_cache, k_idx[i], orb_t[i], row);
}

__global__ void cudaDeterminantRatio(
  std::size_t N, std::size_t particle, std::size_t S,
  const real_t* RESTRICT new_row, const real_t* RESTRICT inv_det,
  real_t* RESTRICT ratio
) {
  const auto [i]{vmc::cudaThreadIdx<1>()};
  if (i >= N) { return; }

  const std::size_t row_offset{particle * S};

  real_t product{inv_det_dot_term(i, row_offset, new_row, inv_det)};

  atomicAdd(ratio, product);
}

}

#else
#include "particles/particles.cuh"
#include "utilities/aligned_soa.cuh"
#include <cstring>
#endif

real_t* SlaterPlaneWave::build_row(std::size_t particle) noexcept {
#ifdef VMC_CUDA_BACKEND
  dim3 buildRowThreads(256);
  dim3 buildRowBlocks(
    vmc::cudaNumBlocks(this->num_orbitals(), buildRowThreads.x)
  );
  cudaBuildRow<<<buildRowBlocks, buildRowThreads>>>(
    this->num_orbitals(), this->trig_row_stride(), particle,
    this->sin_cache(), this->cos_cache(),
    this->orbital_k_index(), this->orbital_type(),
    this->new_row()
  );
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
  return this->new_row();

#else
  const std::size_t N{num_orbitals()};
  const std::size_t ROW_STRIDE{trig_row_stride()};
  const auto& k_index{orbital_k_index()};
  const auto& orb_type{orbital_type()};

  real_t* RESTRICT row{new_row()};
  real_t* RESTRICT sin_cache{this->sin_cache()};
  real_t* RESTRICT cos_cache{this->cos_cache()};

  ASSUME_ALIGNED(row, SIMD_BYTES);
  ASSUME_ALIGNED(sin_cache, SIMD_BYTES);
  ASSUME_ALIGNED(cos_cache, SIMD_BYTES);

  // Not vectorized: loop-carried data dependency
  #pragma omp simd
  for (std::size_t i = 0; i < N; ++i) {
    build_row_cell(i, particle, ROW_STRIDE, sin_cache, cos_cache, k_index[i], orb_type[i], row);
  }

  return row;
#endif
}

real_t SlaterPlaneWave::determinant_ratio(
  std::size_t particle,
  const real_t* new_row
) const noexcept {
#ifdef VMC_CUDA_BACKEND
  dim3 determinantRatioThreads(256);
  dim3 determinantRatioBlocks(
    vmc::cudaNumBlocks(this->num_orbitals(), determinantRatioThreads.x)
  );

  AlignedSoA<real_t> ratio{1, 1};

  cudaDeterminantRatio<<<determinantRatioBlocks, determinantRatioThreads>>>(
    this->num_orbitals(), particle, this->matrix_row_stride(),
    new_row, this->inv_determinant(),
    ratio[0]
  );

  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
  return *ratio[0];
#else
  const std::size_t N{num_orbitals()};
  const std::size_t S{matrix_row_stride()};
  const real_t* RESTRICT inv_det{inv_determinant()};

  ASSUME_ALIGNED(inv_det, SIMD_BYTES);
  real_t ratio{};
  const std::size_t row_offset{particle * S};

  #pragma omp simd reduction(+ : ratio)
  for (std::size_t j = 0; j < N; ++j) {
    ratio += inv_det_dot_term(j, row_offset, new_row, inv_det);
  }

  return ratio;
#endif
}
