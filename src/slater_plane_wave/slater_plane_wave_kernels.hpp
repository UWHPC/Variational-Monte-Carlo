#pragma once

#include <xpu/launch.hpp>
#include <xpu/math.hpp>
#include "../utilities/components.hpp"
#include "../utilities/macros.hpp"

namespace stencil {
namespace slater {

CUDA_CALLABLE
inline void update_trig_cache(
  std::size_t i, std::size_t offset, std::size_t particle,
  xpu::soa_view<real_t, idx(Axis::NUM)> particle_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> k_vector,
  real_t* RESTRICT sin_cache, real_t* RESTRICT cos_cache
) {
  const auto dot{
    k_vector[idx(Axis::X)][i] * particle_pos[idx(Axis::X)][particle] +
    k_vector[idx(Axis::Y)][i] * particle_pos[idx(Axis::Y)][particle] +
    k_vector[idx(Axis::Z)][i] * particle_pos[idx(Axis::Z)][particle]
  };

  const auto idx{offset + i};
  xpu::sincos(dot, &sin_cache[idx], &cos_cache[idx]);
}

CUDA_CALLABLE
inline void k_update_inverse(
  std::size_t i, std::size_t j,
  std::size_t particle,
  std::size_t row_stride, real_t inv_ratio,
  const real_t* RESTRICT inv_d_col,
  const real_t* RESTRICT solution_arr,
  real_t* RESTRICT inv_det
) {
  const auto idx{j * row_stride + i};
  const auto factor{inv_d_col[i] * inv_ratio};

  if (i == particle) {
    inv_det[idx] = factor;
  } else {
    inv_det[idx] -= factor * solution_arr[j];
  }
}

} // namespace stencil::slater
} // namespace stencil

namespace kernel {
namespace slater {

namespace {

#if defined(XPU_CUDA)
__global__
void cudaUpdateTrigRow(
  std::size_t offset, std::size_t particle,
  xpu::soa_view<real_t, idx(Axis::NUM)> particle_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> k_vector,
  real_t* RESTRICT sin_cache, real_t* RESTRICT cos_cache
) {
  const auto [i]{xpu::global_index<1>()};
  if (i >= k_vector.count()) { return; }

  stencil::slater::update_trig_cache(
    i, offset, particle,
    particle_pos, k_vector,
    sin_cache, cos_cache
  );
}
#endif

} // namespace

inline void update_trig_cache(
  std::size_t offset, std::size_t particle,
  xpu::soa_view<real_t, idx(Axis::NUM)> particle_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> k_vector,
  real_t* RESTRICT sin_cache, real_t* RESTRICT cos_cache
) {
  const auto num_k{k_vector.count()};
#if defined(XPU_CUDA)
  dim3 updateTrigRowThreads{256u};
  dim3 updateTrigRowBlocks{
    xpu::block_per_dim(num_k, updateTrigRowThreads.x)
  };
  cudaUpdateTrigRow<<<
    updateTrigRowBlocks, updateTrigRowThreads
  >>>(
    offset, particle,
    particle_pos, k_vector,
    sin_cache, cos_cache
  );
  xpu::cu_check(cudaGetLastError());
  xpu::cu_check(cudaDeviceSynchronize());
#else
  #pragma omp simd
  for (auto i = 0uz; i < num_k; ++i) {
    stencil::slater::update_trig_cache(
      i, offset, particle,
      particle_pos, k_vector,
      sin_cache, cos_cache
    );
  }
#endif
}

inline void k_update_inverse(
  std::size_t num_orbitals, std::size_t particle,
  std::size_t row_stride, real_t inv_ratio,
  const real_t* RESTRICT inv_d_col,
  const real_t* RESTRICT solution_arr,
  real_t* RESTRICT inv_det
) {
#if defined(XPU_CUDA)
  dim3 kUpdateInverseThreads{16u, 16u};
  dim3 kUpdateInverseBlocks(
    xpu::block_per_dim(num_orbitals, kUpdateInverseThreads.x),
    xpu::block_per_dim(num_orbitals, kUpdateInverseThreads.y)
  );
  kUpdateInverse<<<
    kUpdateInverseBlocks, kUpdateInverseThreads
  >>>(
    N, S, particle, inv_d_col, s_k_array, inv_ratio, inv_det
  );
#else
#endif
  for (auto i = 0uz; i < num_orbitals; ++i) {
    for(auto j = 0uz; j < num_orbitals; ++j) {
      stencil::slater::k_update_inverse(
        i, j,
        particle, row_stride, inv_ratio,
        inv_d_col, solution_arr, inv_det
      );
    }
  }
}

} // namespace kernel::slater
} // namespace kernel
