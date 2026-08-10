#pragma once

#include <xpu/xpu.hpp>
#include "../utilities/components.hpp"
#include "../utilities/macros.hpp"

namespace stencil {
namespace wavefunction {

CUDA_CALLABLE
inline void derivative_sum(
  std::size_t i,
  xpu::soa_view<real_t, idx(Derivatives::NUM)> log_derivatives,
  xpu::soa_view<real_t, idx(Derivatives::NUM)> jastrow_derivatives
) {
  for (auto d{idx(Derivatives::GRAD_X)}; d < idx(Derivatives::NUM); ++d) {
    log_derivatives[d][i] += jastrow_derivatives[d][i];
  }
}

} // namespace stencil::wavefunction
} // namespace stencil

namespace kernel {
namespace wavefunction {

namespace {

#if defined(XPU_CUDA)
__global__
void cudaDerivativeSums(
  xpu::soa_view<real_t, idx(Derivatives::NUM)> log_derivatives,
  xpu::soa_view<real_t, idx(Derivatives::NUM)> jastrow_derivatives
) {
  const auto [i]{xpu::global_index<1>()};
  if (i >= log_derivatives.stride()) { return; }

  stencil::wavefunction::derivative_sum(i,log_derivatives, jastrow_derivatives);
}
#endif

} // namespace

inline void derivative_sums(
  xpu::soa_view<real_t, idx(Derivatives::NUM)> log_derivatives,
  xpu::soa_view<real_t, idx(Derivatives::NUM)> jastrow_derivatives
) {
#if defined(XPU_CUDA)
  dim3 derivativeSumsThreads{256};
  dim3 derivativeSumsBlocks{
    xpu::block_per_dim(log_derivatives.stride(), derivativeSumsThreads.x)
  };
  
  cudaDerivativeSums<<<
    derivativeSumsBlocks, derivativeSumsThreads
  >>>(
    log_derivatives,
    jastrow_derivatives
  );

  xpu::cu_check(cudaGetLastError());
  xpu::cu_check(cudaDeviceSynchronize());
#else
  const auto N{log_derivatives.stride()};
  #pragma omp simd
  for (auto i = 0uz; i < N; ++i) {
    stencil::wavefunction::derivative_sum(i, log_derivatives, jastrow_derivatives);
  }
#endif
}

} // namespace kernel::wavefunction
} // namespace kernel