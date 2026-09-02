#pragma once

#include <xpu/xpu.hpp>
#include "wavefunction.hpp"
#include "../jastrow_pade/jastrow_pade_kernels.hpp"
#include "../slater_plane_wave/slater_plane_wave_kernels.hpp"
#include "../utilities/execution.hpp"
#include "../utilities/components.hpp"
#include "../utilities/macros.hpp"

namespace stencil {
namespace wavefunction {

CUDA_CALLABLE
inline void derivative_sum(
  std::size_t i,
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> log_derivatives,
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> jastrow_derivatives
) {
  for (auto d{idx(Derivatives::GRAD_X)}; d < idx(Derivatives::NUM); ++d) {
    log_derivatives[d][i] += jastrow_derivatives[d][i];
  }
}

CUDA_CALLABLE
inline void evaluate_derivatives(
  WaveFunction::View wave_function,
  Particles::View particles
) {
  for (auto derivative{0uz}; derivative < idx(Derivatives::NUM); ++derivative) {
    for (auto particle{execution::thread()}; particle < particles.count; particle += execution::stride()) {
      particles.derivatives[derivative][particle] = 0.0_fp;
      wave_function.jastrow_derivatives[derivative][particle] = 0.0_fp;
    }
  }
  execution::sync();

  const auto slater_elements{
    wave_function.slater.num_orbitals * wave_function.slater.num_orbitals
  };
  for (auto element{execution::thread()}; element < slater_elements; element += execution::stride()) {
    const auto particle{element / wave_function.slater.num_orbitals};
    const auto orbital{element - particle * wave_function.slater.num_orbitals};
    stencil::slater::add_derivatives(
      particle,
      orbital,
      wave_function.slater,
      particles
    );
  }
  execution::sync();

  for (auto particle{execution::thread()}; particle < particles.count; particle += execution::stride()) {
    stencil::slater::accumulate_derivatives(particle, particles.derivatives);
  }
  execution::sync();

  const auto jastrow_elements{particles.count * particles.count};
  for (auto element{execution::thread()}; element < jastrow_elements; element += execution::stride()) {
    const auto particle{element / particles.count};
    const auto other{element - particle * particles.count};
    stencil::jpade::add_derivatives(
      particle,
      other,
      wave_function.jastrow,
      particles,
      wave_function.jastrow_derivatives
    );
  }
  execution::sync();

  for (auto particle{execution::thread()}; particle < particles.count; particle += execution::stride()) {
    derivative_sum(
      particle,
      particles.derivatives,
      wave_function.jastrow_derivatives
    );
  }

  if (execution::thread() == 0uz) {
    *wave_function.jastrow_cache_valid = 1u;
    *wave_function.steps_since_refresh = 0uz;
  }
  execution::sync();
}

} // namespace stencil::wavefunction
} // namespace stencil

namespace kernel {
namespace wavefunction {

namespace {

#if defined(XPU_CUDA)
__global__
void cudaDerivativeSums(
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> log_derivatives,
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> jastrow_derivatives
) {
  const auto [i]{xpu::global_index<1>()};
  if (i >= log_derivatives.stride()) { return; }

  stencil::wavefunction::derivative_sum(i,log_derivatives, jastrow_derivatives);
}
#endif

} // namespace

inline void derivative_sums(
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> log_derivatives,
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> jastrow_derivatives
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

inline void derivative_sums(
  WaveFunction::View wave_function,
  Particles::View particles
) {
  derivative_sums(
    particles.derivatives,
    wave_function.jastrow_derivatives
  );
}

} // namespace kernel::wavefunction
} // namespace kernel
