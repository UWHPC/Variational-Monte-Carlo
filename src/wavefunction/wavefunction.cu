#include <xpu/xpu.hpp>
#include "wavefunction.cuh"

#include <algorithm>

real_t WaveFunction::evaluate_log_psi(const Particles& particles) {
  return this->slater_plane_wave().log_abs_det(particles) + this->jastrow_pade().value(particles);
}

#if defined(VMC_CUDA_BACKEND)

namespace {

__global__
void cudaDerivativeSums(
  std::size_t p_N,
  real_t* log_gx, real_t* log_gy, real_t* log_gz, real_t* log_lap,
  real_t* jp_gx, real_t* jp_gy, real_t* jp_gz, real_t* jp_lap
) {
  const auto [i]{xpu::global_index<1>()};
  if (i >= p_N) { return; }

  log_gx[i] += jp_gx[i];
  log_gy[i] += jp_gy[i];
  log_gz[i] += jp_gz[i];
  log_lap[i] += jp_lap[i];
}

void cudaFillDerivatives(
  std::size_t total_bytes,
  real_t* gx, real_t* gy, real_t* gz,
  real_t* lap
) {
  CUDA_CHECK(cudaMemset(gx, 0, total_bytes));
  CUDA_CHECK(cudaMemset(gy, 0, total_bytes));
  CUDA_CHECK(cudaMemset(gz, 0, total_bytes));
  CUDA_CHECK(cudaMemset(lap, 0, total_bytes));
}

}

#endif

void WaveFunction::evaluate_derivatives(Particles& particles) noexcept {
#if defined(VMC_CUDA_BACKEND)
  const auto total_bytes{
    static_cast<std::size_t>(particles.p_stride() * sizeof(real_t))
  };

  cudaFillDerivatives(
    total_bytes,
    particles.grad_log_psi().x_,
    particles.grad_log_psi().y_,
    particles.grad_log_psi().z_,
    particles.lap_log_psi()
  );

  cudaFillDerivatives(
    total_bytes,
    this->j_grad().x_,
    this->j_grad().y_,
    this->j_grad().z_,
    this->j_lap()
  );

  slater_plane_wave_.add_derivatives(
    particles.grad_log_psi().x_,
    particles.grad_log_psi().y_,
    particles.grad_log_psi().z_,
    particles.lap_log_psi()
  );

  jastrow_pade_.add_derivatives(
    particles,
    this->j_grad().x_,
    this->j_grad().y_,
    this->j_grad().z_,
    this->j_lap()
  );

  dim3 derivativeSumsThreads(256);
  dim3 derivativeSumsBlocks(
    vmc::cudaNumBlocks(particles.p_stride(), derivativeSumsThreads.x)
  );
  
  cudaDerivativeSums<<<derivativeSumsBlocks, derivativeSumsThreads>>>(
    particles.p_stride(),
    particles.grad_log_psi().x_,
    particles.grad_log_psi().y_,
    particles.grad_log_psi().z_,
    particles.lap_log_psi(),
    this->j_grad().x_,
    this->j_grad().y_,
    this->j_grad().z_,
    this->j_lap()
  );

  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  set_jastrow_cache_valid(true);
  set_steps_since_refresh(0);
#else
  const std::size_t padded_stride{particles.p_stride()};

  const auto log_grad{particles.grad_log_psi().align()};
  real_t* RESTRICT log_lap{particles.lap_log_psi()};

  const auto jg{j_grad().align()};
  real_t* RESTRICT jastrow_lap{j_lap()};

  ASSUME_ALIGNED(log_lap, SIMD_BYTES);
  ASSUME_ALIGNED(jastrow_lap, SIMD_BYTES);

  std::fill_n(log_grad.x_, padded_stride, 0.0_r);
  std::fill_n(log_grad.y_, padded_stride, 0.0_r);
  std::fill_n(log_grad.z_, padded_stride, 0.0_r);
  std::fill_n(log_lap, padded_stride, 0.0_r);

  std::fill_n(jg.x_, padded_stride, 0.0_r);
  std::fill_n(jg.y_, padded_stride, 0.0_r);
  std::fill_n(jg.z_, padded_stride, 0.0_r);
  std::fill_n(jastrow_lap, padded_stride, 0.0_r);

  slater_plane_wave_.add_derivatives(log_grad.x_, log_grad.y_, log_grad.z_, log_lap);
  jastrow_pade_.add_derivatives(particles, jg.x_, jg.y_, jg.z_, jastrow_lap);

  #pragma omp simd
  for (std::size_t i = 0; i < padded_stride; ++i) {
    log_grad.x_[i] += jg.x_[i];
    log_grad.y_[i] += jg.y_[i];
    log_grad.z_[i] += jg.z_[i];
    log_lap[i] += jastrow_lap[i];
  }
  set_jastrow_cache_valid(true);
  set_steps_since_refresh(0);
#endif
}

void WaveFunction::evaluate_derivatives(
  Particles& particles,
  bool move_accepted,
  std::size_t moved,
  real_t old_x, real_t old_y, real_t old_z
) noexcept {
  if (!jastrow_cache_valid()) {
    evaluate_derivatives(particles);
    return;
  }
  if (steps_since_refresh() >= 500) {
    evaluate_derivatives(particles);
    return;
  }
  set_steps_since_refresh(steps_since_refresh() + 1);
  if (!move_accepted) {
    return;
  }

#if defined(VMC_CUDA_BACKEND)
  const auto total_bytes{
    static_cast<std::size_t>(particles.p_stride() * sizeof(real_t))
  };

  jastrow_pade_.update_derivatives_for_move(
    particles,
    moved,
    old_x, old_y, old_z,
    this->j_grad().x_,
    this->j_grad().y_,
    this->j_grad().z_,
    this->j_lap()
  );

  cudaFillDerivatives(
    total_bytes,
    particles.grad_log_psi().x_,
    particles.grad_log_psi().y_,
    particles.grad_log_psi().z_,
    particles.lap_log_psi()
  );

  slater_plane_wave_.add_derivatives(
    particles.grad_log_psi().x_,
    particles.grad_log_psi().y_,
    particles.grad_log_psi().z_,
    particles.lap_log_psi()
  );

  dim3 derivativeSumsThreads(256);
  dim3 derivativeSumsBlocks(
    vmc::cudaNumBlocks(particles.p_stride(), derivativeSumsThreads.x)
  );
  
  cudaDerivativeSums<<<derivativeSumsBlocks, derivativeSumsThreads>>>(
    particles.p_stride(),
    particles.grad_log_psi().x_,
    particles.grad_log_psi().y_,
    particles.grad_log_psi().z_,
    particles.lap_log_psi(),
    this->j_grad().x_,
    this->j_grad().y_,
    this->j_grad().z_,
    this->j_lap()
  );

  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
#else
  const std::size_t padded_stride{particles.p_stride()};

  const auto log_grad{particles.grad_log_psi().align()};
  real_t* RESTRICT log_lap{particles.lap_log_psi()};

  const auto jg{j_grad().align()};
  real_t* RESTRICT jastrow_lap{j_lap()};

  ASSUME_ALIGNED(log_lap, SIMD_BYTES);
  ASSUME_ALIGNED(jastrow_lap, SIMD_BYTES);

  jastrow_pade_.update_derivatives_for_move(
    particles,
    moved,
    old_x, old_y, old_z,
    jg.x_, jg.y_, jg.z_,
    jastrow_lap
  );

  std::fill_n(log_grad.x_, padded_stride, 0.0_r);
  std::fill_n(log_grad.y_, padded_stride, 0.0_r);
  std::fill_n(log_grad.z_, padded_stride, 0.0_r);
  std::fill_n(log_lap, padded_stride, 0.0_r);

  slater_plane_wave_.add_derivatives(
    log_grad.x_, log_grad.y_, log_grad.z_,
    log_lap
  );

  #pragma omp simd
  for (std::size_t i = 0; i < padded_stride; ++i) {
    log_grad.x_[i] += jg.x_[i];
    log_grad.y_[i] += jg.y_[i];
    log_grad.z_[i] += jg.z_[i];
    log_lap[i] += jastrow_lap[i];
  }
#endif
}
