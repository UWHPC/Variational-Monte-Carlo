#include "wavefunction.hpp"

#include <algorithm>

real_t WaveFunction::evaluate_log_psi(const Particles& particles) {
  const real_t log_det{slater_plane_wave().log_abs_det(particles)};
  const real_t jastrow{jastrow_pade().value(particles)};

  return log_det + jastrow;
}

void WaveFunction::evaluate_derivatives(Particles& particles) noexcept {
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
  std::size_t steps{steps_since_refresh()};
  if (steps >= 500) {
    evaluate_derivatives(particles);
    return;
  }
  set_steps_since_refresh(steps + 1);
  if (!move_accepted) {
    return;
  }
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
}
