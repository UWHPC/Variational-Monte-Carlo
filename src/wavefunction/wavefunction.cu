#include <xpu/xpu.hpp>
#include "wavefunction.hpp"
#include "wavefunction_kernels.hpp"

namespace {

void reset_derivatives(
  xpu::soa_view<real_t, idx(Derivatives::NUM)> derivative
) {
  const auto total_bytes{
    static_cast<std::size_t>(derivative.stride() * idx(Derivatives::NUM) * sizeof(real_t))
  };
  xpu::memset(derivative[0], 0, total_bytes);
}

} // namespace

void WaveFunction::evaluate_derivatives(Particles& particles) noexcept {
  reset_derivatives(particles.derivatives());
  reset_derivatives(this->j_derivatives());

  slater_plane_wave_.add_derivatives(particles.derivatives());
  jastrow_pade_.add_derivatives(particles, this->j_derivatives());

  kernel::wavefunction::derivative_sums(particles.derivatives(), this->j_derivatives());

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
  if (steps_since_refresh() >= 500) {
    evaluate_derivatives(particles);
    return;
  }
  set_steps_since_refresh(steps_since_refresh() + 1);
  if (!move_accepted) {
    return;
  }

  jastrow_pade_.update_derivatives_for_move(
    particles,
    moved,
    old_x, old_y, old_z,
    this->j_derivatives()
  );

  reset_derivatives(particles.derivatives());

  slater_plane_wave_.add_derivatives(
    particles.derivatives()
  );
  
  kernel::wavefunction::derivative_sums(particles.derivatives(), this->j_derivatives());
}
