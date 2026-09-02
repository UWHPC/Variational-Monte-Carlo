#include <xpu/xpu.hpp>
#include "wavefunction.hpp"
#include "wavefunction_kernels.hpp"

namespace {

void reset_derivatives(
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> derivative
) {
  const auto total_bytes{
    static_cast<std::size_t>(derivative.stride() * idx(Derivatives::NUM) * sizeof(fp_t))
  };
  xpu::memset(derivative[0], 0, total_bytes);
}

} // namespace

void WaveFunction::evaluate_derivatives(Particles::View particles, std::size_t walker) noexcept {
  reset_derivatives(particles.derivatives);
  reset_derivatives(this->j_derivatives(walker));

  slater_plane_wave_.add_derivatives(particles, walker);
  jastrow_pade_.add_derivatives(particles, this->j_derivatives(walker));

  kernel::wavefunction::derivative_sums(this->view(walker), particles);

  set_jastrow_cache_valid(true, walker);
  set_steps_since_refresh(0uz, walker);
}

void WaveFunction::evaluate_derivatives(
  Particles::View particles,
  bool move_accepted,
  std::size_t moved,
  xpu::array<fp_t, idx(Axis::NUM)> old_pos,
  std::size_t walker
) noexcept {
  if (!jastrow_cache_valid(walker)) {
    evaluate_derivatives(particles, walker);
    return;
  }
  if (steps_since_refresh(walker) >= 500) {
    evaluate_derivatives(particles, walker);
    return;
  }
  set_steps_since_refresh(steps_since_refresh(walker) + 1uz, walker);
  if (!move_accepted) {
    return;
  }

  jastrow_pade_.update_derivatives_for_move(
    moved,
    old_pos,
    particles,
    this->j_derivatives(walker)
  );

  reset_derivatives(particles.derivatives);

  slater_plane_wave_.add_derivatives(
    particles,
    walker
  );
  
  kernel::wavefunction::derivative_sums(this->view(walker), particles);
}
