#include "wavefunction.hpp"

#include <algorithm>

void WaveFunction::evaluate_log_psi(Particles& particles, const PeriodicBoundaryCondition& pbc) {
    const double log_det{slater_plane_wave_ptr().log_abs_det(particles)};
    const double jastrow_pade{jastrow_pade_ptr().value(particles, pbc)};

    particles.log_psi_ptr()[0] = log_det + jastrow_pade;
}

void WaveFunction::evaluate_derivatives(Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept {
    const std::size_t padded_stride{particles.padding_stride_ptr()};

    std::fill_n(particles.grad_log_psi_x_ptr(), padded_stride, 0.0);
    std::fill_n(particles.grad_log_psi_y_ptr(), padded_stride, 0.0);
    std::fill_n(particles.grad_log_psi_z_ptr(), padded_stride, 0.0);
    std::fill_n(particles.lap_log_psi_ptr(), padded_stride, 0.0);

    slater_plane_wave_.add_derivatives(particles, particles.grad_log_psi_x_ptr(), particles.grad_log_psi_y_ptr(),
                                       particles.grad_log_psi_z_ptr(), particles.lap_log_psi_ptr());

    jastrow_pade_.add_derivatives(particles, pbc, particles.grad_log_psi_x_ptr(), particles.grad_log_psi_y_ptr(),
                                  particles.grad_log_psi_z_ptr(), particles.lap_log_psi_ptr());
}
