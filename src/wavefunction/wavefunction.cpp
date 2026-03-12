#include "wavefunction.hpp"

#include <algorithm>

double WaveFunction::evaluate_log_psi(const Particles& particles) {
    const double log_det{slater_plane_wave_get().log_abs_det(particles)};
    const double jastrow_pade{jastrow_pade_get().value(particles)};

    return log_det + jastrow_pade;
}

void WaveFunction::evaluate_derivatives(Particles& particles) noexcept {
    const std::size_t padded_stride{particles.padding_stride_get()};

    std::fill_n(particles.grad_log_psi_x_get(), padded_stride, 0.0);
    std::fill_n(particles.grad_log_psi_y_get(), padded_stride, 0.0);
    std::fill_n(particles.grad_log_psi_z_get(), padded_stride, 0.0);
    std::fill_n(particles.lap_log_psi_get(), padded_stride, 0.0);

    const std::size_t N{particles.num_particles_get()};
    std::fill_n(jastrow_grad_x_get(), N, 0.0);
    std::fill_n(jastrow_grad_y_get(), N, 0.0);
    std::fill_n(jastrow_grad_z_get(), N, 0.0);
    std::fill_n(jastrow_lap_get(), N, 0.0);

    slater_plane_wave_.add_derivatives(particles, particles.grad_log_psi_x_get(), particles.grad_log_psi_y_get(),
                                       particles.grad_log_psi_z_get(), particles.lap_log_psi_get());

    jastrow_pade_.add_derivatives(particles, jastrow_grad_x_get(), jastrow_grad_y_get(), jastrow_grad_z_get(),
                                  jastrow_lap_get());

    for (std::size_t i{}; i < N; ++i) {
        particles.grad_log_psi_x_get()[i] += jastrow_grad_x_get()[i];
        particles.grad_log_psi_y_get()[i] += jastrow_grad_y_get()[i];
        particles.grad_log_psi_z_get()[i] += jastrow_grad_z_get()[i];
        particles.lap_log_psi_get()[i] += jastrow_lap_get()[i];
    }
    jastrow_cache_valid_set(true);
    steps_since_refresh_set(0);
}

void WaveFunction::evaluate_derivatives(Particles& particles, bool move_accepted, std::size_t moved, double old_x,
                                        double old_y, double old_z) noexcept {
    if (!jastrow_cache_valid_get()) {
        evaluate_derivatives(particles);
        return;
    }
    std::size_t steps = steps_since_refresh_get();
    if (steps >= 500) {
        evaluate_derivatives(particles);
        return;
    }
    steps_since_refresh_set(steps + 1);
    if (move_accepted) {
        jastrow_pade_.update_derivatives_for_move(particles, moved, old_x, old_y, old_z, jastrow_grad_x_get(),
                                                  jastrow_grad_y_get(), jastrow_grad_z_get(), jastrow_lap_get());
    }
    const std::size_t padded_stride{particles.padding_stride_get()};

    std::fill_n(particles.grad_log_psi_x_get(), padded_stride, 0.0);
    std::fill_n(particles.grad_log_psi_y_get(), padded_stride, 0.0);
    std::fill_n(particles.grad_log_psi_z_get(), padded_stride, 0.0);
    std::fill_n(particles.lap_log_psi_get(), padded_stride, 0.0);

    const std::size_t N{particles.num_particles_get()};

    slater_plane_wave_.add_derivatives(particles, particles.grad_log_psi_x_get(), particles.grad_log_psi_y_get(),
                                       particles.grad_log_psi_z_get(), particles.lap_log_psi_get());

    for (std::size_t i{}; i < N; ++i) {
        particles.grad_log_psi_x_get()[i] += jastrow_grad_x_get()[i];
        particles.grad_log_psi_y_get()[i] += jastrow_grad_y_get()[i];
        particles.grad_log_psi_z_get()[i] += jastrow_grad_z_get()[i];
        particles.lap_log_psi_get()[i] += jastrow_lap_get()[i];
    }
}
