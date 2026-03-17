#include "wavefunction.hpp"

#include <algorithm>

double WaveFunction::evaluate_log_psi(const Particles& particles) {
    const double log_det{slater_plane_wave_get().log_abs_det(particles)};
    const double jastrow_pade{jastrow_pade_get().value(particles)};

    return log_det + jastrow_pade;
}

void WaveFunction::evaluate_derivatives(Particles& particles) noexcept {
    const std::size_t padded_stride{particles.padding_stride_get()};

    double* RESTRICT log_grad_x{particles.grad_log_psi_x_get()};
    double* RESTRICT log_grad_y{particles.grad_log_psi_y_get()};
    double* RESTRICT log_grad_z{particles.grad_log_psi_z_get()};
    double* RESTRICT log_lap{particles.lap_log_psi_get()};

    double* RESTRICT jastrow_grad_x{jastrow_grad_x_get()};
    double* RESTRICT jastrow_grad_y{jastrow_grad_y_get()};
    double* RESTRICT jastrow_grad_z{jastrow_grad_z_get()};
    double* RESTRICT jastrow_lap{jastrow_lap_get()};

    ASSUME_ALIGNED(log_grad_x, SIMD_BYTES);
    ASSUME_ALIGNED(log_grad_y, SIMD_BYTES);
    ASSUME_ALIGNED(log_grad_z, SIMD_BYTES);
    ASSUME_ALIGNED(log_lap, SIMD_BYTES);

    ASSUME_ALIGNED(jastrow_grad_x, SIMD_BYTES);
    ASSUME_ALIGNED(jastrow_grad_y, SIMD_BYTES);
    ASSUME_ALIGNED(jastrow_grad_z, SIMD_BYTES);
    ASSUME_ALIGNED(jastrow_lap, SIMD_BYTES);

    std::fill_n(log_grad_x, padded_stride, 0.0);
    std::fill_n(log_grad_y, padded_stride, 0.0);
    std::fill_n(log_grad_z, padded_stride, 0.0);
    std::fill_n(log_lap, padded_stride, 0.0);

    std::fill_n(jastrow_grad_x, padded_stride, 0.0);
    std::fill_n(jastrow_grad_y, padded_stride, 0.0);
    std::fill_n(jastrow_grad_z, padded_stride, 0.0);
    std::fill_n(jastrow_lap, padded_stride, 0.0);

    slater_plane_wave_.add_derivatives(log_grad_x, log_grad_y, log_grad_z, log_lap);
    jastrow_pade_.add_derivatives(particles, jastrow_grad_x, jastrow_grad_y, jastrow_grad_z,
                                  jastrow_lap);

#pragma omp simd
    for (std::size_t i = 0; i < padded_stride; ++i) {
        log_grad_x[i] += jastrow_grad_x[i];
        log_grad_y[i] += jastrow_grad_y[i];
        log_grad_z[i] += jastrow_grad_z[i];
        log_lap[i] += jastrow_lap[i];
    }
    jastrow_cache_valid_set(true);
    steps_since_refresh_set(0);
}

void WaveFunction::evaluate_derivatives(Particles& particles, bool move_accepted, std::size_t moved,
                                        double old_x, double old_y, double old_z) noexcept {
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
    if (!move_accepted) {
        return;
    }
    const std::size_t padded_stride{particles.padding_stride_get()};

    double* RESTRICT log_grad_x{particles.grad_log_psi_x_get()};
    double* RESTRICT log_grad_y{particles.grad_log_psi_y_get()};
    double* RESTRICT log_grad_z{particles.grad_log_psi_z_get()};
    double* RESTRICT log_lap{particles.lap_log_psi_get()};

    double* RESTRICT jastrow_grad_x{jastrow_grad_x_get()};
    double* RESTRICT jastrow_grad_y{jastrow_grad_y_get()};
    double* RESTRICT jastrow_grad_z{jastrow_grad_z_get()};
    double* RESTRICT jastrow_lap{jastrow_lap_get()};

    ASSUME_ALIGNED(log_grad_x, SIMD_BYTES);
    ASSUME_ALIGNED(log_grad_y, SIMD_BYTES);
    ASSUME_ALIGNED(log_grad_z, SIMD_BYTES);
    ASSUME_ALIGNED(log_lap, SIMD_BYTES);

    ASSUME_ALIGNED(jastrow_grad_x, SIMD_BYTES);
    ASSUME_ALIGNED(jastrow_grad_y, SIMD_BYTES);
    ASSUME_ALIGNED(jastrow_grad_z, SIMD_BYTES);
    ASSUME_ALIGNED(jastrow_lap, SIMD_BYTES);

    jastrow_pade_.update_derivatives_for_move(particles, moved, old_x, old_y, old_z, jastrow_grad_x,
                                              jastrow_grad_y, jastrow_grad_z, jastrow_lap);

    std::fill_n(log_grad_x, padded_stride, 0.0);
    std::fill_n(log_grad_y, padded_stride, 0.0);
    std::fill_n(log_grad_z, padded_stride, 0.0);
    std::fill_n(log_lap, padded_stride, 0.0);

    slater_plane_wave_.add_derivatives(log_grad_x, log_grad_y, log_grad_z, log_lap);

#pragma omp simd
    for (std::size_t i = 0; i < padded_stride; ++i) {
        log_grad_x[i] += jastrow_grad_x[i];
        log_grad_y[i] += jastrow_grad_y[i];
        log_grad_z[i] += jastrow_grad_z[i];
        log_lap[i] += jastrow_lap[i];
    }
}
