#pragma once

#include "../jastrow_pade/jastrow_pade.hpp"
#include "../particles/particles.hpp"
#include "../slater_plane_wave/slater_plane_wave.hpp"
#include <cstddef>
#include <vector>

class WaveFunction {
private:
    JastrowPade jastrow_pade_;
    SlaterPlaneWave slater_plane_wave_;

    std::vector<double> jastrow_grad_x_;
    std::vector<double> jastrow_grad_y_;
    std::vector<double> jastrow_grad_z_;
    std::vector<double> jastrow_lap_;
    bool jastrow_cache_valid_;        // whether cache is filled
    std::size_t steps_since_refresh_; // for drift correction

public:
    explicit WaveFunction(std::size_t num_particles, double box_length, double a = 0.5, double b = 1.0)
        : jastrow_pade_{box_length, a, b}, slater_plane_wave_{num_particles, box_length},
          jastrow_grad_x_(num_particles, 0.0), jastrow_grad_y_(num_particles, 0.0), jastrow_grad_z_(num_particles, 0.0),
          jastrow_lap_(num_particles, 0.0), jastrow_cache_valid_{false}, steps_since_refresh_{0} {}

    [[nodiscard]] JastrowPade& jastrow_pade_get() { return jastrow_pade_; }
    [[nodiscard]] const JastrowPade& jastrow_pade_get() const { return jastrow_pade_; }

    [[nodiscard]] SlaterPlaneWave& slater_plane_wave_get() { return slater_plane_wave_; }
    [[nodiscard]] const SlaterPlaneWave& slater_plane_wave_get() const { return slater_plane_wave_; }

    [[nodiscard]] double* jastrow_grad_x_get() noexcept { return jastrow_grad_x_.data(); }
    [[nodiscard]] double* jastrow_grad_y_get() noexcept { return jastrow_grad_y_.data(); }
    [[nodiscard]] double* jastrow_grad_z_get() noexcept { return jastrow_grad_z_.data(); }
    [[nodiscard]] double* jastrow_lap_get() noexcept { return jastrow_lap_.data(); }

    [[nodiscard]] bool jastrow_cache_valid_get() const noexcept { return jastrow_cache_valid_; }
    void jastrow_cache_valid_set(bool value) noexcept { jastrow_cache_valid_ = value; }

    [[nodiscard]] std::size_t steps_since_refresh_get() const noexcept { return steps_since_refresh_; }
    void steps_since_refresh_set(std::size_t value) noexcept { steps_since_refresh_ = value; }

    void evaluate_derivatives(Particles& particles) noexcept;
    void evaluate_derivatives(Particles& particles, bool move_accepted, std::size_t moved, double old_x, double old_y, double old_z) noexcept;
    double evaluate_log_psi(const Particles& particles);
};
