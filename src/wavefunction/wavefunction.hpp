#pragma once

#include "../jastrow_pade/jastrow_pade.hpp"
#include "../particles/particles.hpp"
#include "../slater_plane_wave/slater_plane_wave.hpp"
#include "../utilities/aligned_soa.hpp"

#include <cstddef>
#include <vector>

class WaveFunction {
private:
    JastrowPade jastrow_pade_;
    SlaterPlaneWave slater_plane_wave_;

    bool jastrow_cache_valid_;        // whether cache is filled
    std::size_t steps_since_refresh_; // for drift correction

    // Aligned SOA for the gradieent components and laplacian:
    enum ArrayIndex : std::size_t { GRAD_X_, GRAD_Y_, GRAD_Z_, LAP_, NUM_ARRAYS_ };
    AlignedSoA<double> derivatives_;

public:
    explicit WaveFunction(const Particles& particles, double box_length, double a = 0.5, double b = 1.0)
        : jastrow_pade_{box_length, a, b}, slater_plane_wave_{particles, box_length}, jastrow_cache_valid_{},
          steps_since_refresh_{}, derivatives_{particles.num_particles_get(), NUM_ARRAYS_} {}

    [[nodiscard]] JastrowPade& jastrow_pade_get() { return jastrow_pade_; }
    [[nodiscard]] const JastrowPade& jastrow_pade_get() const { return jastrow_pade_; }

    [[nodiscard]] SlaterPlaneWave& slater_plane_wave_get() { return slater_plane_wave_; }
    [[nodiscard]] const SlaterPlaneWave& slater_plane_wave_get() const { return slater_plane_wave_; }

    [[nodiscard]] double* jastrow_grad_x_get() noexcept { return derivatives_[GRAD_X_]; }
    [[nodiscard]] double* jastrow_grad_y_get() noexcept { return derivatives_[GRAD_Y_]; }
    [[nodiscard]] double* jastrow_grad_z_get() noexcept { return derivatives_[GRAD_Z_]; }
    [[nodiscard]] double* jastrow_lap_get() noexcept { return derivatives_[LAP_]; }

    [[nodiscard]] bool jastrow_cache_valid_get() const noexcept { return jastrow_cache_valid_; }
    void jastrow_cache_valid_set(bool value) noexcept { jastrow_cache_valid_ = value; }

    [[nodiscard]] std::size_t steps_since_refresh_get() const noexcept { return steps_since_refresh_; }
    void steps_since_refresh_set(std::size_t value) noexcept { steps_since_refresh_ = value; }

    void evaluate_derivatives(Particles& particles) noexcept;
    void evaluate_derivatives(Particles& particles, bool move_accepted, std::size_t moved, double old_x, double old_y,
                              double old_z) noexcept;
    double evaluate_log_psi(const Particles& particles);
};
