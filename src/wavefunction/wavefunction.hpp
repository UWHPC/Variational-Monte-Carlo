#pragma once

#include "../jastrow_pade/jastrow_pade.hpp"
#include "../particles/particles.hpp"
#include "../slater_plane_wave/slater_plane_wave.hpp"

class WaveFunction {
private:
    JastrowPade jastrow_pade_;
    SlaterPlaneWave slater_plane_wave_;

public:
    explicit WaveFunction(std::size_t num_particles, double box_length, double a = 0.5, double b = 1.0)
        : jastrow_pade_{box_length, a, b}, slater_plane_wave_{num_particles, box_length} {}

    [[nodiscard]] JastrowPade& jastrow_pade_get() { return jastrow_pade_; }
    [[nodiscard]] const JastrowPade& jastrow_pade_get() const { return jastrow_pade_; }

    [[nodiscard]] SlaterPlaneWave& slater_plane_wave_get() { return slater_plane_wave_; }
    [[nodiscard]] const SlaterPlaneWave& slater_plane_wave_get() const { return slater_plane_wave_; }

    void evaluate_derivatives(Particles& particles) const noexcept;
    double evaluate_log_psi(const Particles& particles);
};
