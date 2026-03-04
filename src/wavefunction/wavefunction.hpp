#pragma once

#include "../jastrow_pade/jastrow_pade.hpp"
#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "../slater_plane_wave/slater_plane_wave.hpp"

class WaveFunction {
private:
    const JastrowPade& jastrow_pade_;
    SlaterPlaneWave& slater_plane_wave_; // holds a non-const reference because logsAbsDet catches LU/inv
public:
    explicit WaveFunction(const JastrowPade& jastrow_pade, SlaterPlaneWave& slater_plane_wave) noexcept
        : jastrow_pade_{jastrow_pade}, slater_plane_wave_{slater_plane_wave} {}

    void evaluate_log_psi(
        Particles& particles,
        const PeriodicBoundaryCondition& pbc); // not const noexcept because it calls slaterPlaneWave_.logAbsDet();

    [[nodiscard]] const JastrowPade& jastrow_pade_ptr() const { return jastrow_pade_; }
    [[nodiscard]] SlaterPlaneWave& slater_plane_wave_ptr() const { return slater_plane_wave_; }

    void evaluateDerivatives(Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept;
};
