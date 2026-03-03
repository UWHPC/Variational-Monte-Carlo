#pragma once

#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "jastrow_pade.hpp"
#include "slater_plane_wave.hpp"

class WaveFunction {
private:
    const JastrowPade& jastrowPade_;
    SlaterPlaneWave& slaterPlaneWave_; // holds a non-const reference because logsAbsDet catches LU/inv
public:
    explicit WaveFunction(const JastrowPade& jastrowPade, SlaterPlaneWave& slaterPlaneWave) noexcept
    : jastrowPade_{jastrowPade}, slaterPlaneWave_{slaterPlaneWave} {};

    void evaluateLogPsi(Particles& particles, const periodicBoundaryCondition& pbc); // not const noexcept because it calls slaterPlaneWave_.logAbsDet();

    void evaluateDerivatives(Particles& particles, const periodicBoundaryCondition& pbc) const noexcept;
};