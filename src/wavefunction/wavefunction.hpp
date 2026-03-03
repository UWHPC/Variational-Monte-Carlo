#pragma once

#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "jastrow_pade.hpp"

class WaveFunction {
private:
    JastrowPade jastrowPade_;

public:
    explicit WaveFunction(const JastrowPade& jastrowPade) noexcept
    : jastrowPade_{jastrowPade}
    { }

    void evaluateLogPsi(
        Particles& particles, 
        const PeriodicBoundaryCondition& pbc
    ) const noexcept;

    void evaluateDerivatives(
        Particles& particles, 
        const PeriodicBoundaryCondition& pbc
    ) const noexcept;
};