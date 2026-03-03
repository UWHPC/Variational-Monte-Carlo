#pragma once
#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "jastrow_pade.hpp"

class WaveFunction {
private:
    JastrowPade j_;

public:
    explicit WaveFunction(JastrowPade j) noexcept : j_(j) {};

    void evaluateLogPsi(Particles<>& particles, const periodicBoundaryCondition& pbc) const noexcept;

    void evaluateDerivatives(Particles<>& particles, const periodicBoundaryCondition& pbc) const noexcept;
};