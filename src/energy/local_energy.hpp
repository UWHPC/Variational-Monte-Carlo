#pragma once

#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"

class EnergyTracker {
private:
    double kinetic_energy(const Particles& particles) const noexcept;
    double potential_energy(const Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept;

public:
    double eval_total_energy(const Particles& p, const PeriodicBoundaryCondition& pbc) const noexcept;
};