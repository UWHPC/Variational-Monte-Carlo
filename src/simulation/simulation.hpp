#pragma once

#include "../config/config.hpp"
#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "../wavefunction/jastrow_pade.hpp"
#include "../wavefunction/wavefunction.hpp"

class Simulation {
private:
    Config config_;
    Particles particles_;
    PeriodicBoundaryCondition pbc_;
    JastrowPade jastrowPade_;
    WaveFunction waveFunction_;

public:
    explicit Simulation(Config config) noexcept;
};