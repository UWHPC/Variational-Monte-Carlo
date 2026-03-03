#pragma once

#include "../config/config.hpp"
#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "../wavefunction/jastrow_pade.hpp"
#include "../wavefunction/wavefunction.hpp"

#include <random>

class Simulation {
private:
    Config config_;
    Particles particles_;
    PeriodicBoundaryCondition pbc_;
    JastrowPade jastrowPade_;
    WaveFunction waveFunction_;

    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> uniform01_{0.0, 1.0};
    std::uniform_real_distribution<double> proposal_;
    std::uniform_int_distribution<std::size_t> pickParticle_;

    std::size_t proposed_{0};
    std::size_t accepted_{0};
    double logPsiCurrent_{0.0};

public:
    explicit Simulation(Config cfg);
    void run();

private:
    void initializePositions();
    bool metropolisStep();
    void warmup();
    void measure();
};