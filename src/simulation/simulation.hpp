#pragma once

#include "../particles/particles.hpp"
#include "../config/config.hpp"
#include "../pbc/pbc.hpp"
#include "../jastrow_pade/jastrow_pade.hpp"
#include "../slater_plane_wave/slater_plane_wave.hpp"
#include "../wavefunction/wavefunction.hpp"

#include <random>

class Simulation {
private:
    Config config_;
    Particles particles_;
    PeriodicBoundaryCondition pbc_;

    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> uniform01_{0.0, 1.0};
    std::uniform_real_distribution<double> proposal_;
    std::uniform_int_distribution<std::size_t> pickParticle_;

    std::size_t proposed_{};
    std::size_t accepted_{};
    double logPsiCurrent_{};

public:
    explicit Simulation(Config cfg) noexcept;
    void run();

private:
    void initializePositions();
    bool metropolisStep();
    void warmup();
    void measure();
};