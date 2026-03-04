#pragma once

#include "../config/config.hpp"
#include "../jastrow_pade/jastrow_pade.hpp"
#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "../slater_plane_wave/slater_plane_wave.hpp"
#include "../wavefunction/wavefunction.hpp"

#include <random>

class Simulation {
private:
    Config config_;
    Particles particles_;
    PeriodicBoundaryCondition pbc_;

    JastrowPade jastrow_pade_;
    SlaterPlaneWave slater_plane_wave_;
    WaveFunction wave_function_;

    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> uniform01_{0.0, 1.0};
    std::uniform_real_distribution<double> proposal_;
    std::uniform_int_distribution<std::size_t> pick_particle_;

    std::size_t proposed_{};
    std::size_t accepted_{};
    double log_psi_current_{};

public:
    explicit Simulation(Config cfg) noexcept;

    [[nodiscard]] std::mt19937_64& rng() { return rng_; }
    [[nodiscard]] std::uniform_real_distribution<double> uniform01() { return uniform01_; }
    [[nodiscard]] std::uniform_real_distribution<double>& proposal() { return proposal_; }
    [[nodiscard]] std::uniform_int_distribution<std::size_t>& pick_particle() { return pick_particle_; }
    [[nodiscard]] PeriodicBoundaryCondition pbc() { return pbc_; }
    [[nodiscard]] Particles& particles() { return particles_; }

    void run();

private:
    void initialize_positions();
    bool metropolis_step();
    void warmup();
    void measure();
   
};