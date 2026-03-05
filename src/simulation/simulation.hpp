#pragma once

#include "../config/config.hpp"
#include "../energy/local_energy.hpp"
#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "../wavefunction/wavefunction.hpp"

#include <random>

class Simulation {
private:
    // Configuration of sim:
    Config config_;

    // Objects and physical handlers of sim:
    Particles particles_;
    PeriodicBoundaryCondition pbc_;
    WaveFunction wave_function_;
    EnergyTracker energy_tracker_;

    // Physical quantities of sim:
    std::size_t proposed_;
    std::size_t accepted_;
    double log_psi_current_;

    // Random num generation:
    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> uniform01_{0.0, 1.0};
    std::uniform_real_distribution<double> proposal_;
    std::uniform_int_distribution<std::size_t> pick_particle_;

    // Getters:
    [[nodiscard]] std::mt19937_64& rng() { return rng_; }
    [[nodiscard]] std::uniform_real_distribution<double>& uniform01() { return uniform01_; }
    [[nodiscard]] std::uniform_real_distribution<double>& proposal() { return proposal_; }
    [[nodiscard]] std::uniform_int_distribution<std::size_t>& pick_particle() { return pick_particle_; }

    // Randomly generated uniform, proposal, and particle:
    [[nodiscard]] double rand_uniform_double() { return uniform01()(rng()); }
    [[nodiscard]] double rand_proposal_double() { return proposal()(rng()); }
    [[nodiscard]] std::size_t rand_particle_get() { return pick_particle()(rng()); }

public:
    explicit Simulation(Config cfg) noexcept;

    // Getters for objects:
    [[nodiscard]] PeriodicBoundaryCondition& pbc() { return pbc_; }
    [[nodiscard]] Particles& particles() { return particles_; }
    [[nodiscard]] WaveFunction& wave_function() { return wave_function_; }

    // Getters:
    [[nodiscard]] double& log_psi_set() { return log_psi_current_; }
    [[nodiscard]] double log_psi_get() const { return log_psi_current_; }

    void run();

private:
    void initialize_positions();
    bool metropolis_step();
    void warmup();
    void measure();
};
