#pragma once

#include "../blocking_analysis/blocking_analysis.hpp"
#include "../config/config.hpp"
#include "../output_writer/output_writer.hpp"
#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "../wavefunction/wavefunction.hpp"

#include <memory>
#include <optional>
#include <random>
#include <vector>

class Simulation {
private:
    // Configuration of sim:
    Config config_;

    // Objects and physical handlers of sim:
    Particles particles_;
    PeriodicBoundaryCondition pbc_;
    WaveFunction wave_function_;
    BlockingAnalysis blocking_analysis_;
    std::unique_ptr<OutputWriter> output_writer_;

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

    [[nodiscard]] PeriodicBoundaryCondition& pbc_get() { return pbc_; }
    [[nodiscard]] Particles& particles_get() { return particles_; }
    [[nodiscard]] WaveFunction& wave_function_get() { return wave_function_; }
    [[nodiscard]] BlockingAnalysis& blocking_analysis_get() { return blocking_analysis_; }

    // Randomly generated uniform, proposal, and particle:
    [[nodiscard]] double rand_uniform_double() { return uniform01()(rng()); }
    [[nodiscard]] double rand_proposal_double() { return proposal()(rng()); }
    [[nodiscard]] std::size_t rand_particle_get() { return pick_particle()(rng()); }

    [[nodiscard]] double& log_psi_set() { return log_psi_current_; }
    [[nodiscard]] double log_psi_get() const { return log_psi_current_; }

    [[nodiscard]] double acceptance_rate() const {
        if (proposed_ == 0U) {
            return 0.0;
        }
        return static_cast<double>(accepted_) / static_cast<double>(proposed_);
    }

    struct MeasurementSummary {
        double mean_energy;
        std::optional<double> standard_error;
    };

    [[nodiscard]] std::vector<double> positions_snapshot() const;

public:
    explicit Simulation(Config cfg, std::unique_ptr<OutputWriter> output_writer = nullptr) noexcept;
    void run();

private:
    double kinetic_energy(const Particles& particles) const noexcept;
    double potential_energy(const Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept;
    double eval_total_energy(const Particles& p, const PeriodicBoundaryCondition& pbc) const noexcept;

    void initialize_positions();
    bool metropolis_step();
    void warmup();
    MeasurementSummary measure();
};
