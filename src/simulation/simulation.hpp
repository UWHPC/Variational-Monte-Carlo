#pragma once

#include "../blocking_analysis/blocking_analysis.hpp"
#include "../config/config.hpp"
#include "../energy_tracking/energy_tracking.hpp"
#include "../output_writer/output_writer.hpp"
#include "../particles/particles.hpp"
#include "../wavefunction/wavefunction.hpp"

#include <memory>
#include <optional>
#include <random>
#include <vector>

class Simulation {
private:
  Config config_;

  Particles particles_;
  WaveFunction wave_function_;
  BlockingAnalysis blocking_analysis_;
  EnergyTracker energy_tracker_;
  std::unique_ptr<OutputWriter> output_writer_;

  std::size_t proposed_;
  std::size_t accepted_;
  real_t log_psi_current_;

  std::mt19937_64 rng_;
  std::uniform_real_distribution<real_t> uniform01_{0.0_r, 1.0_r};
  std::uniform_real_distribution<real_t> proposal_;
  std::uniform_int_distribution<std::size_t> pick_particle_;

  [[nodiscard]] real_t rand_uniform() { return uniform01_(rng_); }
  [[nodiscard]] real_t rand_proposal() { return proposal_(rng_); }
  [[nodiscard]] std::size_t rand_particle() { return pick_particle_(rng_); }

  [[nodiscard]] real_t acceptance_rate() const {
    if (proposed_ == 0U) {
      return 0.0_r;
    }
    return static_cast<real_t>(accepted_) / static_cast<real_t>(proposed_);
  }

  struct StepResult {
    bool accepted;
    std::size_t moved_particle;
    real_t old_x;
    real_t old_y;
    real_t old_z;
  };

  [[nodiscard]] std::vector<real_t> positions_snapshot() const;

public:
  struct MeasurementSummary {
    real_t mean_energy;
    std::optional<real_t> standard_error;
    real_t acceptance_rate;
  };

  explicit Simulation(Config cfg, std::unique_ptr<OutputWriter> output_writer = nullptr);
  MeasurementSummary run();

private:
  void initialize_positions();
  StepResult metropolis_step();
  void warmup();
  MeasurementSummary measure();
};
