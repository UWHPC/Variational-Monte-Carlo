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
  double log_psi_current_;

  std::mt19937_64 rng_;
  std::uniform_real_distribution<double> uniform01_{0.0, 1.0};
  std::uniform_real_distribution<double> proposal_;
  std::uniform_int_distribution<std::size_t> pick_particle_;

  [[nodiscard]] double rand_uniform() { return uniform01_(rng_); }
  [[nodiscard]] double rand_proposal() { return proposal_(rng_); }
  [[nodiscard]] std::size_t rand_particle() { return pick_particle_(rng_); }

  [[nodiscard]] double acceptance_rate() const {
    if (proposed_ == 0U) {
      return 0.0;
    }
    return static_cast<double>(accepted_) / static_cast<double>(proposed_);
  }

  struct StepResult {
    bool accepted;
    std::size_t moved_particle;
    double old_x;
    double old_y;
    double old_z;
  };

  [[nodiscard]] std::vector<double> positions_snapshot() const;

public:
  struct MeasurementSummary {
    double mean_energy;
    std::optional<double> standard_error;
    double acceptance_rate;
  };

  explicit Simulation(Config cfg, std::unique_ptr<OutputWriter> output_writer = nullptr);
  MeasurementSummary run();

private:
  void initialize_positions();
  StepResult metropolis_step();
  void warmup();
  MeasurementSummary measure();
};
