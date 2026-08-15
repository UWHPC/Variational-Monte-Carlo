#pragma once

#include "../blocking_analysis/blocking_analysis.hpp"
#include "../config/config.hpp"
#include "../energy_tracking/energy_tracking.hpp"
#include "../output_writer/output_writer.hpp"
#include "../particles/particles.hpp"
#include "../wavefunction/wavefunction.hpp"
#include "../utilities/random.hpp"
#include <xpu/soa.hpp>

#include <memory>
#include <optional>

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

  xpu::soa<real_t, 3> positions_;


  WalkerRNG walker_rng_;

  CUDA_CALLABLE [[nodiscard]] real_t rand_uniform() { return walker_rng_.rand_uniform(); }
  CUDA_CALLABLE [[nodiscard]] real_t rand_proposal() { return walker_rng_.rand_proposal(); }
  CUDA_CALLABLE [[nodiscard]] std::size_t rand_particle() { return walker_rng_.rand_particle(); }

  [[nodiscard]] real_t acceptance_rate() const {
    if (proposed_ == 0uz) {
      return 0.0_r;
    } else {
      return static_cast<real_t>(accepted_) / static_cast<real_t>(proposed_);
    }
  }

  struct StepResult {
    bool accepted;
    std::size_t moved_particle;
    xpu::array<real_t, idx(Axis::NUM)> old_pos;
  };

  [[nodiscard]] const xpu::soa<real_t, 3>& positions_snapshot();

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
  CUDA_CALLABLE StepResult metropolis_step();
  void warmup();
  MeasurementSummary measure();
};
