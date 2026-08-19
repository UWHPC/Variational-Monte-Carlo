#pragma once

#include "../blocking_analysis/blocking_analysis.hpp"
#include "../config/config.hpp"
#include "../energy_tracking/energy_tracking.hpp"
#include "../output_writer/output_writer.hpp"
#include "../particles/particles.hpp"
#include "../wavefunction/wavefunction.hpp"
#include <xpu/buffer.hpp>
#include <xpu/random.hpp>

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

namespace simulation {

struct StepResult {
  bool accepted;
  std::size_t moved_particle;
  xpu::array<fp_t, idx(Axis::NUM)> old_pos;
  xpu::array<fp_t, idx(Axis::NUM)> new_pos;
  fp_t log_psi_delta;
  fp_t real_energy_delta;
  fp_t reciprocal_energy;
};

} // namespace simulation

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
  fp_t log_psi_current_;

  std::array<std::vector<fp_t>, idx(Axis::NUM)> positions_;

  [[nodiscard]] fp_t acceptance_rate() const {
    if (proposed_ == 0uz) {
      return 0.0_fp;
    } else {
      return static_cast<fp_t>(accepted_) / static_cast<fp_t>(proposed_);
    }
  }

#if defined(XPU_CUDA)
  xpu::buffer<xpu::random::generator> walker_rng_;
  xpu::buffer<simulation::StepResult> step_result_;
#else
  xpu::random::generator walker_rng_;
#endif

  [[nodiscard]] const std::array<std::vector<fp_t>, idx(Axis::NUM)>& positions_snapshot();

public:
  struct MeasurementSummary {
    fp_t mean_energy;
    std::optional<fp_t> standard_error;
    fp_t acceptance_rate;
  };

  explicit Simulation(
    Config cfg,
    std::unique_ptr<OutputWriter> output_writer = nullptr,
    std::uint64_t walker_id = 0
  );
  MeasurementSummary run();

private:
  void initialize_positions();
  simulation::StepResult metropolis_step();
  void warmup();
  MeasurementSummary measure();
};
