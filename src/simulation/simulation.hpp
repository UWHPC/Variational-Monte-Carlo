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

class Simulation {
public:
  struct RandomProposal {
    std::size_t particle;
    xpu::array<fp_t, idx(Axis::NUM)> displacement;
    fp_t acceptance;
  };

  struct StepResult {
    bool accepted;
    std::size_t moved_particle;
    xpu::array<fp_t, idx(Axis::NUM)> old_pos;
    xpu::array<fp_t, idx(Axis::NUM)> new_pos;
    fp_t log_psi_delta;
    fp_t real_energy_delta;
    fp_t reciprocal_energy;
  };

  struct MetropolisScratch {
    RandomProposal proposal;
    StepResult result;
    fp_t slater_ratio;
    fp_t jastrow_delta;
    fp_t real_energy_delta;
    fp_t reciprocal_sum;
  };

  struct SweepResult {
    std::size_t proposed;
    std::size_t accepted;

    [[nodiscard]] fp_t acceptance_rate() const noexcept {
      if (proposed == 0uz) { return 0.0_fp; }
      return static_cast<fp_t>(accepted) / static_cast<fp_t>(proposed);
    }
  };

  struct View {
    Particles::View particles{};
    WaveFunction::View wave_function{};
    EnergyTracker::View energy_tracker{};
    xpu::random::generator* generator{};
    StepResult* step_result{};
  };

  struct MeasurementSummary {
    fp_t mean_energy;
    std::optional<fp_t> standard_error;
    fp_t acceptance_rate;
  };

private:
  Config config_;

  Particles particles_;
  WaveFunction wave_function_;
  BlockingAnalysis blocking_analysis_;
  EnergyTracker energy_tracker_;
  std::unique_ptr<OutputWriter> output_writer_;

  std::size_t proposed_;
  std::size_t accepted_;

  std::array<std::vector<fp_t>, idx(Axis::NUM)> positions_;

  [[nodiscard]] fp_t acceptance_rate() const {
    if (proposed_ == 0uz) {
      return 0.0_fp;
    } else {
      return static_cast<fp_t>(accepted_) / static_cast<fp_t>(proposed_);
    }
  }

  xpu::buffer<xpu::random::generator> walker_rng_;
  xpu::buffer<StepResult> step_result_;
  xpu::buffer<View> walker_views_;
  xpu::buffer<SweepResult> sweep_result_;

  [[nodiscard]] const std::array<std::vector<fp_t>, idx(Axis::NUM)>& positions_snapshot();

public:
  explicit Simulation(
    Config cfg,
    std::unique_ptr<OutputWriter> output_writer = nullptr,
    std::uint64_t walker_id = 0
  );

  [[nodiscard]] CUDA_CALLABLE
  View view(std::size_t walker = 0uz) noexcept {
    return {
      particles_.view(walker),
      wave_function_.view(walker),
      energy_tracker_.view(walker),
      walker_rng_.data() + walker,
      step_result_.data() + walker
    };
  }

  MeasurementSummary run();
  StepResult metropolis_step();
  SweepResult metropolis_sweep();

private:
  void initialize_positions();
  void warmup();
  MeasurementSummary measure();
};
