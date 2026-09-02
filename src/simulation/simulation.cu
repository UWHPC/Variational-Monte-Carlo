#include <xpu/xpu.hpp>
#include "simulation.hpp"
#include "simulation_kernels.hpp"

#include <cmath>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

Simulation::Simulation(
  Config config,
  std::unique_ptr<OutputWriter> output_writer,
  std::uint64_t walker_id
)
: config_{std::move(config)}
, particles_{config_.num_particles, config_.num_walkers}
, wave_function_{particles_, config_.box_length, config_.jastrow_a, config_.jastrow_b}
, blocking_analysis_{config_.block_size}
, energy_tracker_{config_.box_length, particles_}
, output_writer_{std::move(output_writer)}
, proposed_{}
, accepted_{}
, positions_{
    std::vector<fp_t>(particles_.count()),
    std::vector<fp_t>(particles_.count()),
    std::vector<fp_t>(particles_.count())
  }
, walker_rng_{config_.num_walkers}
, step_result_{config_.num_walkers}
, local_energies_{config_.num_walkers}
, walker_views_{config_.num_walkers}
, sweep_result_{1uz}
{
  std::vector<View> views{};
  views.reserve(config_.num_walkers);

  for (auto walker{0uz}; walker < config_.num_walkers; ++walker) {
    kernel::simulation::seed_generator(
      walker_rng_.data() + walker,
      config_.master_seed,
      walker_id * config_.num_walkers + walker
    );
    views.emplace_back(this->view(walker));
  }

  xpu::copy_n(walker_views_.data(), views.data(), views.size());
}

const std::array<std::vector<fp_t>, idx(Axis::NUM)>& Simulation::positions_snapshot() {
  const std::size_t N{particles_.count()};
  auto [p_x, p_y, p_z]{particles_.pos().pointers()};

  xpu::copy_n(positions_[idx(Axis::X)].data(), p_x, N);
  xpu::copy_n(positions_[idx(Axis::Y)].data(), p_y, N);
  xpu::copy_n(positions_[idx(Axis::Z)].data(), p_z, N);

  return positions_;
}

void Simulation::initialize_positions() {
  const fp_t length{config_.box_length};

  constexpr std::size_t MAX_INIT_ATTEMPTS{100};

  for (auto walker{0uz}; walker < particles_.walker_count(); ++walker) {
    for (std::size_t attempt = 0; attempt < MAX_INIT_ATTEMPTS; ++attempt) {
      kernel::simulation::initialize_positions(
        walker_rng_.data() + walker,
        length,
        particles_.pos(walker)
      );

      const auto log_psi{wave_function_.evaluate_log_psi(
        particles_.view(walker), walker
      )};
      if (std::isfinite(log_psi)) { break; }

      if (attempt == MAX_INIT_ATTEMPTS - 1) {
        throw std::runtime_error("Failed to find non-singular initial configuration");
      }
    }

    energy_tracker_.initialize_structure_factors(particles_.view(walker), walker);
    energy_tracker_.initialize_reciprocal_energy(walker);
    energy_tracker_.initialize_real_energy(particles_.view(walker), walker);
  }
}

Simulation::StepResult Simulation::metropolis_step() {
  const Simulation::StepResult result{
    kernel::simulation::metropolis_step(
      this->view(),
      config_.step_size
    )
  };

  return result;
}

Simulation::SweepResult Simulation::metropolis_sweep() {
  kernel::simulation::metropolis_sweep(
    walker_views_.data(),
    particles_.walker_count(),
    particles_.count(),
    config_.step_size,
    sweep_result_.data()
  );

  SweepResult result{};
  xpu::copy_n(&result, sweep_result_.data(), 1uz);
  return result;
}

void Simulation::measure_walkers() {
  kernel::simulation::measure_walkers(
    walker_views_.data(),
    particles_.walker_count()
  );
}

void Simulation::warmup() {
  auto& step_size{config_.step_size};

  constexpr auto target_rate{0.50_fp};
  constexpr auto initial_gain{0.25_fp};

  const auto max_step_size{0.5_fp * config_.box_length};
  step_size = xpu::min(step_size, max_step_size);

  for (auto sweep{0uz}; sweep < config_.warmup_sweeps; ++sweep) {
    const auto result{metropolis_sweep()};
    const auto acceptance_rate{result.acceptance_rate()};
    const auto gain{initial_gain * xpu::rsqrt(static_cast<fp_t>(sweep + 1uz))};

    step_size *= xpu::exp(gain * (acceptance_rate - target_rate));
    step_size = xpu::min(step_size, max_step_size);
  }
}

Simulation::MeasurementSummary Simulation::measure() {
  auto& particles{particles_};
  auto& blocking_analysis{blocking_analysis_};
  proposed_ = 0U;
  accepted_ = 0U;

  fp_t running_energy_sum{};
  fp_t final_mean_energy{};
  std::size_t sample_count{};
  std::optional<fp_t> final_standard_error{};
  std::vector<fp_t> local_energies(particles.walker_count());

  for (auto sweep{0uz}; sweep < config_.measure_sweeps; ++sweep) {
    const auto result{metropolis_sweep()};
    proposed_ += result.proposed;
    accepted_ += result.accepted;

    measure_walkers();
    xpu::copy_n(
      local_energies.data(),
      local_energies_.data(),
      local_energies.size()
    );

    auto sweep_energy_sum{0.0_fp};
    for (auto walker{0uz}; walker < particles.walker_count(); ++walker) {
      const auto local_energy{local_energies[walker]};
      running_energy_sum += local_energy;
      sweep_energy_sum += local_energy;
      ++sample_count;
      blocking_analysis.add(local_energy);
    }

    const auto local_energy{
      sweep_energy_sum / static_cast<fp_t>(particles.walker_count())
    };
    const fp_t running_mean{
      running_energy_sum / static_cast<fp_t>(sample_count)
    };
    final_mean_energy = running_mean;

    std::optional<fp_t> frame_standard_error{};
    if (blocking_analysis.ready()) {
      const auto [blocked_mean, standard_error]{blocking_analysis.mean_and_standard_error()};
      final_mean_energy = blocked_mean;
      final_standard_error = standard_error;
      frame_standard_error = standard_error;
    }

    if (output_writer_) {
      const auto& snapshot{positions_snapshot()};
      const std::size_t N{particles_.count()};

      std::vector<fp_t> flat_positions(N * 3U);
      for (std::size_t p{}; p < N; ++p) {
        flat_positions[p * 3U] = snapshot[0][p];
        flat_positions[p * 3U + 1U] = snapshot[1][p];
        flat_positions[p * 3U + 2U] = snapshot[2][p];
      }

      output_writer_->write_frame(FrameData{
        .step = sweep + 1uz,
        .accepted = accepted_,
        .proposed = proposed_,
        .acceptance_rate = acceptance_rate(),
        .local_energy = local_energy,
        .mean_energy = running_mean,
        .standard_error = frame_standard_error,
        .positions = std::move(flat_positions)
      });
    }

    if (config_.is_master_thread) {
      if ((sweep & 127uz) == 0uz || sweep + 1uz == config_.measure_sweeps) {
        std::cout << "\rProgress: "
                  << ((sweep + 1uz) * 100uz / config_.measure_sweeps)
                  << "%" << std::flush;
      }
    }
  }
  if (config_.is_master_thread) {
    std::cout << "\r" << std::string(20, ' ') << "\r";
  }

  return MeasurementSummary{
    .mean_energy = final_mean_energy,
    .standard_error = final_standard_error,
    .acceptance_rate = acceptance_rate()
  };
}

Simulation::MeasurementSummary Simulation::run() {
  initialize_positions();

  if (output_writer_) {
    output_writer_->write_init(InitData{
      .run_id = "vmc-seed-" + std::to_string(config_.master_seed),
      .num_particles = config_.num_particles,
      .box_length = config_.box_length,
      .warmup_steps = config_.warmup_steps,
      .measure_steps = config_.measure_steps,
      .step_size = config_.step_size,
      .seed = config_.master_seed,
      .block_size = config_.block_size
    });
  }

  warmup();
  const MeasurementSummary summary{measure()};

  if (output_writer_) {
    output_writer_->write_done(DoneData{
      .total_accepted = accepted_,
      .total_proposed = proposed_,
      .final_acceptance_rate = acceptance_rate(),
      .final_mean_energy = summary.mean_energy,
      .final_standard_error = summary.standard_error
    });
  }

  return summary;
}
