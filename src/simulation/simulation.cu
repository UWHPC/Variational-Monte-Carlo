#include <xpu/xpu.hpp>
#include "simulation.hpp"
#include "simulation_kernels.hpp"

#include <cmath>
#include <iostream>
#include <limits>
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
, particles_{config_.num_particles}
, wave_function_{particles_, config_.box_length, config_.jastrow_a, config_.jastrow_b}
, blocking_analysis_{config_.block_size}
, energy_tracker_{config_.box_length, config_.num_particles}
, output_writer_{std::move(output_writer)}
, proposed_{}
, accepted_{}
, log_psi_current_{}
, positions_{
    std::vector<fp_t>(particles_.count()),
    std::vector<fp_t>(particles_.count()),
    std::vector<fp_t>(particles_.count())
  }
#if defined(XPU_CUDA)
, walker_rng_{1uz}
, step_result_{1uz}
#else
, walker_rng_{}
#endif
{
#if defined(XPU_CUDA)
  kernel::simulation::seed_generator(
    walker_rng_.data(),
    config_.master_seed,
    walker_id
  );
#else
  walker_rng_.seed(
    config_.master_seed,
    walker_id
  );
#endif
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

  for (std::size_t attempt = 0; attempt < MAX_INIT_ATTEMPTS; ++attempt) {
#if defined(XPU_CUDA)
    kernel::simulation::initialize_positions(
      walker_rng_.data(),
      length,
      particles_.pos()
    );
#else
    const std::size_t N{particles_.count()};
    auto [p_x, p_y, p_z]{particles_.pos().pointers()};

    for (std::size_t i = 0; i < N; i++) {
      p_x[i] = walker_rng_.uniform<fp_t>() * length;
      p_y[i] = walker_rng_.uniform<fp_t>() * length;
      p_z[i] = walker_rng_.uniform<fp_t>() * length;
    }
#endif

    log_psi_current_ = wave_function_.evaluate_log_psi(particles_);
    if (std::isfinite(log_psi_current_)) { break; }

    if (attempt == MAX_INIT_ATTEMPTS - 1) {
      throw std::runtime_error("Failed to find non-singular initial configuration");
    }
  }

  energy_tracker_.initialize_structure_factors(particles_);
  energy_tracker_.initialize_reciprocal_energy();
  energy_tracker_.initialize_real_energy(particles_);
}

simulation::StepResult Simulation::metropolis_step() {
#if defined(XPU_CUDA)
  auto& slater{wave_function_.slater_plane_wave()};
  const auto& jastrow{wave_function_.jastrow_pade()};

  const simulation::StepResult result{
    kernel::simulation::metropolis_step(
      walker_rng_.data(),
      config_.step_size,
      config_.box_length,
      jastrow.a(),
      jastrow.b(),
      particles_.pos(),
      slater.view(),
      energy_tracker_.view(),
      step_result_.data()
    )
  };
#else
  auto [p_x, p_y, p_z]{particles_.pos().pointers()};

  const auto random{
    stencil::simulation::generate_random_proposal(
      walker_rng_,
      particles_.count(),
      config_.step_size
    )
  };

  const std::size_t moved{random.particle};
  const fp_t L{config_.box_length};
  const fp_t inv_L{1.0_fp / L};

  const xpu::array<fp_t, idx(Axis::NUM)> old_pos{
    p_x[moved], p_y[moved], p_z[moved]
  };

  p_x[moved] += random.displacement[idx(Axis::X)];
  p_y[moved] += random.displacement[idx(Axis::Y)];
  p_z[moved] += random.displacement[idx(Axis::Z)];

  // Branchless wrapping for [0, L)
  p_x[moved] -= L * xpu::floor(p_x[moved] * inv_L);
  p_y[moved] -= L * xpu::floor(p_y[moved] * inv_L);
  p_z[moved] -= L * xpu::floor(p_z[moved] * inv_L);

  auto& slater{wave_function_.slater_plane_wave()};

  slater.update_trig_cache(moved, particles_);

  const fp_t* new_row{slater.build_row(moved)};
  const fp_t slater_ratio{slater.determinant_ratio(moved, new_row)};

  const fp_t delta_jastrow{
    wave_function_.jastrow_pade().delta_value(
      moved,
      old_pos,
      particles_.pos()
    )
  };
  const fp_t log_ratio_sq{2.0_fp * xpu::log(xpu::abs(slater_ratio)) + 2.0_fp * delta_jastrow};

  const fp_t u{xpu::max(random.acceptance, std::numeric_limits<fp_t>::min())};
  const fp_t log_u{xpu::log(u)};
  const fp_t min_term{xpu::min(0.0_fp, log_ratio_sq)};

  const bool accepted{log_u < min_term};

  const xpu::array<fp_t, idx(Axis::NUM)> new_pos{
    p_x[moved], p_y[moved], p_z[moved]
  };

  const simulation::StepResult result{
    accepted,
    moved,
    old_pos,
    new_pos,
    xpu::log(xpu::abs(slater_ratio)) + delta_jastrow,
    0.0_fp,
    0.0_fp
  };

  if (accepted) {
    slater.accept_move(moved, new_row, slater_ratio);
  } else {
    slater.restore_trig_row(moved);
    p_x[moved] = old_pos[idx(Axis::X)];
    p_y[moved] = old_pos[idx(Axis::Y)];
    p_z[moved] = old_pos[idx(Axis::Z)];
  }
#endif

  if (result.accepted) {
    log_psi_current_ += result.log_psi_delta;

#if defined(XPU_CUDA)
    energy_tracker_.accept_move(
      result.real_energy_delta,
      result.reciprocal_energy
    );
#else
    energy_tracker_.update_structure_factors(
      result.old_pos, result.new_pos
    );
    energy_tracker_.update_real_energy(
      result.moved_particle,
      result.old_pos,
      particles_
    );
#endif
  }

  return result;
}

void Simulation::warmup() {
  auto& step_size{config_.step_size};

  constexpr auto target_rate{0.50_fp};
  constexpr auto initial_gain{0.25_fp};

  const auto proposals_per_sweep{particles_.count()};

  const auto max_step_size{0.5_fp * config_.box_length};
  step_size = xpu::min(step_size, max_step_size);

  for (auto sweep{0uz}; sweep < config_.warmup_sweeps; ++sweep) {
    auto accepted_proposals{0uz};

    for (auto proposal{0uz}; proposal < proposals_per_sweep; ++proposal) {
      if (metropolis_step().accepted) { ++accepted_proposals; }
    }

    const auto acceptance_rate{static_cast<fp_t>(accepted_proposals) / static_cast<fp_t>(proposals_per_sweep)};
    const auto gain{initial_gain * xpu::rsqrt(static_cast<fp_t>(sweep + 1uz))};

    step_size *= xpu::exp(gain * (acceptance_rate - target_rate));
    step_size = xpu::min(step_size, max_step_size);
  }
}

Simulation::MeasurementSummary Simulation::measure() {
  const std::size_t measure_steps{config_.measure_steps};

  auto& wavefunction{wave_function_};
  auto& particles{particles_};
  auto& blocking_analysis{blocking_analysis_};
  auto& energy_tracker{energy_tracker_};
  proposed_ = 0U;
  accepted_ = 0U;

  fp_t running_energy_sum{};
  fp_t final_mean_energy{};
  std::optional<fp_t> final_standard_error{};

  for (std::size_t i = 0; i < measure_steps; ++i) {
    ++proposed_;
    const simulation::StepResult result{metropolis_step()};
    if (result.accepted) {
      ++accepted_;
    }

    wavefunction.evaluate_derivatives(
      particles, result.accepted, result.moved_particle,
      result.old_pos
    );

    const fp_t E_local{energy_tracker.eval_total_energy(particles)};
    running_energy_sum += E_local;
    blocking_analysis.add(E_local);

    const fp_t running_mean{running_energy_sum / static_cast<fp_t>(i + 1U)};
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
        .step = i + 1U,
        .accepted = accepted_,
        .proposed = proposed_,
        .acceptance_rate = acceptance_rate(),
        .local_energy = E_local,
        .mean_energy = running_mean,
        .standard_error = frame_standard_error,
        .positions = std::move(flat_positions)
      });
    }

    if (config_.is_master_thread) {
      if ((i & 127) == 0 || i == measure_steps) {
        std::cout << "\rProgress: " << (i * 100 / measure_steps) << "%" << std::flush;
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
