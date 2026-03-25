#include "simulation.hpp"

#include <cmath>
#include <limits>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <stdexcept>

Simulation::Simulation(Config config, std::unique_ptr<OutputWriter> output_writer)
    : config_{std::move(config)}, particles_{config_.num_particles},
      wave_function_{particles_, config_.box_length}, blocking_analysis_{config_.block_size},
      energy_tracker_{config_.box_length, static_cast<double>(config_.num_particles)},
      output_writer_{std::move(output_writer)}, proposed_{}, accepted_{}, log_psi_current_{},
      rng_{config_.seed}, proposal_{-config_.step_size, config_.step_size},
      pick_particle_{0, config_.num_particles - 1} {}

std::vector<double> Simulation::positions_snapshot() const {
    const std::size_t N{particles_.num_particles_get()};
    const double* RESTRICT p_x{particles_.pos_x_get()};
    const double* RESTRICT p_y{particles_.pos_y_get()};
    const double* RESTRICT p_z{particles_.pos_z_get()};

    ASSUME_ALIGNED(p_x, SIMD_BYTES);
    ASSUME_ALIGNED(p_y, SIMD_BYTES);
    ASSUME_ALIGNED(p_z, SIMD_BYTES);

    std::vector<double> positions{};
    positions.reserve(N * 3U);

    for (std::size_t i{}; i < N; ++i) {
        positions.push_back(p_x[i]);
        positions.push_back(p_y[i]);
        positions.push_back(p_z[i]);
    }
    return positions;
}

/// @brief Intializes Random Positions for each particle, making sure to not exceed box length
void Simulation::initialize_positions() {
    // X, Y, Z Position Blocks of Memory
    double* RESTRICT p_x{particles_.pos_x_get()};
    double* RESTRICT p_y{particles_.pos_y_get()};
    double* RESTRICT p_z{particles_.pos_z_get()};

    ASSUME_ALIGNED(p_x, SIMD_BYTES);
    ASSUME_ALIGNED(p_y, SIMD_BYTES);
    ASSUME_ALIGNED(p_z, SIMD_BYTES);

    const std::size_t N{particles_.num_particles_get()};
    const double length{config_.box_length};

    constexpr std::size_t MAX_INIT_ATTEMPTS{100};

    for (std::size_t attempt = 0; attempt < MAX_INIT_ATTEMPTS; ++attempt) {
        for (std::size_t i = 0; i < N; i++) {
            p_x[i] = rand_uniform_double() * length;
            p_y[i] = rand_uniform_double() * length;
            p_z[i] = rand_uniform_double() * length;
        }

        log_psi_set() = wave_function_get().evaluate_log_psi(particles_get());
        if (std::isfinite(log_psi_get())) break;

        if (attempt == MAX_INIT_ATTEMPTS - 1) {
            throw std::runtime_error("Failed to find non-singular initial configuration");
        }
    }

    energy_tracker_get().initialize_structure_factors(particles_get());
    energy_tracker_get().initialize_reciprocal_energy();
    energy_tracker_get().initialize_real_energy(particles_get());
}

/// @brief Randomly selects a particle and proposes a new position, then checks if the proposed move
/// was valid then either accepts or rejects it
Simulation::StepResult Simulation::metropolis_step() {
    // Position pointers:
    double* RESTRICT p_x{particles_.pos_x_get()};
    double* RESTRICT p_y{particles_.pos_y_get()};
    double* RESTRICT p_z{particles_.pos_z_get()};

    ASSUME_ALIGNED(p_x, SIMD_BYTES);
    ASSUME_ALIGNED(p_y, SIMD_BYTES);
    ASSUME_ALIGNED(p_z, SIMD_BYTES);

    // Local random vars:
    const std::size_t rand_particle{rand_particle_get()};
    const double L{config_.box_length};
    const double inv_L{1.0 / L};

    // Old positions:
    const double old_x{p_x[rand_particle]};
    const double old_y{p_y[rand_particle]};
    const double old_z{p_z[rand_particle]};

    // Add randomness:
    p_x[rand_particle] += rand_proposal_double();
    p_y[rand_particle] += rand_proposal_double();
    p_z[rand_particle] += rand_proposal_double();

    // Branchless wrapping for [0, L)
    p_x[rand_particle] -= L * std::floor(p_x[rand_particle] * inv_L);
    p_y[rand_particle] -= L * std::floor(p_y[rand_particle] * inv_L);
    p_z[rand_particle] -= L * std::floor(p_z[rand_particle] * inv_L);

    // Build new Slater row for moved particle and compute determinant ratio - O(N):
    auto& slater{wave_function_get().slater_plane_wave_get()};

    // Save trig in case need to revert:
    slater.save_trig_row(rand_particle);

    // Update it to proceed:
    slater.update_trig_cache(rand_particle, particles_get());

    const double* new_row{slater.build_row(rand_particle)};
    const double slater_ratio{slater.determinant_ratio(rand_particle, new_row)};

    // Compute new Jastrow value:
    const double delta_jastrow{wave_function_get().jastrow_pade_get().delta_value(
        particles_get(), rand_particle, old_x, old_y, old_z)};
    const double log_ratio_sq{2.0 * std::log(std::abs(slater_ratio)) + 2.0 * delta_jastrow};

    const double u{std::max(rand_uniform_double(), std::numeric_limits<double>::min())};
    const double log_u{std::log(u)};
    const double min_term{std::min(0.0, log_ratio_sq)};

    const bool accepted{log_u < min_term};

    if (accepted) {
        // Update log_psi:
        log_psi_set() += std::log(std::abs(slater_ratio)) + delta_jastrow;

        // Sherman-Morrison inverse update:
        slater.accept_move(rand_particle, new_row, slater_ratio);

        // Update Ewald structure factors, and potential energies:
        energy_tracker_get().update_structure_factors(old_x, old_y, old_z, p_x[rand_particle],
                                                      p_y[rand_particle], p_z[rand_particle]);
        energy_tracker_get().update_real_energy(rand_particle, old_x, old_y, old_z,
                                                particles_get());

        return StepResult{true, rand_particle, old_x, old_y, old_z};
    }

    // Restore positions and trig if rejected:
    slater.restore_trig_row(rand_particle);
    p_x[rand_particle] = old_x;
    p_y[rand_particle] = old_y;
    p_z[rand_particle] = old_z;

    return StepResult{false, rand_particle, old_x, old_y, old_z};
}

/// @brief Warmup the simulation by processing a small warmup sweep on particles
void Simulation::warmup() {
    double& step_size{config_.step_size};

    const std::size_t warmup_steps{config_.warmup_steps};
    const std::size_t warmup_batch_size{particles_get().num_particles_get()};

    std::size_t window_proposed{};
    std::size_t window_accepted{};

    double acceptance_rate_window{};
    const double acceptance_target{0.5}; // Currently targeting a 50% acceptance rate
    const double gain{0.25};             // This limits making large changes to step size

    for (std::size_t i{}; i < warmup_steps; i++) {
        window_proposed++;
        const bool accepted{metropolis_step().accepted};
        if (accepted)
            ++window_accepted;

        if (window_proposed % warmup_batch_size == 0) {
            acceptance_rate_window =
                static_cast<double>(window_accepted) / static_cast<double>(window_proposed);
            step_size *= exp(gain * (acceptance_rate_window - acceptance_target));

            proposal().param(
                std::uniform_real_distribution<double>::param_type(-step_size, step_size));

            window_accepted = 0;
            window_proposed = 0;
        }
    }
}

Simulation::MeasurementSummary Simulation::measure() {
    const std::size_t measure_steps{config_.measure_steps};

    auto& wavefunction{wave_function_get()};
    auto& particles{particles_get()};
    auto& blocking_analysis{blocking_analysis_get()};
    auto& energy_tracker{energy_tracker_get()};
    proposed_ = 0U;
    accepted_ = 0U;

    double running_energy_sum{};
    double final_mean_energy{};
    std::optional<double> final_standard_error{};

    for (std::size_t i = 0; i < measure_steps; ++i) {
        ++proposed_;
        const StepResult result{metropolis_step()};
        if (result.accepted) {
            ++accepted_;
        }

        wavefunction.evaluate_derivatives(particles, result.accepted, result.moved_particle,
                                          result.old_x, result.old_y, result.old_z);

        const double E_local{energy_tracker.eval_total_energy(particles)};
        running_energy_sum += E_local;
        blocking_analysis.add(E_local);

        const double running_mean{running_energy_sum / static_cast<double>(i + 1U)};
        final_mean_energy = running_mean;

        std::optional<double> frame_standard_error{};
        if (blocking_analysis.ready()) {
            const auto [blocked_mean, standard_error]{blocking_analysis.mean_and_standard_error()};
            final_mean_energy = blocked_mean;
            final_standard_error = standard_error;
            frame_standard_error = standard_error;
        }

        if (output_writer_) {
            output_writer_->write_frame(FrameData{.step = i + 1U,
                                                  .accepted = accepted_,
                                                  .proposed = proposed_,
                                                  .acceptance_rate = acceptance_rate(),
                                                  .local_energy = E_local,
                                                  .mean_energy = running_mean,
                                                  .standard_error = frame_standard_error,
                                                  .positions = positions_snapshot()});
        }

        if (config_.is_master_thread) {
            if ((i & 127) == 0 || i == measure_steps) {
               std::cout << "\rProgress: " << (i * 100 / measure_steps) << "%" << std::flush;
            }
        }
        
    }
    if (config_.is_master_thread) {
        // Clear the progress bar:
        std::cout << "\r" << std::string(20, ' ') << "\r";
    }

    return MeasurementSummary{.mean_energy = final_mean_energy,
                              .standard_error = final_standard_error,
                              .acceptance_rate = acceptance_rate()};
}

Simulation::MeasurementSummary Simulation::run() {
    initialize_positions();

    if (output_writer_) {
        output_writer_->write_init(InitData{.run_id = "vmc-seed-" + std::to_string(config_.seed),
                                            .num_particles = config_.num_particles,
                                            .box_length = config_.box_length,
                                            .warmup_steps = config_.warmup_steps,
                                            .measure_steps = config_.measure_steps,
                                            .step_size = config_.step_size,
                                            .seed = config_.seed,
                                            .block_size = config_.block_size});
    }

    warmup();
    const MeasurementSummary summary{measure()};

    if (output_writer_) {
        output_writer_->write_done(DoneData{.total_accepted = accepted_,
                                            .total_proposed = proposed_,
                                            .final_acceptance_rate = acceptance_rate(),
                                            .final_mean_energy = summary.mean_energy,
                                            .final_standard_error = summary.standard_error});
    }

    return summary;
}
