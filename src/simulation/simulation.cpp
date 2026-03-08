#include "simulation.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <utility>

#ifdef VMC_PROFILE_MODE
namespace {

class ScopedTimer {
public:
    explicit ScopedTimer(double& accumulator)
        : start_{std::chrono::steady_clock::now()}, accumulator_{accumulator} {}

    ~ScopedTimer() {
        const auto end{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed{end - start_};
        accumulator_ += elapsed.count();
    }

private:
    std::chrono::steady_clock::time_point start_;
    double& accumulator_;
};

} // namespace
#endif

Simulation::Simulation(Config config, std::unique_ptr<OutputWriter> output_writer)
    : config_{std::move(config)}, particles_{config_.num_particles}, pbc_{config_.box_length},
      wave_function_{config_.num_particles, config_.box_length}, blocking_analysis_{config_.block_size},
      energy_tracker_{config_.box_length, static_cast<double>(config_.num_particles)},
      output_writer_{std::move(output_writer)}, proposed_{}, accepted_{}, log_psi_current_{}, rng_{config_.seed},
      proposal_{-config_.step_size, config_.step_size}, pick_particle_{0, config_.num_particles - 1} {}

std::vector<double> Simulation::positions_snapshot() const {
    const std::size_t N{particles_.num_particles_get()};
    const double* RESTRICT p_x{particles_.pos_x_get()};
    const double* RESTRICT p_y{particles_.pos_y_get()};
    const double* RESTRICT p_z{particles_.pos_z_get()};

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

    const std::size_t N{particles_.num_particles_get()};
    const double length{pbc_.L_get()};

    // Generate Random Starting Positions
    for (std::size_t i{}; i < N; i++) {
        p_x[i] = rand_uniform_double() * length;
        p_y[i] = rand_uniform_double() * length;
        p_z[i] = rand_uniform_double() * length;
    }
    log_psi_set() = wave_function_get().evaluate_log_psi(particles_get(), pbc_get());
}

/// @brief Randomly selects a particle and proposes a new position, then checks if the proposed move was valid then
/// either accepts or rejects it
/// @return bool - Whether the proposed change was accepted
bool Simulation::metropolis_step() {
    // Position pointers:
    double* RESTRICT p_x{particles_.pos_x_get()};
    double* RESTRICT p_y{particles_.pos_y_get()};
    double* RESTRICT p_z{particles_.pos_z_get()};

    // Local random vars:
    const std::size_t rand_particle{rand_particle_get()};

    // Old positions:
    const double old_x{p_x[rand_particle]};
    const double old_y{p_y[rand_particle]};
    const double old_z{p_z[rand_particle]};

    // Add randomness:
    p_x[rand_particle] += rand_proposal_double();
    p_y[rand_particle] += rand_proposal_double();
    p_z[rand_particle] += rand_proposal_double();

    // Wrap to ensure particles dont drift outside [0,L):
    pbc_get().wrap3(p_x[rand_particle], p_y[rand_particle], p_z[rand_particle]);

    const double new_log_psi{wave_function_get().evaluate_log_psi(particles_get(), pbc_get())};
    const double delta_log_psi{new_log_psi - log_psi_get()};

    const double log_u{log(rand_uniform_double())};
    const double min_term{std::min(0.0, 2.0 * delta_log_psi)};

    bool accepted{log_u < min_term};

    if (accepted) {
        log_psi_set() = new_log_psi;
        return true;
    }

    p_x[rand_particle] = old_x;
    p_y[rand_particle] = old_y;
    p_z[rand_particle] = old_z;
    return false;
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
    const double gain{0.05};             // This limits making large changes to step size

    for (std::size_t i{}; i < warmup_steps; i++) {
        window_proposed++;
        bool accepted{};
#ifdef VMC_PROFILE_MODE
        {
            ScopedTimer timer{profile_stats_.metropolis_warmup_seconds};
            accepted = metropolis_step();
        }
        ++profile_stats_.warmup_metropolis_calls;
#else
        accepted = metropolis_step();
#endif
        if (accepted)
            ++window_accepted;

        if (window_proposed % warmup_batch_size == 0) {
            acceptance_rate_window = static_cast<double>(window_accepted) / static_cast<double>(window_proposed);
            step_size *= exp(gain * (acceptance_rate_window - acceptance_target));

            proposal().param(std::uniform_real_distribution<double>::param_type(-step_size, step_size));

            window_accepted = 0;
            window_proposed = 0;
        }
    }
}

Simulation::MeasurementSummary Simulation::measure() {
    const std::size_t measure_steps{config_.measure_steps};

    auto& wavefunction{wave_function_get()};
    auto& particles{particles_get()};
    auto& pbc{pbc_get()};
    auto& blocking_analysis{blocking_analysis_get()};
    auto& energy_tracker{energy_tracker_get()};
    proposed_ = 0U;
    accepted_ = 0U;

    double running_energy_sum{};
    double final_mean_energy{};
    std::optional<double> final_standard_error{};

    for (std::size_t i = 0; i < measure_steps; ++i) {
#ifdef VMC_PROFILE_MODE
        ++profile_stats_.measure_iterations;
#endif
        ++proposed_;
        bool accepted{};
#ifdef VMC_PROFILE_MODE
        {
            ScopedTimer timer{profile_stats_.metropolis_measure_seconds};
            accepted = metropolis_step();
        }
        ++profile_stats_.measure_metropolis_calls;
#else
        accepted = metropolis_step();
#endif
        if (accepted) {
            ++accepted_;
        }

#ifdef VMC_PROFILE_MODE
        {
            ScopedTimer timer{profile_stats_.evaluate_log_psi_seconds};
            wavefunction.evaluate_log_psi(particles, pbc);
        }
        {
            ScopedTimer timer{profile_stats_.evaluate_derivatives_seconds};
            wavefunction.evaluate_derivatives(particles, pbc);
        }
#else
        wavefunction.evaluate_log_psi(particles, pbc);
        wavefunction.evaluate_derivatives(particles, pbc);
#endif

#ifdef VMC_PROFILE_MODE
        double E_local{};
        {
            ScopedTimer timer{profile_stats_.energy_eval_seconds};
            E_local = energy_tracker.eval_total_energy(particles, pbc);
        }
#else
        const double E_local{energy_tracker.eval_total_energy(particles, pbc)};
#endif
        running_energy_sum += E_local;
#ifdef VMC_PROFILE_MODE
        {
            ScopedTimer timer{profile_stats_.blocking_seconds};
            blocking_analysis.add(E_local);
        }
#else
        blocking_analysis.add(E_local);
#endif

        const double running_mean{running_energy_sum / static_cast<double>(i + 1U)};
        final_mean_energy = running_mean;

        std::optional<double> frame_standard_error{};
        if (blocking_analysis.ready()) {
#ifdef VMC_PROFILE_MODE
            ScopedTimer timer{profile_stats_.blocking_seconds};
#endif
            const auto [blocked_mean, standard_error]{blocking_analysis.mean_and_standard_error()};
            final_mean_energy = blocked_mean;
            final_standard_error = standard_error;
            frame_standard_error = standard_error;
        }

        if (output_writer_) {
#ifdef VMC_PROFILE_MODE
            ScopedTimer timer{profile_stats_.output_seconds};
#endif
            output_writer_->write_frame(FrameData{.step = i + 1U,
                                                  .accepted = accepted_,
                                                  .proposed = proposed_,
                                                  .acceptance_rate = acceptance_rate(),
                                                  .local_energy = E_local,
                                                  .mean_energy = running_mean,
                                                  .standard_error = frame_standard_error,
                                                  .positions = positions_snapshot()});
        }
    }

    if (!output_writer_ && blocking_analysis.ready()) {
        // Pair of the mean and standard deviation:
        const auto [mean, stand_error]{blocking_analysis.mean_and_standard_error()};
        std::cout << "Energy: " << mean << " +/- " << stand_error << std::endl;
    }

    if (!output_writer_) {
        std::cout << "Acceptance Rate: " << acceptance_rate() << std::endl;
    }

    return MeasurementSummary{.mean_energy = final_mean_energy, .standard_error = final_standard_error};
}

void Simulation::run() {
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
#ifdef VMC_PROFILE_MODE
    print_profile_summary();
#endif
}

#ifdef VMC_PROFILE_MODE
void Simulation::print_profile_summary() const {
    const double timed_measure_total{
        profile_stats_.metropolis_measure_seconds + profile_stats_.evaluate_log_psi_seconds +
        profile_stats_.evaluate_derivatives_seconds + profile_stats_.energy_eval_seconds + profile_stats_.blocking_seconds +
        profile_stats_.output_seconds};

    const auto pct = [timed_measure_total](double seconds) {
        if (timed_measure_total <= 0.0) {
            return 0.0;
        }
        return (100.0 * seconds) / timed_measure_total;
    };

    const auto us_per = [](double seconds, std::uint64_t count) {
        if (count == 0U) {
            return 0.0;
        }
        return (seconds * 1'000'000.0) / static_cast<double>(count);
    };

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "[profile] warmup metropolis: " << profile_stats_.metropolis_warmup_seconds << " s ("
              << us_per(profile_stats_.metropolis_warmup_seconds, profile_stats_.warmup_metropolis_calls)
              << " us/call)\n";
    std::cout << "[profile] measure metropolis: " << profile_stats_.metropolis_measure_seconds << " s ("
              << pct(profile_stats_.metropolis_measure_seconds) << "%, "
              << us_per(profile_stats_.metropolis_measure_seconds, profile_stats_.measure_metropolis_calls)
              << " us/call)\n";
    std::cout << "[profile] evaluate_log_psi: " << profile_stats_.evaluate_log_psi_seconds << " s ("
              << pct(profile_stats_.evaluate_log_psi_seconds) << "%)\n";
    std::cout << "[profile] evaluate_derivatives: " << profile_stats_.evaluate_derivatives_seconds << " s ("
              << pct(profile_stats_.evaluate_derivatives_seconds) << "%)\n";
    std::cout << "[profile] energy eval: " << profile_stats_.energy_eval_seconds << " s ("
              << pct(profile_stats_.energy_eval_seconds) << "%)\n";
    std::cout << "[profile] blocking stats: " << profile_stats_.blocking_seconds << " s ("
              << pct(profile_stats_.blocking_seconds) << "%)\n";
    std::cout << "[profile] output write: " << profile_stats_.output_seconds << " s (" << pct(profile_stats_.output_seconds)
              << "%)\n";
    std::cout << "[profile] measured iterations: " << profile_stats_.measure_iterations << '\n';
    std::cout << "[profile] timed measure subtotal: " << timed_measure_total << " s\n";
}
#endif
