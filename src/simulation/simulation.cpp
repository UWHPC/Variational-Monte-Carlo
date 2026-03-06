#include "simulation.hpp"

Simulation::Simulation(Config config) noexcept
    : config_{std::move(config)}, particles_{config_.num_particles}, pbc_{config_.box_length},
      wave_function_{config_.num_particles, config_.box_length}, blocking_analysis_{config_.block_size}, proposed_{},
      accepted_{}, log_psi_current_{}, rng_{config_.seed}, proposal_{-config_.step_size, config_.step_size},
      pick_particle_{0, config_.num_particles - 1} {}

double Simulation::kinetic_energy(const Particles& particles) const noexcept {
    const double* RESTRICT grad_x{particles.grad_log_psi_x_get()};
    const double* RESTRICT grad_y{particles.grad_log_psi_y_get()};
    const double* RESTRICT grad_z{particles.grad_log_psi_z_get()};
    const double* RESTRICT lap{particles.lap_log_psi_get()};

    // Kinetic
    double T_sum{};
    const std::size_t N{particles.num_particles_get()};

    for (std::size_t i = 0; i < N; i++) {
        // Computes ||Grad(logPsi)||^2
        const double grad_sq{grad_x[i] * grad_x[i] + grad_y[i] * grad_y[i] + grad_z[i] * grad_z[i]};

        // Accumulate Lapl(LogPsi) + ||Grad(LogPsi)||^2
        T_sum += (lap[i] + grad_sq);
    }

    return -0.5 * T_sum;
}

double Simulation::potential_energy(const Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept {
    const double* RESTRICT px{particles.pos_x_get()};
    const double* RESTRICT py{particles.pos_y_get()};
    const double* RESTRICT pz{particles.pos_z_get()};

    // Potential
    double V_sum{};
    const std::size_t N{particles.num_particles_get()};

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i + 1; j < N; ++j) {
            // To get around if else statement for branching,
            // ensures does not run if r_ij < 1.0e-12
            const double r_ij{pbc.distance(px[i], py[i], pz[i], px[j], py[j], pz[j])};
            const bool degenerate{r_ij < 1.0e-12};

            const double mask{degenerate ? 0.0 : 1.0};
            const double inv_r_ij{degenerate ? 1.0 : 1 / r_ij};

            // if degenerate == true, V_sum = 0.0 / 1.0 = 0
            // so the contribution is 0
            V_sum += mask * inv_r_ij;
        }
    }

    return V_sum;
}

double Simulation::eval_total_energy(const Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept {
    return kinetic_energy(particles) + potential_energy(particles, pbc);
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

    if (log_u < min_term) {
        log_psi_set() = new_log_psi;
        return true;
    } else {
        p_x[rand_particle] = old_x;
        p_y[rand_particle] = old_y;
        p_z[rand_particle] = old_z;
    }

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
        if (metropolis_step())
            window_accepted++;

        if (window_proposed % warmup_batch_size == 0) {
            acceptance_rate_window = static_cast<double>(window_accepted) / static_cast<double>(window_proposed);
            step_size *= exp(gain * (acceptance_rate_window - acceptance_target));

            proposal().param(std::uniform_real_distribution<double>::param_type(-step_size, step_size));

            window_accepted = 0;
            window_proposed = 0;
        }
    }
}

void Simulation::measure() {
    const std::size_t measure_steps{config_.measure_steps};

    auto& wavefunction{wave_function_get()};
    auto& particles{particles_get()};
    auto& pbc{pbc_get()};
    auto& blocking_analysis{blocking_analysis_get()};

    for (std::size_t i = 0; i < measure_steps; i++) {
        metropolis_step();
        wavefunction.evaluate_derivatives(particles, pbc);

        const double E_local{eval_total_energy(particles, pbc)};
        blocking_analysis.add(E_local);
    }

    if (blocking_analysis.ready()) {
        // Pair of the mean and standard deviation:
        const auto [mean, se]{blocking_analysis.mean_and_standard_error()};
        std::cout << "Energy: " << mean << " +/- " << se << '\n';
    }
}

void Simulation::run() {
    initialize_positions();
    warmup();
    measure();
}