#include "simulation.hpp"

Simulation::Simulation(Config config) noexcept
    : config_{std::move(config)}, particles_{config_.num_particles}, pbc_{config_.box_length},
      wave_function_{config_.num_particles, config_.box_length}, proposed_{}, accepted_{}, log_psi_current_{},
      rng_{config_.seed}, proposal_{-config_.step_size, config_.step_size},
      pick_particle_{0, config_.num_particles - 1} {}

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

    wave_function().evaluate_log_psi(particles(), pbc());
    log_psi_current() = particles().log_psi_get()[0];

    return;
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
    const double old_X{p_x[rand_particle]};
    const double old_y{p_y[rand_particle]};
    const double old_z{p_z[rand_particle]};

    // Add randomness:
    p_x[rand_particle] += rand_proposal_double();
    p_y[rand_particle] += rand_proposal_double();
    p_z[rand_particle] += rand_proposal_double();

    // Wrap to ensure particles dont drift outside [0,L):
    pbc().wrap3(p_x[rand_particle], p_y[rand_particle], p_z[rand_particle]);

    wave_function().evaluate_log_psi(particles(), pbc());

    // Log psi for 1st element:
    const double new_log_psi{particles().log_psi_get()[0]};
    const double delta_log_psi{new_log_psi - log_psi_current()};

    const double log_u{log(rand_uniform_double())};
    const double min_term{std::min(0.0, 2.0 * delta_log_psi)};

    if (log_u < min_term) {
        log_psi_current_ = new_log_psi;
        return true;
    }
    p_x[rand_particle] = old_X;
    p_y[rand_particle] = old_y;
    p_z[rand_particle] = old_z;

    return false;
}

/// @brief Warmup the simulation by processing a small warmup sweep on particles
void Simulation::warmup() {
    const std::size_t warmup_steps{config_.warmup_steps};

    for (std::size_t i{}; i < warmup_steps; i++) {
        metropolis_step();
    }

    return;
}

void Simulation::measure() {
    const std::size_t measure_steps{config_.measure_steps};

    for (std::size_t i{}; i < measure_steps; i++) {
        metropolis_step();
    }

    return;
}

void Simulation::run() { std::cout << "hi" << std::endl; }