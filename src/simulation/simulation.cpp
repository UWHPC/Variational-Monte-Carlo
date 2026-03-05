#include "simulation.hpp"

Simulation::Simulation(Config config) noexcept
    : config_{std::move(config)}, particles_{config_.num_particles}, pbc_{config_.box_length}, jastrow_pade_{},
      slater_plane_wave_{particles_.num_particles_ptr(), pbc_.L_ptr()},
      wave_function_{jastrow_pade_, slater_plane_wave_}, rng_{config.seed},
      proposal_{-config_.step_size, config_.step_size}, pick_particle_{0, config_.num_particles - 1} {}

/// @brief Intializes Random Positions for each particle, making sure to not exceed box length
void Simulation::initialize_positions() {
    // X, Y, Z Position Blocks of Memory
    double* RESTRICT p_x{particles_.pos_x_ptr()};
    double* RESTRICT p_y{particles_.pos_y_ptr()};
    double* RESTRICT p_z{particles_.pos_z_ptr()};

    const std::size_t N{particles_.num_particles_ptr()};
    std::mt19937_64 rng_local{rng()};
    double length{pbc_.L_ptr()};

    // Generate Random Starting Positions
    for (std::size_t i{}; i < N; i++) {
        p_x[i] = uniform01_(rng_local) * length;
        p_y[i] = uniform01_(rng_local) * length;
        p_z[i] = uniform01_(rng_local) * length;
    }

    wave_function_.evaluate_log_psi(particles(), pbc());
    log_psi_current_ = *particles().log_psi_ptr();

    return;
}

/// @brief Randomly selects a particle and proposes a new position, then checks if the proposed move was valid then
/// either accepts or rejects it
/// @return bool - Whether the proposed change was accepted
bool Simulation::metropolis_step() {
    // Position pointers:
    double* RESTRICT p_x{particles_.pos_x_ptr()};
    double* RESTRICT p_y{particles_.pos_y_ptr()};
    double* RESTRICT p_z{particles_.pos_z_ptr()};

    // Local random vars:
    const double RAND_PROPOSAL{proposal()(rng())};
    const std::size_t RAND_PARTICLE{pick_particle()(rng())};

    // Log psi for 1st element:
    const double LOG_PSI{particles().log_psi_ptr()[0]};

    // Old positions:
    const double OLD_X{p_x[RAND_PARTICLE]};
    const double OLD_Y{p_y[RAND_PARTICLE]};
    const double OLD_Z{p_z[RAND_PARTICLE]};

    p_x[RAND_PARTICLE] += RAND_PROPOSAL;
    p_y[RAND_PARTICLE] += RAND_PROPOSAL;
    p_z[RAND_PARTICLE] += RAND_PROPOSAL;

    wave_function_.evaluate_log_psi(particles(), pbc());

    const double DELTA_LOG_PSI{LOG_PSI - log_psi_current_};

    const double LOG_U{log(uniform01()(rng()))};
    const double MIN_TERM{std::min(0.0, 2.0 * DELTA_LOG_PSI)};

    if (LOG_U < MIN_TERM) {
        log_psi_current_ = LOG_PSI;
        return true;
    } else {
        p_x[RAND_PARTICLE] = OLD_X;
        p_y[RAND_PARTICLE] = OLD_Y;
        p_z[RAND_PARTICLE] = OLD_Z;
    }

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