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
    const double LENGTH{pbc_.L_get()};

    // Generate Random Starting Positions
    for (std::size_t i{}; i < N; i++) {
        p_x[i] = rand_uniform_double() * LENGTH;
        p_y[i] = rand_uniform_double() * LENGTH;
        p_z[i] = rand_uniform_double() * LENGTH;
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
    const std::size_t RAND_PARTICLE{rand_particle()};

    // Old positions:
    const double OLD_X{p_x[RAND_PARTICLE]};
    const double OLD_Y{p_y[RAND_PARTICLE]};
    const double OLD_Z{p_z[RAND_PARTICLE]};

    // Add randomness:
    p_x[RAND_PARTICLE] += rand_proposal_double();
    p_y[RAND_PARTICLE] += rand_proposal_double();
    p_z[RAND_PARTICLE] += rand_proposal_double();

    // Wrap to ensure particles dont drift outside [0,L):
    pbc().wrap3(p_x[RAND_PARTICLE], p_y[RAND_PARTICLE], p_z[RAND_PARTICLE]);

    wave_function().evaluate_log_psi(particles(), pbc());

    // Log psi for 1st element:
    const double NEW_LOG_PSI{particles().log_psi_get()[0]};
    const double DELTA_LOG_PSI{NEW_LOG_PSI - log_psi_current()};

    const double LOG_U{log(rand_uniform_double())};
    const double MIN_TERM{std::min(0.0, 2.0 * DELTA_LOG_PSI)};

    if (LOG_U < MIN_TERM) {
        log_psi_current_ = NEW_LOG_PSI;
        return true;
    } else {
        p_x[RAND_PARTICLE] = OLD_X;
        p_y[RAND_PARTICLE] = OLD_Y;
        p_z[RAND_PARTICLE] = OLD_Z;

        // revert original log_psi in particles buffer
        *particles().log_psi_get() = log_psi_current(); 
    }

    return false;
}

/// @brief Warmup the simulation by processing a small warmup sweep on particles
void Simulation::warmup() {
    const std::size_t WARMUP_STEPS{config_.warmup_steps};
    const std::size_t WARMUP_BATCH_SIZE{particles().num_particles_get()}; 

    std::size_t window_proposed{};
    std::size_t window_accepted{};

    double acceptance_rate_window{};
    const double acceptance_target{0.5}; // Currently targeting a 50% acceptance rate
    const double gain{0.05}; // This limits making large changes to step size

    for (std::size_t i{}; i < WARMUP_STEPS; i++) {
        window_proposed++;
        if (metropolis_step()) window_accepted++;

        if (window_proposed % WARMUP_BATCH_SIZE == 0) {
            acceptance_rate_window = static_cast<double>(window_accepted)/static_cast<double>(window_proposed);
            config_.step_size *= exp(gain*(acceptance_rate_window - acceptance_target));

            proposal().param(std::uniform_real_distribution<double>::param_type(-config_.step_size, config_.step_size));

            window_accepted = 0;
            window_proposed = 0;
        }
    }

    return;
}

void Simulation::measure() {
    const std::size_t MEASURE_STEPS{config_.measure_steps};

    for (std::size_t i{}; i < MEASURE_STEPS; i++) {
        metropolis_step();
    }

    return;
}
