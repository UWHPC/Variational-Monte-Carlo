#include "simulation.hpp"

Simulation::Simulation(Config config) noexcept
    : config_{std::move(config)}, particles_{config_.num_particles}, pbc_{config_.box_length}, rng_{config.seed},
      proposal_{-config_.step_size, config_.step_size}, pickParticle_{0, config_.num_particles - 1} {}

/// @brief Intializes Random Positions for each particle, making sure to not exceed box length
void Simulation::initializePositions() {
    std::size_t const N{particles_.num_particles_ptr()};

    // X, Y, Z Position Blocks of Memory
    double* RESTRICT p_x{particles_.pos_x_ptr()};
    double* RESTRICT p_y{particles_.pos_y_ptr()};
    double* RESTRICT p_z{particles_.pos_z_ptr()};

    std::mt19937_64 rng_local{rng()};

    const double LENGTH{pbc_.L_ptr()};

    for (std::size_t i{}; i < N; i++) {
        p_x[i] += uniform01_(rng_local) * LENGTH;
        p_y[i] += uniform01_(rng_local) * LENGTH;
        p_z[i] += uniform01_(rng_local) * LENGTH;
    }

    return;
}

bool Simulation::metropolisStep() {
    double* RESTRICT px{particles_.pos_x_ptr()};
    double* RESTRICT py{particles_.pos_y_ptr()};
    double* RESTRICT pz{particles_.pos_z_ptr()};

    std::size_t i{pickParticle()(rng())};

    double oldX{px[i]};
    double oldY{py[i]};
    double oldZ{pz[i]};

    px[i] += proposal()(rng());
    py[i] += proposal()(rng());
    pz[i] += proposal()(rng());

    pbc().wrap3(px[i], py[i], pz[i]);
}
