#include "simulation.hpp"

Simulation::Simulation(Config config) noexcept
    : config_{std::move(config)}, particles_{config_.numParticles}, pbc_{config_.boxLength}, rng_{config.seed},
      proposal_{-config_.stepSize, config_.stepSize}, pickParticle_{0, config_.numParticles - 1} {}

/// @brief Intializes Random Positions for each particle, making sure to not exceed box length
void Simulation::initializePositions() {
    std::size_t const N{particles_.numParticles()};

    // X, Y, Z Position Blocks of Memory
    double* RESTRICT p_x{particles_.posX()};
    double* RESTRICT p_y{particles_.posY()};
    double* RESTRICT p_z{particles_.posZ()};

    std::mt19937_64 rng_local{rng()};

    const double length{pbc_.L()};

    for (std::size_t i{}; i < N; i++) {
        p_x[i] += uniform01_(rng_local) * length;
        p_y[i] += uniform01_(rng_local) * length;
        p_z[i] += uniform01_(rng_local) * length;
    }

    return;
}

bool Simulation::metropolisStep() {
    double* RESTRICT px{particles_.posX()};
    double* RESTRICT py{particles_.posY()};
    double* RESTRICT pz{particles_.posZ()};

    std::size_t i{pickParticle()(rng())};

    double oldX{px[i]};
    double oldY{py[i]};
    double oldZ{pz[i]};

    px[i] += proposal()(rng());
    py[i] += proposal()(rng());
    pz[i] += proposal()(rng());

    pbc().wrap3(px[i], py[i], pz[i]);
}
