#include "simulation.hpp"

Simulation::Simulation(Config config) noexcept
: config_{std::move(config)}
, particles_{config_.numParticles}
, pbc_{config_.boxLength}
, rng_{config.seed}
, proposal_{-config_.stepSize, config_.stepSize}
, pickParticle_{0, config_.numParticles - 1}
{ }

/// @brief Intializes Random Positions for each particle, making sure to not exceed box length
void Simulation::initializePositions() {

  // X, Y, Z Position Blocks of Memory
  double* x = particles_.posX();
  double* y = particles_.posY();
  double* z = particles_.posZ();

  std::size_t N = particles_.numParticles();

  double length{pbc_.L()};

  for (std::size_t i{}; i<N; i++) {
    x[i] += uniform01_(rng_) * length;
    y[i] += uniform01_(rng_) * length;
    z[i] += uniform01_(rng_) * length;
  }
  
  return;
}

bool Simulation::metropolisStep() {
  return false;
}