#include "simulation.hpp"

Simulation::Simulation(Config config) noexcept
: config_{std::move(config)}
, particles_{config_.numParticles}
, pbc_{config_.boxLength}
, slaterPlaneWave_{0, 0}
, rng_{config.seed}
, proposal_{-config_.stepSize, config_.stepSize}
, pickParticle_{0, config_.numParticles - 1}
{ }

/// @brief Intializes Random Positions for each particle, making sure to not exceed box length
void Simulation::initializePositions() {
    std::size_t const N{particles_.numParticles()};

  // X, Y, Z Position Blocks of Memory
  double* RESTRICT x{particles_.posX()};
  double* RESTRICT y{particles_.posY()};
  double* RESTRICT z{particles_.posZ()};

  std::size_t N{particles_.numParticles()};
  const std::mt19937_64 rng_local{rng()};
  double length{pbc_.L()};

  // Generate Random Starting Positions
  for (std::size_t i{}; i<N; i++) {
    x[i] += uniform01_(rng_local) * length;
    y[i] += uniform01_(rng_local) * length;
    z[i] += uniform01_(rng_local) * length;
  }

  double slaterLogDet = slaterPlaneWave_.logAbsDet(particles_, pbc_);
  double jastrowVal = jastrowPade_.value(particles_, pbc_);
  
  logPsiCurrent_ = slaterLogDet + jastrowVal;
  
  return;
}

bool Simulation::metropolisStep() {
  return false;
}