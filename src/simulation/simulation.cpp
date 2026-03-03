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
  // X, Y, Z Position Blocks of Memory
  double* RESTRICT p_x{particles_.posX()};
  double* RESTRICT p_y{particles_.posY()};
  double* RESTRICT p_z{particles_.posZ()};

  const std::size_t N{particles_.numParticles()};
  std::mt19937_64 rng_local{rng()};
  double length{pbc_.L()};

  // Generate Random Starting Positions
  for (std::size_t i{}; i<N; i++) {
    p_x[i] = uniform01_(rng_local) * length;
    p_y[i] = uniform01_(rng_local) * length;
    p_z[i] = uniform01_(rng_local) * length;
  }

  double slaterLogDet = slaterPlaneWave_.logAbsDet(particles_, pbc_);
  double jastrowVal = jastrowPade_.value(particles_, pbc_);
  
  logPsiCurrent_ = slaterLogDet + jastrowVal;
  
  return;
}

bool Simulation::metropolisStep() {
    double* RESTRICT px{particles_.posX()};
    double* RESTRICT py{particles_.posY()};
    double* RESTRICT pz{particles_.posZ()};

    std::size_t i{pickParticle_(rng())};

    double oldX{px[i]};
    double oldY{py[i]};
    double oldZ{pz[i]};

    px[i] += proposal()(rng());
    py[i] += proposal()(rng());
    pz[i] += proposal()(rng());

    pbc().wrap3(*px, *py, *pz);
}

/// @brief Warmup the simulation by processing a small warmup sweep on particles
void Simulation::warmup() {
  const std::size_t warmupSteps{config_.warmupSteps};

  for (std::size_t i{}; i < warmupSteps; i++) {
    metropolisStep();
  }

  return;
}


void Simulation::measure() {
  const std::size_t measureSteps{config_.warmupSteps};

  for (std::size_t i{}; i < measureSteps; i++) {
    metropolisStep();
  }

  slaterPlaneWave_.addDerivatives(particles_, pbc(), 0, 0, 0, 0);
  jastrowPade_.addDerivatives(particles_, pbc(), 0, 0, 0, 0);

  return;
}