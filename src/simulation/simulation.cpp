#include "simulation.hpp"

Simulation::Simulation(Config config) noexcept
: config_{std::move(config)}
, particles_{config_.numParticles}
, pbc_{config_.boxLength}
, jastrow_pade_{}
, slater_plane_wave_{particles_.numParticles(), pbc_.L()}
, wave_function_{jastrow_pade_, slater_plane_wave_}
, rng_{config.seed}
, proposal_{-config_.stepSize, config_.stepSize}
, pick_particle_{0, config_.numParticles - 1}
{ }

/// @brief Intializes Random Positions for each particle, making sure to not exceed box length
void Simulation::initialize_positions() {
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

  wave_function_.evaluateLogPsi(particles(), pbc());
  log_psi_current_ = *particles().logPsi();
  
  return;
}

/// @brief Randomly selects a particle and proposes a new position, then checks if the proposed move was valid then either accepts or rejects it
/// @return bool - If the proposed change was accepted 
bool Simulation::metropolis_step() {
    double* RESTRICT px{particles_.posX()};
    double* RESTRICT py{particles_.posY()};
    double* RESTRICT pz{particles_.posZ()};

    std::size_t i{pick_particle()(rng())};

    double old_x{px[i]};
    double old_y{py[i]};
    double old_z{pz[i]};

    px[i] += proposal()(rng());
    py[i] += proposal()(rng());
    pz[i] += proposal()(rng());

    wave_function_.evaluateLogPsi(particles(), pbc());

    double delta_log_psi{*particles().logPsi() - log_psi_current_};

    double log_u{log(uniform01()(rng()))};
    double minterm{std::min(0.0, 2 * delta_log_psi)};

    if (log_u < minterm) {
      log_psi_current_ = *particles().logPsi();
      return true;
    } else {
      px[i] = old_x;
      py[i] = old_y;
      pz[i] = old_z;
    }

    return false;
}

/// @brief Warmup the simulation by processing a small warmup sweep on particles
void Simulation::warmup() {
  const std::size_t warmupSteps{config_.warmupSteps};

  for (std::size_t i{}; i < warmupSteps; i++) {
    metropolis_step();
  }

  return;
}


void Simulation::measure() {
  const std::size_t measureSteps{config_.warmupSteps};

  for (std::size_t i{}; i < measureSteps; i++) {
    metropolis_step();
  }

  return;
}