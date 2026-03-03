#include "simulation.hpp"

Simulation::Simulation(Config config) noexcept
: config_{std::move(config)}
, particles_{config_.numParticles}
, pbc_{config_.boxLength}
, jastrowPade_{}
, waveFunction_{jastrowPade_}
{ }