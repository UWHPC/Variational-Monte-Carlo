#include "simulation.hpp"

Simulation::Simulation(Config config) noexcept
: config_{std::move(config)}
, particles_{config.numParticles}
, pbc_{config.boxLength}
, jastrowPade_{}
, waveFunction_{}
{ }