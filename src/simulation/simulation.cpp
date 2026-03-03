#include "simulation.hpp"

Simulation::Simulation(Config config) noexcept
: config_{std::move(config)}
, particles_{config_.numParticles}
, pbc_{config_.boxLength}
, jastrowPade_{}
, waveFunction_{jastrowPade_}
, rng_{config.seed}
, proposal_{-config_.stepSize, config_.stepSize}
, pickParticle_{0, config_.numParticles - 1}
{ }



bool Simulation::metropolisStep() {

}