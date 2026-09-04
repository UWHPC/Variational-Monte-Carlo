#pragma once

#include "../support/benchmarks.hpp"
#include "config/config.hpp"
#include "simulation/simulation.hpp"

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

TEST_CASE(
  "Finite-cell VMC energies remain near the PZ81 thermodynamic-limit reference",
  "[validation][literature]"
) {
  for (const auto& benchmark : validation::literature_benchmarks) {
    Config config{};
    config.num_particles = validation::literature_particle_count;
    config.num_walkers = validation::literature_walker_count;
    config.warmup_sweeps = validation::literature_warmup_sweeps;
    config.measure_sweeps = validation::literature_measure_sweeps;
    config.box_length = validation::box_length_from_r_s(
      benchmark.r_s,
      validation::literature_particle_count
    );
    config.block_size = validation::literature_block_size;
    config.master_seed = benchmark.seed;
    config.step_size = benchmark.step_size;
    config.is_master_thread = false;

    Simulation simulation{config};
    const auto summary{simulation.run()};
    const auto energy_per_particle{
      summary.mean_energy /
      scast<fp_t>(validation::literature_particle_count)
    };
    const auto reference_energy{validation::pz81_total_energy(benchmark.r_s)};
    const auto difference{
      std::abs(energy_per_particle - reference_energy)
    };

    CAPTURE(
      validation::literature_citation,
      validation::literature_provenance,
      benchmark.r_s,
      reference_energy,
      energy_per_particle,
      difference,
      benchmark.comparison_tolerance,
      summary.standard_error,
      summary.acceptance_rate
    );
    REQUIRE(std::isfinite(energy_per_particle));
    REQUIRE(difference <= benchmark.comparison_tolerance);
  }
}
