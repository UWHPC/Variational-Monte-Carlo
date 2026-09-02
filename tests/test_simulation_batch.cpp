#include <catch2/catch_test_macros.hpp>

#include "config/config.hpp"
#include "simulation/simulation.hpp"

#include <array>
#include <cmath>

TEST_CASE("Batched simulation matches independent walkers", "[simulation][walkers]") {
  Config config{};
  config.num_particles = 3uz;
  config.num_walkers = 2uz;
  config.warmup_sweeps = 0uz;
  config.measure_sweeps = 4uz;
  config.box_length = 6.0_fp;
  config.block_size = 100uz;
  config.master_seed = 20260901uz;

  Simulation batched{config};
  const auto batched_summary{batched.run()};

  config.num_walkers = 1uz;
  std::array<Simulation::MeasurementSummary, 2uz> independent{};
  for (auto walker{0uz}; walker < independent.size(); ++walker) {
    Simulation simulation{config, nullptr, walker};
    independent[walker] = simulation.run();
  }

  const auto expected_energy{
    0.5_fp * (independent[0uz].mean_energy + independent[1uz].mean_energy)
  };
  const auto expected_acceptance{
    0.5_fp * (independent[0uz].acceptance_rate + independent[1uz].acceptance_rate)
  };

  REQUIRE(std::isfinite(batched_summary.mean_energy));
  REQUIRE(std::abs(batched_summary.mean_energy - expected_energy) <= 1e-12_fp);
  REQUIRE(std::abs(batched_summary.acceptance_rate - expected_acceptance) <= 1e-12_fp);
}
