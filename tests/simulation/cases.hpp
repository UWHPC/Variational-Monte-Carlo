#pragma once

#include "../support/checks.hpp"
#include "../support/simulation.hpp"

#include "simulation/simulation.hpp"

#include <array>
#include <cmath>
#include <memory>

TEST_CASE("Batched simulation matches independent walkers", "[simulation][walkers]") {
  Config config{make_config(3uz, 6.0_fp, 0uz, 4uz, 0.6_fp, 20260901uz, 3uz, 2uz)};

  Simulation batched{config};
  const auto batched_summary{batched.run()};

  config.num_walkers = 1uz;
  config.measure_steps = config.num_particles * config.measure_sweeps;
  std::array<Simulation::MeasurementSummary, 2uz> independent{};
  for (auto walker{0uz}; walker < independent.size(); ++walker) {
    Simulation simulation{config, nullptr, walker};
    independent[walker] = simulation.run();
  }

  const auto expected_energy{0.5_fp * (independent[0uz].mean_energy + independent[1uz].mean_energy)};
  const auto expected_acceptance{0.5_fp * (independent[0uz].acceptance_rate + independent[1uz].acceptance_rate)};

  REQUIRE(std::isfinite(batched_summary.mean_energy));
  REQUIRE_FALSE(batched_summary.standard_error.has_value());
  require_near(batched_summary.mean_energy, expected_energy);
  require_near(batched_summary.acceptance_rate, expected_acceptance);
}

TEST_CASE("Simulation reports consistent frames and final summary", "[simulation][output]") {
  const Config config{make_config(3uz, 6.0_fp, 0uz, 4uz, 0.35_fp, 12345uz, 2uz)};
  auto writer{std::make_unique<RecordingOutputWriter>()};
  auto* const output{writer.get()};

  Simulation simulation{config, std::move(writer)};
  const auto summary{simulation.run()};

  REQUIRE(output->init.has_value());
  REQUIRE(output->done.has_value());
  REQUIRE(output->frames.size() == config.measure_sweeps);

  for (auto frame_index{0uz}; frame_index < output->frames.size(); ++frame_index) {
    const auto& frame{output->frames[frame_index]};
    const auto expected_step{frame_index + 1uz};
    REQUIRE(frame.step == expected_step);
    REQUIRE(frame.proposed == expected_step * config.num_particles);
    REQUIRE(frame.accepted <= frame.proposed);
    require_near(frame.acceptance_rate, scast<fp_t>(frame.accepted) / scast<fp_t>(frame.proposed));
    REQUIRE(frame.positions.size() == 3uz * config.num_particles);
  }

  const auto& done{*output->done};
  REQUIRE(done.total_proposed == output->frames.back().proposed);
  REQUIRE(done.total_accepted == output->frames.back().accepted);
  require_near(done.final_mean_energy, summary.mean_energy);
  require_near(done.final_acceptance_rate, summary.acceptance_rate);
  REQUIRE(done.final_standard_error.has_value());
}

TEST_CASE("Zero-size proposals preserve energy and are always accepted", "[simulation]") {
  const Config config{make_config(3uz, 6.0_fp, 0uz, 6uz, 0.0_fp, 314159uz, 3uz)};
  auto writer{std::make_unique<RecordingOutputWriter>()};
  auto* const output{writer.get()};

  Simulation simulation{config, std::move(writer)};
  simulation.run();

  const auto first_energy{output->frames.front().local_energy};
  for (const auto& frame : output->frames) {
    REQUIRE(frame.accepted == frame.proposed);
    require_near(frame.local_energy, first_energy);
  }
}
