#pragma once

#include "../support/checks.hpp"

#include <catch2/catch_test_macros.hpp>

#include "energy_tracking/energy_tracking.hpp"
#include "particles/particles.hpp"

#include <array>
#include <cmath>

namespace test {

inline constexpr auto energy_particle_count{3uz};
inline constexpr auto energy_walker_count{2uz};
inline constexpr auto energy_box_length{8.5_fp};

inline void set_energy_state(
  Particles& target,
  std::size_t target_walker,
  std::size_t source_walker
) {
  auto positions{target.pos(target_walker)};
  auto derivatives{target.derivatives(target_walker)};
  std::array<fp_t, energy_particle_count> host_values{};

  for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
    for (auto particle{0uz}; particle < energy_particle_count; ++particle) {
      host_values[particle] = 0.3_fp + 0.7_fp * scast<fp_t>(source_walker) + 0.4_fp * scast<fp_t>(axis) + 1.1_fp * scast<fp_t>(particle);
    }
    xpu::copy_n(positions[axis], host_values.data(), host_values.size());
  }

  for (auto component{idx(Derivatives::GRAD_X)};
       component < idx(Derivatives::NUM);
       ++component) {
    for (auto particle{0uz}; particle < energy_particle_count; ++particle) {
      host_values[particle] = 0.02_fp * scast<fp_t>(1uz + source_walker + component + particle);
    }
    xpu::copy_n(derivatives[component], host_values.data(), host_values.size());
  }
}

} // namespace test

TEST_CASE("Batched energy walkers match independent calculations", "[energy][walkers]") {
  Particles particles{test::energy_particle_count, test::energy_walker_count};
  EnergyTracker batched{test::energy_box_length, particles};
  std::array<fp_t, test::energy_walker_count> batched_energies{};

  for (auto walker{0uz}; walker < test::energy_walker_count; ++walker) {
    test::set_energy_state(particles, walker, walker);
    batched.initialize_structure_factors(particles.view(walker), walker);
    batched.initialize_reciprocal_energy(walker);
    batched.initialize_real_energy(particles.view(walker), walker);

    Particles independent_particles{test::energy_particle_count};
    test::set_energy_state(independent_particles, 0uz, walker);
    EnergyTracker independent{test::energy_box_length, independent_particles};
    independent.initialize_structure_factors(independent_particles.view());
    independent.initialize_reciprocal_energy();
    independent.initialize_real_energy(independent_particles.view());

    batched_energies[walker] = batched.eval_total_energy(particles.view(walker), walker);
    const auto independent_energy{independent.eval_total_energy(independent_particles.view())};

    CAPTURE(walker, batched_energies[walker], independent_energy);
    REQUIRE(std::isfinite(batched_energies[walker]));
    require_near(batched_energies[walker], independent_energy);
  }

  REQUIRE(std::abs(batched_energies[0uz] - batched_energies[1uz]) > test_tolerance);
}
