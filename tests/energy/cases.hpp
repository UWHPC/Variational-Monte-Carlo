#pragma once

#include "../support/checks.hpp"

#include <catch2/catch_test_macros.hpp>

#include "energy_tracking/energy_tracking.hpp"
#include "particles/particles.hpp"

#include <array>
#include <cmath>
#include <vector>

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
  }
  batched.initialize(particles, 2uz);

  for (auto walker{0uz}; walker < test::energy_walker_count; ++walker) {
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

TEST_CASE("Energy reinitialization and incremental moves agree with reconstruction", "[energy][walkers]") {
  constexpr auto walker_count{3uz};
  constexpr auto box_length{11.0_fp};

  for (const auto particle_count : {1uz, 37uz}) {
    Particles particles{particle_count, walker_count};
    EnergyTracker tracker{box_length, particles};
    std::vector<fp_t> values(particle_count);

    for (auto pass{0uz}; pass < 2uz; ++pass) {
      for (auto walker{0uz}; walker < walker_count; ++walker) {
        for (auto axis{0uz}; axis < idx(Axis::NUM); ++axis) {
          for (auto particle{0uz}; particle < particle_count; ++particle) {
            values[particle] = std::fmod(
              0.37_fp + scast<fp_t>(particle + 1uz) * (0.71_fp + 0.13_fp * scast<fp_t>(axis)) +
              0.019_fp * scast<fp_t>((pass + walker) * (particle + axis + 1uz)), box_length
            );
          }
          xpu::copy_n(particles.pos(walker)[axis], values.data(), particle_count);
        }
        for (auto derivative{0uz}; derivative < idx(Derivatives::NUM); ++derivative) {
          xpu::zero_n(particles.derivatives(walker)[derivative], particle_count);
        }
      }
      tracker.initialize(particles, 2uz);

      for (auto walker{0uz}; walker < walker_count; ++walker) {
        EnergyTracker rebuilt{box_length, particle_count};
        rebuilt.initialize_structure_factors(particles.view(walker));
        rebuilt.initialize_reciprocal_energy();
        rebuilt.initialize_real_energy(particles.view(walker));
        const auto expected{rebuilt.eval_total_energy(particles.view(walker))};
        const auto actual{tracker.eval_total_energy(particles.view(walker), walker)};
        const auto tolerance{100.0_fp * test_tolerance * xpu::max(1.0_fp, std::abs(expected))};
        CAPTURE(particle_count, pass, walker);
        REQUIRE(std::isfinite(actual));
        require_near(actual, expected, tolerance);

        constexpr auto moved{0uz};
        xpu::array<fp_t, idx(Axis::NUM)> old_pos{};
        for (auto axis{0uz}; axis < idx(Axis::NUM); ++axis) {
          xpu::copy_n(&old_pos[axis], particles.pos(walker)[axis] + moved, 1uz);
        }
        auto new_pos{old_pos};
        new_pos[idx(Axis::X)] += 0.07_fp;
        xpu::copy_n(particles.pos(walker)[idx(Axis::X)] + moved, &new_pos[idx(Axis::X)], 1uz);
        tracker.update_structure_factors(old_pos, new_pos, walker);
        tracker.update_real_energy(moved, old_pos, particles.view(walker), walker);

        rebuilt.initialize_structure_factors(particles.view(walker));
        rebuilt.initialize_reciprocal_energy();
        rebuilt.initialize_real_energy(particles.view(walker));
        require_near(
          tracker.eval_total_energy(particles.view(walker), walker),
          rebuilt.eval_total_energy(particles.view(walker)), tolerance
        );

        auto reciprocal{0.0_fp};
        xpu::copy_n(&reciprocal, tracker.view(walker).reciprocal_energy, 1uz);
        tracker.accept_move(0.25_fp, reciprocal, walker);
        require_near(
          tracker.eval_total_energy(particles.view(walker), walker),
          rebuilt.eval_total_energy(particles.view(walker)) + 0.25_fp, tolerance
        );
      }
    }
  }
}
