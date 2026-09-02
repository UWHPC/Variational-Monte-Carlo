#include <catch2/catch_test_macros.hpp>

#include "energy_tracking/energy_tracking.hpp"
#include "particles/particles.hpp"

#include <array>
#include <cmath>

TEST_CASE("Batched energy walkers match independent calculations", "[energy][walkers]") {
  constexpr auto num_particles{3uz};
  constexpr auto num_walkers{2uz};
  constexpr auto box_length{8.5_fp};
  constexpr auto tolerance{1e-10_fp};

  const auto set_state{
    [](Particles& target, std::size_t target_walker, std::size_t source_walker) {
      auto pos{target.pos(target_walker)};
      for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
        for (auto particle{0uz}; particle < num_particles; ++particle) {
          pos[axis][particle] =
            0.3_fp +
            0.7_fp * static_cast<fp_t>(source_walker) +
            0.4_fp * static_cast<fp_t>(axis) +
            1.1_fp * static_cast<fp_t>(particle);
        }
      }

      auto derivatives{target.derivatives(target_walker)};
      for (auto component{idx(Derivatives::GRAD_X)};
           component < idx(Derivatives::NUM);
           ++component) {
        for (auto particle{0uz}; particle < num_particles; ++particle) {
          derivatives[component][particle] =
            0.02_fp * static_cast<fp_t>(1uz + source_walker + component + particle);
        }
      }
    }
  };

  Particles particles{num_particles, num_walkers};
  EnergyTracker batched_tracker{box_length, particles};
  std::array<fp_t, num_walkers> batched_energies{};

  for (auto walker{0uz}; walker < num_walkers; ++walker) {
    set_state(particles, walker, walker);
    batched_tracker.initialize_structure_factors(particles.view(walker), walker);
    batched_tracker.initialize_reciprocal_energy(walker);
    batched_tracker.initialize_real_energy(particles.view(walker), walker);

    Particles independent_particles{num_particles};
    set_state(independent_particles, 0uz, walker);
    EnergyTracker independent_tracker{box_length, independent_particles};
    independent_tracker.initialize_structure_factors(independent_particles.view());
    independent_tracker.initialize_reciprocal_energy();
    independent_tracker.initialize_real_energy(independent_particles.view());

    batched_energies[walker] = batched_tracker.eval_total_energy(particles.view(walker), walker);
    const auto independent_energy{
      independent_tracker.eval_total_energy(independent_particles.view())
    };

    CAPTURE(walker, batched_energies[walker], independent_energy);
    REQUIRE(std::isfinite(batched_energies[walker]));
    REQUIRE(std::abs(batched_energies[walker] - independent_energy) <= tolerance);
  }

  REQUIRE(std::abs(batched_energies[0uz] - batched_energies[1uz]) > tolerance);
}
