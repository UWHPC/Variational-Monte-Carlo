#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"
#include "wavefunction/wavefunction.hpp"

#include <algorithm>
#include <array>
#include <cmath>

TEST_CASE("Batched wavefunctions match independent walkers", "[wavefunction][walkers]") {
  constexpr auto num_particles{3uz};
  constexpr auto num_walkers{2uz};
  constexpr auto box_length{11.0_fp};
  constexpr std::array walker_positions{
    std::array{
      std::array{0.3_fp, 0.4_fp, 0.5_fp},
      std::array{1.7_fp, 0.2_fp, 2.1_fp},
      std::array{2.2_fp, 1.8_fp, 0.9_fp}
    },
    std::array{
      std::array{0.8_fp, 1.1_fp, 0.3_fp},
      std::array{2.7_fp, 0.4_fp, 1.9_fp},
      std::array{4.2_fp, 3.5_fp, 2.6_fp}
    }
  };

  const auto set_positions{
    [&](Particles& target, std::size_t target_walker, std::size_t source_walker) {
      auto pos{target.pos(target_walker)};
      for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
        for (auto particle{0uz}; particle < num_particles; ++particle) {
          pos[axis][particle] = walker_positions[source_walker][particle][axis];
        }
      }
    }
  };

  Particles particles{num_particles, num_walkers};
  for (auto walker{0uz}; walker < num_walkers; ++walker) {
    set_positions(particles, walker, walker);
  }

  WaveFunction batched{particles, box_length};
  std::array<fp_t, num_walkers> batched_log_psi{};

  for (auto walker{0uz}; walker < num_walkers; ++walker) {
    Particles independent_particles{num_particles};
    set_positions(independent_particles, 0uz, walker);
    WaveFunction independent{independent_particles, box_length};

    batched_log_psi[walker] = batched.evaluate_log_psi(particles.view(walker), walker);
    const auto independent_log_psi{independent.evaluate_log_psi(independent_particles.view())};

    batched.evaluate_derivatives(particles.view(walker), false, 0uz, {}, walker);
    independent.evaluate_derivatives(independent_particles.view());

    const auto batched_derivatives{particles.derivatives(walker)};
    const auto independent_derivatives{independent_particles.derivatives()};
    auto max_derivative_error{0.0_fp};
    for (auto derivative{idx(Derivatives::GRAD_X)}; derivative < idx(Derivatives::NUM); ++derivative) {
      for (auto particle{0uz}; particle < num_particles; ++particle) {
        max_derivative_error = std::max(
          max_derivative_error,
          std::abs(batched_derivatives[derivative][particle]
            - independent_derivatives[derivative][particle])
        );
      }
    }

    CAPTURE(walker, batched_log_psi[walker], independent_log_psi, max_derivative_error);
    REQUIRE(std::abs(batched_log_psi[walker] - independent_log_psi) <= 1e-12_fp);
    REQUIRE(max_derivative_error <= 1e-12_fp);
  }

  REQUIRE(std::abs(batched_log_psi[0uz] - batched_log_psi[1uz]) > 1e-12_fp);
}
