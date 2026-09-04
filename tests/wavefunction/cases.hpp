#pragma once

#include "../support/checks.hpp"

#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"
#include "wavefunction/wavefunction.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace test {

inline constexpr auto wavefunction_particle_count{3uz};
inline constexpr auto wavefunction_walker_count{2uz};
inline constexpr auto wavefunction_box_length{11.0_fp};
inline constexpr std::array wavefunction_positions{
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

inline void set_wavefunction_positions(
  Particles& target,
  std::size_t target_walker,
  std::size_t source_walker
) {
  auto positions{target.pos(target_walker)};
  std::array<fp_t, wavefunction_particle_count> host_positions{};

  for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
    for (auto particle{0uz}; particle < wavefunction_particle_count; ++particle) {
      host_positions[particle] = wavefunction_positions[source_walker][particle][axis];
    }
    xpu::copy_n(positions[axis], host_positions.data(), host_positions.size());
  }
}

} // namespace test

TEST_CASE("Batched wavefunctions match independent walkers", "[wavefunction][walkers]") {
  Particles particles{test::wavefunction_particle_count, test::wavefunction_walker_count};
  for (auto walker{0uz}; walker < test::wavefunction_walker_count; ++walker) {
    test::set_wavefunction_positions(particles, walker, walker);
  }

  WaveFunction batched{particles, test::wavefunction_box_length};
  std::array<fp_t, test::wavefunction_walker_count> batched_log_psi{};

  for (auto walker{0uz}; walker < test::wavefunction_walker_count; ++walker) {
    Particles independent_particles{test::wavefunction_particle_count};
    test::set_wavefunction_positions(independent_particles, 0uz, walker);
    WaveFunction independent{independent_particles, test::wavefunction_box_length};

    batched_log_psi[walker] = batched.evaluate_log_psi(particles.view(walker), walker);
    const auto independent_log_psi{independent.evaluate_log_psi(independent_particles.view())};

    batched.evaluate_derivatives(particles.view(walker), walker);
    independent.evaluate_derivatives(independent_particles.view());

    auto max_derivative_error{0.0_fp};
    std::array<fp_t, test::wavefunction_particle_count> batched_derivatives{};
    std::array<fp_t, test::wavefunction_particle_count> independent_derivatives{};
    for (auto derivative{idx(Derivatives::GRAD_X)};
         derivative < idx(Derivatives::NUM);
         ++derivative) {
      xpu::copy_n(
        batched_derivatives.data(),
        particles.derivatives(walker)[derivative],
        batched_derivatives.size()
      );
      xpu::copy_n(
        independent_derivatives.data(),
        independent_particles.derivatives()[derivative],
        independent_derivatives.size()
      );

      for (auto particle{0uz}; particle < test::wavefunction_particle_count; ++particle) {
        max_derivative_error = std::max(
          max_derivative_error,
          std::abs(batched_derivatives[particle] - independent_derivatives[particle])
        );
      }
    }

    CAPTURE(walker, batched_log_psi[walker], independent_log_psi, max_derivative_error);
    require_near(batched_log_psi[walker], independent_log_psi);
    REQUIRE(max_derivative_error <= test_tolerance);
  }

  REQUIRE(
    std::abs(batched_log_psi[0uz] - batched_log_psi[1uz]) >
    test_tolerance
  );
}
