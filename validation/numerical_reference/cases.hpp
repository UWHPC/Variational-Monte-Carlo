#pragma once

#include "../support/physics.hpp"
#include "energy_tracking/energy_tracking.hpp"

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <vector>

#if defined(FP_64)

inline constexpr auto ewald_tolerance{8e-6_fp};
inline constexpr auto slater_update_tolerance{1e-8_fp};

#else

inline constexpr auto ewald_tolerance{8e-3_fp};
inline constexpr auto slater_update_tolerance{2e-2_fp};

#endif

namespace validation {

template <std::size_t particle_count>
[[nodiscard]] fp_t direct_ewald_energy(
  const std::array<std::array<fp_t, idx(Axis::NUM)>, particle_count>& positions,
  fp_t box_length
) {
  constexpr auto reciprocal_tolerance{1e-6_fp};
  constexpr auto alpha_coefficient{6.0_fp};
  const auto alpha{alpha_coefficient / box_length};
  fp_t real_energy{};

  for (auto first{0uz}; first < particle_count; ++first) {
    for (auto second{first + 1uz}; second < particle_count; ++second) {
      fp_t distance_squared{};
      for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
        auto displacement{positions[first][axis] - positions[second][axis]};
        displacement -= box_length * std::nearbyint(displacement / box_length);
        distance_squared += displacement * displacement;
      }
      const auto distance{std::sqrt(distance_squared)};
      real_energy += std::erfc(alpha * distance) / distance;
    }
  }

  const auto two_pi_over_length{
    2.0_fp * std::numbers::pi_v<fp_t> / box_length
  };
  const auto maximum_wave_number_squared{
    -4.0_fp * alpha * alpha * std::log(reciprocal_tolerance)
  };
  const auto maximum_index{
    scast<int>(std::ceil(
      std::sqrt(maximum_wave_number_squared) / two_pi_over_length
    )) + 1
  };
  fp_t reciprocal_energy{};

  for (auto x{-maximum_index}; x <= maximum_index; ++x) {
    for (auto y{-maximum_index}; y <= maximum_index; ++y) {
      for (auto z{-maximum_index}; z <= maximum_index; ++z) {
        if (x < 0 || (x == 0 && y < 0) || (x == 0 && y == 0 && z <= 0)) {
          continue;
        }
        const std::array wave_vector{
          two_pi_over_length * scast<fp_t>(x),
          two_pi_over_length * scast<fp_t>(y),
          two_pi_over_length * scast<fp_t>(z)
        };
        const auto wave_number_squared{
          wave_vector[0uz] * wave_vector[0uz] +
          wave_vector[1uz] * wave_vector[1uz] +
          wave_vector[2uz] * wave_vector[2uz]
        };
        if (wave_number_squared > maximum_wave_number_squared) {
          continue;
        }

        fp_t real_sum{};
        fp_t imaginary_sum{};
        for (const auto& position : positions) {
          const auto phase{
            wave_vector[0uz] * position[0uz] +
            wave_vector[1uz] * position[1uz] +
            wave_vector[2uz] * position[2uz]
          };
          real_sum += std::cos(phase);
          imaginary_sum += std::sin(phase);
        }
        const auto weight{
          8.0_fp * std::numbers::pi_v<fp_t> * std::numbers::pi_v<fp_t> *
          std::exp(-wave_number_squared / (4.0_fp * alpha * alpha)) /
          wave_number_squared
        };
        reciprocal_energy += weight * (
          real_sum * real_sum + imaginary_sum * imaginary_sum
        );
      }
    }
  }

  reciprocal_energy /= 2.0_fp * std::numbers::pi_v<fp_t> *
    box_length * box_length * box_length;
  const auto count{scast<fp_t>(particle_count)};
  const auto self{-alpha * count / std::sqrt(std::numbers::pi_v<fp_t>)};
  const auto background{
    -std::numbers::pi_v<fp_t> * count * count /
    (2.0_fp * alpha * alpha * box_length * box_length * box_length)
  };
  return real_energy + reciprocal_energy + self + background;
}

} // namespace validation

TEST_CASE(
  "Energy tracker agrees with an independent truncated Ewald sum",
  "[validation][numerical-reference][ewald]"
) {
  constexpr auto particle_count{3uz};
  constexpr auto box_length{8.25_fp};
  constexpr std::array positions{
    std::array{0.55_fp, 1.20_fp, 2.10_fp},
    std::array{3.15_fp, 2.45_fp, 1.05_fp},
    std::array{5.60_fp, 0.85_fp, 4.45_fp}
  };

  Particles particles{particle_count};
  validation::set_positions(particles, positions);
  validation::zero_derivatives(particles);
  EnergyTracker tracker{box_length, particles};
  tracker.initialize_structure_factors(particles.view());
  tracker.initialize_reciprocal_energy();
  tracker.initialize_real_energy(particles.view());

  const auto actual{tracker.eval_total_energy(particles.view())};
  const auto expected{validation::direct_ewald_energy(positions, box_length)};
  CAPTURE(actual, expected);
  REQUIRE(std::abs(actual - expected) <= ewald_tolerance);
}

TEST_CASE(
  "Repeated Slater updates agree with fresh determinant reconstruction",
  "[validation][numerical-reference][slater]"
) {
  constexpr auto particle_count{3uz};
  constexpr auto box_length{10.0_fp};
  constexpr auto update_count{32uz};
  constexpr std::array initial_positions{
    std::array{0.7_fp, 1.4_fp, 2.1_fp},
    std::array{2.8_fp, 0.9_fp, 1.6_fp},
    std::array{4.3_fp, 3.2_fp, 0.5_fp}
  };

  auto positions{initial_positions};
  Particles particles{particle_count};
  validation::set_positions(particles, positions);
  SlaterPlaneWave maintained{particles, box_length};
  REQUIRE(std::isfinite(maintained.log_abs_det(particles.view())));

  for (auto update{0uz}; update < update_count; ++update) {
    const auto moved{update % particle_count};
    positions[moved][0uz] = validation::wrap_coordinate(
      positions[moved][0uz] + 0.003_fp * scast<fp_t>(update % 5uz + 1uz),
      box_length
    );
    positions[moved][1uz] = validation::wrap_coordinate(
      positions[moved][1uz] - 0.0025_fp * scast<fp_t>(update % 7uz + 1uz),
      box_length
    );
    positions[moved][2uz] = validation::wrap_coordinate(
      positions[moved][2uz] + 0.002_fp * scast<fp_t>(update % 3uz + 1uz),
      box_length
    );
    validation::set_positions(particles, positions);
    maintained.update_trig_cache(moved, particles.view());
    const auto* const new_row{maintained.build_row(moved)};
    const auto ratio{maintained.determinant_ratio(moved, new_row)};
    REQUIRE(std::isfinite(ratio));
    maintained.accept_move(moved, new_row, ratio);

    SlaterPlaneWave rebuilt{particles, box_length};
    REQUIRE(std::isfinite(rebuilt.log_abs_det(particles.view())));
    const auto maintained_matrix{validation::determinant_matrix(maintained)};
    const auto rebuilt_matrix{validation::determinant_matrix(rebuilt)};
    fp_t maximum_difference{};
    for (auto element{0uz}; element < maintained_matrix.size(); ++element) {
      maximum_difference = std::max(
        maximum_difference,
        std::abs(maintained_matrix[element] - rebuilt_matrix[element])
      );
    }

    CAPTURE(update, moved, ratio, maximum_difference);
    REQUIRE(maximum_difference <= slater_update_tolerance);
  }
}
