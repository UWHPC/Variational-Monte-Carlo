#pragma once

#include "../support/physics.hpp"
#include "jastrow_pade/jastrow_pade.hpp"
#include "wavefunction/wavefunction.hpp"

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <numbers>

#if defined(FP_64)

inline constexpr auto analytical_tolerance{1e-8_fp};
inline constexpr auto periodic_tolerance{1e-9_fp};

#else

inline constexpr auto analytical_tolerance{2e-3_fp};
inline constexpr auto periodic_tolerance{2e-3_fp};

#endif

TEST_CASE(
  "Closed-shell free gas has the analytical kinetic energy",
  "[validation][analytical][free-gas]"
) {
  constexpr auto particle_count{7uz};
  constexpr auto box_length{6.0_fp};
  constexpr auto two{2.0_fp};
  constexpr auto three{3.0_fp};
  constexpr std::array positions{
    std::array{0.7_fp, 1.1_fp, 2.2_fp},
    std::array{1.4_fp, 2.8_fp, 0.6_fp},
    std::array{2.3_fp, 0.9_fp, 4.7_fp},
    std::array{3.1_fp, 4.2_fp, 1.5_fp},
    std::array{4.6_fp, 2.0_fp, 3.4_fp},
    std::array{5.2_fp, 3.7_fp, 4.9_fp},
    std::array{0.4_fp, 5.1_fp, 3.0_fp}
  };

  Particles particles{particle_count};
  validation::set_positions(particles, positions);
  WaveFunction wave_function{particles, box_length, 0.0_fp, 1.0_fp};
  scast<void>(wave_function.evaluate_log_psi(particles.view()));
  wave_function.evaluate_derivatives(particles.view());

  const auto wave_number{two * std::numbers::pi_v<fp_t> / box_length};
  const auto expected{three * wave_number * wave_number};
  const auto orbital_sum{validation::exact_free_gas_kinetic_energy(
    wave_function.slater_plane_wave()
  )};
  const auto local{validation::local_kinetic_energy(particles)};

  CAPTURE(expected, orbital_sum, local);
  REQUIRE(std::abs(orbital_sum - expected) <= analytical_tolerance);
  REQUIRE(std::abs(local - expected) <= analytical_tolerance);
}

TEST_CASE(
  "Closed-shell free gas local energy is position independent",
  "[validation][analytical][free-gas]"
) {
  constexpr auto particle_count{7uz};
  constexpr auto box_length{5.5_fp};
  constexpr std::array first_positions{
    std::array{0.7_fp, 1.1_fp, 2.2_fp},
    std::array{1.4_fp, 2.8_fp, 0.6_fp},
    std::array{2.3_fp, 0.9_fp, 4.7_fp},
    std::array{3.1_fp, 4.2_fp, 1.5_fp},
    std::array{4.6_fp, 2.0_fp, 3.4_fp},
    std::array{5.2_fp, 3.7_fp, 4.9_fp},
    std::array{0.4_fp, 5.1_fp, 3.0_fp}
  };
  constexpr std::array second_positions{
    std::array{0.2_fp, 3.4_fp, 1.7_fp},
    std::array{1.8_fp, 0.5_fp, 4.1_fp},
    std::array{2.7_fp, 4.8_fp, 0.9_fp},
    std::array{3.6_fp, 1.4_fp, 2.8_fp},
    std::array{4.9_fp, 2.3_fp, 3.7_fp},
    std::array{5.3_fp, 3.9_fp, 1.2_fp},
    std::array{0.8_fp, 5.0_fp, 4.6_fp}
  };

  Particles particles{particle_count};
  WaveFunction wave_function{particles, box_length, 0.0_fp, 1.0_fp};

  validation::set_positions(particles, first_positions);
  scast<void>(wave_function.evaluate_log_psi(particles.view()));
  wave_function.evaluate_derivatives(particles.view());
  const auto first_energy{validation::local_kinetic_energy(particles)};

  validation::set_positions(particles, second_positions);
  scast<void>(wave_function.evaluate_log_psi(particles.view()));
  wave_function.evaluate_derivatives(particles.view());
  const auto second_energy{validation::local_kinetic_energy(particles)};

  CAPTURE(first_energy, second_energy);
  REQUIRE(std::abs(first_energy - second_energy) <= analytical_tolerance);
}

TEST_CASE(
  "Fully polarized Pade Jastrow satisfies the same-spin cusp",
  "[validation][analytical][jastrow]"
) {
  constexpr auto particle_count{2uz};
  constexpr auto box_length{50.0_fp};
  constexpr auto separation{1e-3_fp};
  constexpr auto cusp{0.25_fp};
  constexpr std::array positions{
    std::array{0.0_fp, 0.0_fp, 0.0_fp},
    std::array{separation, 0.0_fp, 0.0_fp}
  };

  Particles particles{particle_count};
  validation::set_positions(particles, positions);
  validation::zero_derivatives(particles);
  const JastrowPade jastrow{box_length, cusp, 1.0_fp};
  jastrow.add_derivatives(particles.view(), particles.derivatives());

  std::array<fp_t, particle_count> gradient_x{};
  xpu::copy_n(
    gradient_x.data(),
    particles.derivatives()[idx(Derivatives::GRAD_X)],
    gradient_x.size()
  );
  const auto expected{cusp / ((1.0_fp + separation) * (1.0_fp + separation))};

  CAPTURE(gradient_x, expected);
  REQUIRE(std::abs(gradient_x[1uz] - expected) <= analytical_tolerance);
  REQUIRE(std::abs(gradient_x[0uz] + expected) <= analytical_tolerance);
}

TEST_CASE(
  "Wavefunction is invariant under periodic box translations",
  "[validation][analytical][periodic]"
) {
  constexpr auto particle_count{3uz};
  constexpr auto box_length{10.0_fp};
  constexpr std::array baseline_positions{
    std::array{1.1_fp, 2.2_fp, 3.3_fp},
    std::array{4.4_fp, 5.5_fp, 6.6_fp},
    std::array{7.7_fp, 1.8_fp, 2.9_fp}
  };
  constexpr std::array translated_positions{
    std::array{1.1_fp, 2.2_fp, 3.3_fp},
    std::array{14.4_fp, -4.5_fp, 16.6_fp},
    std::array{7.7_fp, 1.8_fp, 2.9_fp}
  };

  Particles baseline_particles{particle_count};
  Particles translated_particles{particle_count};
  validation::set_positions(baseline_particles, baseline_positions);
  validation::set_positions(translated_particles, translated_positions);
  WaveFunction baseline{baseline_particles, box_length};
  WaveFunction translated{translated_particles, box_length};

  const auto baseline_log_psi{baseline.evaluate_log_psi(baseline_particles.view())};
  const auto translated_log_psi{
    translated.evaluate_log_psi(translated_particles.view())
  };
  baseline.evaluate_derivatives(baseline_particles.view());
  translated.evaluate_derivatives(translated_particles.view());

  CAPTURE(baseline_log_psi, translated_log_psi);
  REQUIRE(std::abs(baseline_log_psi - translated_log_psi) <= periodic_tolerance);

  std::array<fp_t, particle_count> baseline_values{};
  std::array<fp_t, particle_count> translated_values{};
  for (auto component{idx(Derivatives::GRAD_X)};
       component < idx(Derivatives::NUM);
       ++component) {
    xpu::copy_n(
      baseline_values.data(),
      baseline_particles.derivatives()[component],
      baseline_values.size()
    );
    xpu::copy_n(
      translated_values.data(),
      translated_particles.derivatives()[component],
      translated_values.size()
    );
    for (auto particle{0uz}; particle < particle_count; ++particle) {
      CAPTURE(component, particle, baseline_values[particle], translated_values[particle]);
      REQUIRE(
        std::abs(baseline_values[particle] - translated_values[particle]) <=
        periodic_tolerance
      );
    }
  }
}
