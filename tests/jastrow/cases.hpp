#pragma once

#include "../support/checks.hpp"

#include "jastrow_pade/jastrow_pade.hpp"
#include "particles/particles.hpp"

#include <array>

namespace test {

inline void set_two_particle_positions(
  Particles& particles,
  const std::array<std::array<fp_t, 3uz>, 2uz>& values
) {
  auto positions{particles.pos()};
  std::array<fp_t, 2uz> host_positions{};

  for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
    for (auto particle{0uz}; particle < host_positions.size(); ++particle) {
      host_positions[particle] = values[particle][axis];
    }
    xpu::copy_n(positions[axis], host_positions.data(), host_positions.size());
  }
}

} // namespace test

TEST_CASE("Jastrow value uses minimum-image distance", "[jastrow]") {
  constexpr std::array positions{
    std::array{0.1_fp, 0.0_fp, 0.0_fp},
    std::array{9.9_fp, 0.0_fp, 0.0_fp}
  };
  Particles particles{2uz};
  test::set_two_particle_positions(particles, positions);
  const JastrowPade jastrow{10.0_fp, 0.25_fp, 1.0_fp};

  constexpr auto distance{0.2_fp};
  constexpr auto expected{0.25_fp * distance / (1.0_fp + distance)};
  require_near(jastrow.value(particles.view()), expected);
}

TEST_CASE("Jastrow value excludes coincident pairs", "[jastrow]") {
  constexpr std::array positions{
    std::array{1.0_fp, 2.0_fp, 3.0_fp},
    std::array{1.0_fp, 2.0_fp, 3.0_fp}
  };
  Particles particles{2uz};
  test::set_two_particle_positions(particles, positions);
  const JastrowPade jastrow{10.0_fp, 0.25_fp, 1.0_fp};

  require_near(jastrow.value(particles.view()), 0.0_fp);
}

TEST_CASE("Jastrow derivatives match the two-particle analytic result", "[jastrow]") {
  constexpr std::array positions{
    std::array{0.0_fp, 0.0_fp, 0.0_fp},
    std::array{1.0_fp, 0.0_fp, 0.0_fp}
  };
  Particles particles{2uz};
  test::set_two_particle_positions(particles, positions);
  const JastrowPade jastrow{100.0_fp, 0.25_fp, 1.0_fp};
  xpu::soa<fp_t, idx(Derivatives::NUM)> derivatives{2uz};
  xpu::zero_n(derivatives[0uz], derivatives.storage_size());

  jastrow.add_derivatives(
    particles.view(),
    derivatives.view<idx(Derivatives::NUM), 0uz>()
  );

  std::array<fp_t, 2uz> host_values{};
  constexpr std::array expected{
    std::array{-0.0625_fp, 0.0625_fp},
    std::array{0.0_fp, 0.0_fp},
    std::array{0.0_fp, 0.0_fp},
    std::array{0.0625_fp, 0.0625_fp}
  };

  for (auto derivative{idx(Derivatives::GRAD_X)};
       derivative < idx(Derivatives::NUM);
       ++derivative) {
    xpu::copy_n(host_values.data(), derivatives[derivative], host_values.size());
    for (auto particle{0uz}; particle < host_values.size(); ++particle) {
      require_near(host_values[particle], expected[derivative][particle]);
    }
  }
}
