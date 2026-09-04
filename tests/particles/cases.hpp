#pragma once

#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"

#include <array>

TEST_CASE("Particle walkers keep independent physical state", "[particles][walkers]") {
  constexpr auto particle_count{3uz};
  constexpr auto walker_count{2uz};
  Particles particles{particle_count, walker_count};

  const std::array positions{1.0_fp, 2.0_fp, 3.0_fp};
  const std::array laplacians{4.0_fp, 5.0_fp, 6.0_fp};
  xpu::copy_n(
    particles.pos(0uz)[idx(Axis::X)],
    positions.data(),
    positions.size()
  );
  xpu::copy_n(
    particles.derivatives(0uz)[idx(Derivatives::LAP)],
    laplacians.data(),
    laplacians.size()
  );

  std::array<fp_t, particle_count> second_positions{};
  std::array<fp_t, particle_count> second_laplacians{};
  xpu::copy_n(
    second_positions.data(),
    particles.pos(1uz)[idx(Axis::X)],
    second_positions.size()
  );
  xpu::copy_n(
    second_laplacians.data(),
    particles.derivatives(1uz)[idx(Derivatives::LAP)],
    second_laplacians.size()
  );

  REQUIRE(particles.count() == particle_count);
  REQUIRE(particles.walker_count() == walker_count);
  for (auto particle{0uz}; particle < particle_count; ++particle) {
    REQUIRE(second_positions[particle] == 0.0_fp);
    REQUIRE(second_laplacians[particle] == 0.0_fp);
  }
}
