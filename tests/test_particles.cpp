#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>

TEST_CASE("Particles allocates aligned padded arrays", "[particles]") {
  constexpr auto num_particles{3uz};
  Particles particles{num_particles};

  const auto stride{particles.stride()};
  const auto elements_per_alignment{xpu::simd_bytes / sizeof(fp_t)};
  const auto pos{particles.pos()};

  REQUIRE(particles.count() == num_particles);
  REQUIRE(particles.walker_count() == 1uz);
  REQUIRE(stride >= num_particles);
  REQUIRE(stride % elements_per_alignment == 0uz);

  const auto base_address{reinterpret_cast<std::uintptr_t>(pos[idx(Axis::X)])};
  REQUIRE(base_address % xpu::simd_bytes == 0uz);
}

TEST_CASE("Particles exposes non-overlapping component arrays", "[particles]") {
  Particles particles{2uz};

  const auto stride{static_cast<std::ptrdiff_t>(particles.stride())};
  auto pos{particles.pos()};
  auto derivatives{particles.derivatives()};

  REQUIRE(pos[idx(Axis::Y)] - pos[idx(Axis::X)] == stride);
  REQUIRE(pos[idx(Axis::Z)] - pos[idx(Axis::Y)] == stride);
  REQUIRE(
    derivatives[idx(Derivatives::GRAD_Y)]
    - derivatives[idx(Derivatives::GRAD_X)] == stride
  );
  REQUIRE(
    derivatives[idx(Derivatives::GRAD_Z)]
    - derivatives[idx(Derivatives::GRAD_Y)] == stride
  );
  REQUIRE(
    derivatives[idx(Derivatives::LAP)]
    - derivatives[idx(Derivatives::GRAD_Z)] == stride
  );
}

TEST_CASE("Particles keeps walker storage isolated", "[particles]") {
  constexpr auto num_particles{3uz};
  constexpr auto num_walkers{2uz};
  Particles particles{num_particles, num_walkers};

  auto first_pos{particles.pos(0uz)};
  auto second_pos{particles.pos(1uz)};
  auto first_derivatives{particles.derivatives(0uz)};
  auto second_derivatives{particles.derivatives(1uz)};

  first_pos[idx(Axis::X)][0uz] = 1.0_fp;
  first_derivatives[idx(Derivatives::LAP)][1uz] = 2.0_fp;

  REQUIRE(particles.count() == num_particles);
  REQUIRE(particles.walker_count() == num_walkers);
  REQUIRE(second_pos[idx(Axis::X)][0uz] == 0.0_fp);
  REQUIRE(second_derivatives[idx(Derivatives::LAP)][1uz] == 0.0_fp);

  const auto& const_particles{std::as_const(particles)};
  REQUIRE(const_particles.pos(0uz)[idx(Axis::X)][0uz] == 1.0_fp);
  REQUIRE(
    const_particles.derivatives(0uz)[idx(Derivatives::LAP)][1uz] == 2.0_fp
  );
}
