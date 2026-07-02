#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"

#include <cstddef>
#include <cstdint>

TEST_CASE("Particles allocates aligned padded blocks and zero-initializes them", "[particles]") {
  constexpr std::size_t numParticles{3U};
  Particles particles{numParticles};

  const std::size_t stride{particles.p_stride()};
  const std::size_t doublesPerAlignment{SIMD_BYTES / sizeof(real_t)};

  REQUIRE(particles.size() == numParticles);
  REQUIRE(stride >= numParticles);
  REQUIRE(stride % doublesPerAlignment == 0U);

  const auto baseAddress{reinterpret_cast<std::uintptr_t>(particles.pos().x_)};
  REQUIRE(baseAddress % SIMD_BYTES == 0U);

  for (std::size_t i = 0; i < stride; ++i) {
    REQUIRE(particles.pos().x_[i] == 0.0_r);
    REQUIRE(particles.pos().y_[i] == 0.0_r);
    REQUIRE(particles.pos().z_[i] == 0.0_r);
    REQUIRE(particles.grad_log_psi().x_[i] == 0.0_r);
    REQUIRE(particles.grad_log_psi().y_[i] == 0.0_r);
    REQUIRE(particles.grad_log_psi().z_[i] == 0.0_r);
    REQUIRE(particles.lap_log_psi()[i] == 0.0_r);
  }
}

TEST_CASE("Particles exposes non-overlapping slices for each component", "[particles]") {
  Particles particles{2U};

  const std::ptrdiff_t stride{static_cast<std::ptrdiff_t>(particles.p_stride())};
  REQUIRE(particles.pos().y_ - particles.pos().x_ == stride);
  REQUIRE(particles.pos().z_ - particles.pos().y_ == stride);
  REQUIRE(particles.grad_log_psi().x_ - particles.pos().z_ == stride);
  REQUIRE(particles.grad_log_psi().y_ - particles.grad_log_psi().x_ == stride);
  REQUIRE(particles.grad_log_psi().z_ - particles.grad_log_psi().y_ == stride);
  REQUIRE(particles.lap_log_psi() - particles.grad_log_psi().z_ == stride);

  particles.pos().x_[0] = 1.0_r;
  particles.pos().y_[0] = 2.0_r;
  particles.pos().z_[0] = 3.0_r;

  REQUIRE(particles.pos().x_[0] == 1.0_r);
  REQUIRE(particles.pos().y_[0] == 2.0_r);
  REQUIRE(particles.pos().z_[0] == 3.0_r);
  REQUIRE(particles.grad_log_psi().x_[0] == 0.0_r);
}
