#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"

#include <cstddef>
#include <cstdint>

TEST_CASE("Particles allocates aligned padded blocks and zero-initializes them", "[particles]") {
    constexpr std::size_t numParticles{3U};
    Particles particles{numParticles};

    const std::size_t stride{particles.padding_stride_get()};
    const std::size_t doublesPerAlignment{SIMD_BYTES / sizeof(double)};

    REQUIRE(particles.num_particles_get() == numParticles);
    REQUIRE(stride >= numParticles);
    REQUIRE(stride % doublesPerAlignment == 0U);

    const auto baseAddress{reinterpret_cast<std::uintptr_t>(particles.pos_x_get())};
    REQUIRE(baseAddress % SIMD_BYTES == 0U);

    for (std::size_t i = 0; i < stride; ++i) {
        REQUIRE(particles.pos_x_get()[i] == 0.0);
        REQUIRE(particles.pos_y_get()[i] == 0.0);
        REQUIRE(particles.pos_z_get()[i] == 0.0);
        REQUIRE(particles.grad_log_psi_x_get()[i] == 0.0);
        REQUIRE(particles.grad_log_psi_y_get()[i] == 0.0);
        REQUIRE(particles.grad_log_psi_z_get()[i] == 0.0);
        REQUIRE(particles.lap_log_psi_get()[i] == 0.0);
    }
}

TEST_CASE("Particles exposes non-overlapping slices for each component", "[particles]") {
    Particles particles{2U};

    const std::ptrdiff_t stride{static_cast<std::ptrdiff_t>(particles.padding_stride_get())};
    REQUIRE(particles.pos_y_get() - particles.pos_x_get() == stride);
    REQUIRE(particles.pos_z_get() - particles.pos_y_get() == stride);
    REQUIRE(particles.grad_log_psi_x_get() - particles.pos_z_get() == stride);
    REQUIRE(particles.grad_log_psi_y_get() - particles.grad_log_psi_x_get() == stride);
    REQUIRE(particles.grad_log_psi_z_get() - particles.grad_log_psi_y_get() == stride);
    REQUIRE(particles.lap_log_psi_get() - particles.grad_log_psi_z_get() == stride);

    particles.pos_x_get()[0] = 1.0;
    particles.pos_y_get()[0] = 2.0;
    particles.pos_z_get()[0] = 3.0;

    REQUIRE(particles.pos_x_get()[0] == 1.0);
    REQUIRE(particles.pos_y_get()[0] == 2.0);
    REQUIRE(particles.pos_z_get()[0] == 3.0);
    REQUIRE(particles.grad_log_psi_x_get()[0] == 0.0);
}
