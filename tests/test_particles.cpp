#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"

#include <cstddef>
#include <cstdint>

TEST_CASE("Particles allocates aligned padded blocks and zero-initializes them", "[particles]") {
    constexpr std::size_t numParticles{3U};
    Particles particles{numParticles};

    const std::size_t stride{particles.paddingStride()};
    const std::size_t doublesPerAlignment{SIMD_BYTES / sizeof(double)};

    REQUIRE(particles.numParticles() == numParticles);
    REQUIRE(stride >= numParticles);
    REQUIRE(stride % doublesPerAlignment == 0U);

    const auto baseAddress{reinterpret_cast<std::uintptr_t>(particles.posX())};
    REQUIRE(baseAddress % SIMD_BYTES == 0U);

    for (std::size_t i = 0; i < stride; ++i) {
        REQUIRE(particles.posX()[i] == 0.0);
        REQUIRE(particles.posY()[i] == 0.0);
        REQUIRE(particles.posZ()[i] == 0.0);
        REQUIRE(particles.gradLogPsiX()[i] == 0.0);
        REQUIRE(particles.gradLogPsiY()[i] == 0.0);
        REQUIRE(particles.gradLogPsiZ()[i] == 0.0);
        REQUIRE(particles.logPsi()[i] == 0.0);
        REQUIRE(particles.laplogPsi()[i] == 0.0);
    }
}

TEST_CASE("Particles exposes non-overlapping slices for each component", "[particles]") {
    Particles particles{2U};

    const std::ptrdiff_t stride{static_cast<std::ptrdiff_t>(particles.paddingStride())};
    REQUIRE(particles.posY() - particles.posX() == stride);
    REQUIRE(particles.posZ() - particles.posY() == stride);
    REQUIRE(particles.gradLogPsiX() - particles.posZ() == stride);
    REQUIRE(particles.gradLogPsiY() - particles.gradLogPsiX() == stride);
    REQUIRE(particles.gradLogPsiZ() - particles.gradLogPsiY() == stride);
    REQUIRE(particles.logPsi() - particles.gradLogPsiZ() == stride);
    REQUIRE(particles.laplogPsi() - particles.logPsi() == stride);

    particles.posX()[0] = 1.0;
    particles.posY()[0] = 2.0;
    particles.posZ()[0] = 3.0;

    REQUIRE(particles.posX()[0] == 1.0);
    REQUIRE(particles.posY()[0] == 2.0);
    REQUIRE(particles.posZ()[0] == 3.0);
    REQUIRE(particles.gradLogPsiX()[0] == 0.0);
}
