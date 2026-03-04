#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"
#include "pbc/pbc.hpp"
#include "slater_plane_wave/slater_plane_wave.hpp"

#include <cmath>
#include <numbers>

namespace {

void requireNearSlater(double actual, double expected, double tolerance = 1e-10) {
    REQUIRE(std::abs(actual - expected) <= tolerance);
}

std::size_t roundUpToSimd(std::size_t n) {
    const std::size_t doublesPerAlignment{SIMD_BYTES / sizeof(double)};
    return (n + doublesPerAlignment - 1) & ~(doublesPerAlignment - 1);
}

} // namespace

TEST_CASE("SlaterPlaneWave constructor rounds strides to SIMD boundaries", "[slater]") {
    constexpr std::size_t N{3U};
    SlaterPlaneWave slater{N, 5.0};

    REQUIRE(slater.N() == N);
    REQUIRE(slater.L() == 5.0);
    REQUIRE(slater.vecStride() == roundUpToSimd(N));
    REQUIRE(slater.matStride() == roundUpToSimd(N * N));
}

TEST_CASE("logAbsDet handles the N=1 constant orbital case", "[slater]") {
    SlaterPlaneWave slater{1U, 10.0};
    Particles particles{1U};
    const PeriodicBoundaryCondition pbc{10.0};

    slater.kx()[0] = 0.0;
    slater.ky()[0] = 0.0;
    slater.kz()[0] = 0.0;

    particles.posX()[0] = 3.25;
    particles.posY()[0] = 1.50;
    particles.posZ()[0] = 7.75;

    const double logDet{slater.logAbsDet(particles, pbc)};

    requireNearSlater(logDet, 0.0);
    requireNearSlater(slater.D()[0], 1.0);
    requireNearSlater(slater.invD()[0], 1.0);
}

TEST_CASE("logAbsDet returns negative infinity for singular matrices", "[slater]") {
    SlaterPlaneWave slater{2U, 10.0};
    Particles particles{2U};
    const PeriodicBoundaryCondition pbc{10.0};

    slater.kx()[0] = 0.0;
    slater.ky()[0] = 0.0;
    slater.kz()[0] = 0.0;
    slater.kx()[1] = 0.0;
    slater.ky()[1] = 0.0;
    slater.kz()[1] = 0.0;

    particles.posX()[0] = 0.0;
    particles.posY()[0] = 0.0;
    particles.posZ()[0] = 0.0;
    particles.posX()[1] = 1.0;
    particles.posY()[1] = 0.0;
    particles.posZ()[1] = 0.0;

    const double logDet{slater.logAbsDet(particles, pbc)};

    REQUIRE(std::isinf(logDet));
    REQUIRE(logDet < 0.0);
}

TEST_CASE("logAbsDet computes determinant and inverse for a known 2x2 system", "[slater]") {
    SlaterPlaneWave slater{2U, 10.0};
    Particles particles{2U};
    const PeriodicBoundaryCondition pbc{10.0};

    slater.kx()[0] = 0.0;
    slater.ky()[0] = 0.0;
    slater.kz()[0] = 0.0;

    slater.kx()[1] = 1.0;
    slater.ky()[1] = 0.0;
    slater.kz()[1] = 0.0;

    particles.posX()[0] = 0.0;
    particles.posY()[0] = 0.0;
    particles.posZ()[0] = 0.0;

    particles.posX()[1] = std::numbers::pi_v<double> * 0.5;
    particles.posY()[1] = 0.0;
    particles.posZ()[1] = 0.0;

    const double logDet{slater.logAbsDet(particles, pbc)};

    requireNearSlater(logDet, 0.0, 1e-9);

    requireNearSlater(slater.D()[0], 1.0);
    requireNearSlater(slater.D()[1], 1.0);
    requireNearSlater(slater.D()[2], 1.0);
    requireNearSlater(slater.D()[3], 0.0, 1e-12);

    requireNearSlater(slater.invD()[0], 0.0, 1e-9);
    requireNearSlater(slater.invD()[1], 1.0, 1e-9);
    requireNearSlater(slater.invD()[2], 1.0, 1e-9);
    requireNearSlater(slater.invD()[3], -1.0, 1e-9);
}
