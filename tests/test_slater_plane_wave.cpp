#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"
#include "pbc/pbc.hpp"
#include "slater_plane_wave/slater_plane_wave.hpp"

#include <cmath>
#include <numbers>
#include <vector>

namespace {

void requireNearSlater(double actual, double expected, double tolerance = 1e-10) {
    REQUIRE(std::abs(actual - expected) <= tolerance);
}

std::size_t roundUpToSimd(std::size_t n) {
    const std::size_t doublesPerAlignment{SIMD_BYTES / sizeof(double)};
    return (n + doublesPerAlignment - 1) & ~(doublesPerAlignment - 1);
}

std::size_t matrixIndex(std::size_t row, std::size_t col, std::size_t n) { return row * n + col; }

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

TEST_CASE("logAbsDet computes an inverse satisfying D*invD = I", "[slater]") {
    constexpr std::size_t N{3U};
    SlaterPlaneWave slater{N, 11.0};
    Particles particles{N};
    const PeriodicBoundaryCondition pbc{11.0};

    slater.kx()[0] = 0.0;
    slater.ky()[0] = 0.0;
    slater.kz()[0] = 0.0;

    slater.kx()[1] = 1.0;
    slater.ky()[1] = 0.5;
    slater.kz()[1] = 0.0;

    slater.kx()[2] = -0.25;
    slater.ky()[2] = 1.1;
    slater.kz()[2] = 0.7;

    particles.posX()[0] = 0.3;
    particles.posY()[0] = 0.4;
    particles.posZ()[0] = 0.5;

    particles.posX()[1] = 1.7;
    particles.posY()[1] = 0.2;
    particles.posZ()[1] = 2.1;

    particles.posX()[2] = 2.2;
    particles.posY()[2] = 1.8;
    particles.posZ()[2] = 0.9;

    const double logDet{slater.logAbsDet(particles, pbc)};
    REQUIRE(std::isfinite(logDet));

    for (std::size_t row = 0; row < N; ++row) {
        for (std::size_t col = 0; col < N; ++col) {
            double value{};
            for (std::size_t k = 0; k < N; ++k) {
                value += slater.D()[matrixIndex(row, k, N)] * slater.invD()[matrixIndex(k, col, N)];
            }
            const double expected{row == col ? 1.0 : 0.0};
            requireNearSlater(value, expected, 1e-9);
        }
    }
}

TEST_CASE("Slater derivatives match analytic N=1 formulas", "[slater]") {
    SlaterPlaneWave slater{1U, 10.0};
    Particles particles{1U};
    const PeriodicBoundaryCondition pbc{10.0};

    slater.kx()[0] = 1.2;
    slater.ky()[0] = -0.4;
    slater.kz()[0] = 0.7;

    particles.posX()[0] = 0.3;
    particles.posY()[0] = 0.5;
    particles.posZ()[0] = 0.9;

    const double logDet{slater.logAbsDet(particles, pbc)};
    REQUIRE(std::isfinite(logDet));

    const std::size_t stride{particles.paddingStride()};
    std::vector<double> gradX(stride, 2.0);
    std::vector<double> gradY(stride, 2.0);
    std::vector<double> gradZ(stride, 2.0);
    std::vector<double> lap(stride, -3.0);

    slater.addDerivatives(particles, pbc, gradX.data(), gradY.data(), gradZ.data(), lap.data());

    const double kdotr{slater.kx()[0] * particles.posX()[0] + slater.ky()[0] * particles.posY()[0] +
                       slater.kz()[0] * particles.posZ()[0]};
    const double tangent{std::tan(kdotr)};
    const double kSquared{slater.kx()[0] * slater.kx()[0] + slater.ky()[0] * slater.ky()[0] +
                          slater.kz()[0] * slater.kz()[0]};
    const double cosine{std::cos(kdotr)};

    const double expectedGradX{-tangent * slater.kx()[0]};
    const double expectedGradY{-tangent * slater.ky()[0]};
    const double expectedGradZ{-tangent * slater.kz()[0]};
    const double expectedLap{-kSquared / (cosine * cosine)};

    requireNearSlater(gradX[0], 2.0 + expectedGradX, 1e-10);
    requireNearSlater(gradY[0], 2.0 + expectedGradY, 1e-10);
    requireNearSlater(gradZ[0], 2.0 + expectedGradZ, 1e-10);
    requireNearSlater(lap[0], -3.0 + expectedLap, 1e-10);

    for (std::size_t i = 1; i < stride; ++i) {
        requireNearSlater(gradX[i], 2.0);
        requireNearSlater(gradY[i], 2.0);
        requireNearSlater(gradZ[i], 2.0);
        requireNearSlater(lap[i], -3.0);
    }
}
