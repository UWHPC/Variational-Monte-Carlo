#include <catch2/catch_test_macros.hpp>

#include "jastrow_pade/jastrow_pade.hpp"

#include <cmath>
#include <vector>

namespace {

void requireNearJastrow(double actual, double expected, double tolerance = 1e-12) {
    REQUIRE(std::abs(actual - expected) <= tolerance);
}

} // namespace

TEST_CASE("Jastrow value uses minimum-image pair distances", "[jastrow]") {
    const JastrowPade jastrow{0.5, 1.0};
    const PeriodicBoundaryCondition pbc{10.0};
    Particles particles{2U};

    particles.posX()[0] = 0.1;
    particles.posY()[0] = 0.0;
    particles.posZ()[0] = 0.0;

    particles.posX()[1] = 9.9;
    particles.posY()[1] = 0.0;
    particles.posZ()[1] = 0.0;

    const double r{0.2};
    const double expected{(0.5 * r) / (1.0 + r)};
    requireNearJastrow(jastrow.value(particles, pbc), expected);
}

TEST_CASE("Jastrow value skips degenerate pairs", "[jastrow]") {
    const JastrowPade jastrow{0.5, 1.0};
    const PeriodicBoundaryCondition pbc{10.0};
    Particles particles{2U};

    particles.posX()[0] = 1.0;
    particles.posY()[0] = 2.0;
    particles.posZ()[0] = 3.0;

    particles.posX()[1] = 1.0;
    particles.posY()[1] = 2.0;
    particles.posZ()[1] = 3.0;

    requireNearJastrow(jastrow.value(particles, pbc), 0.0);
}

TEST_CASE("Jastrow derivatives match the analytic two-particle result", "[jastrow]") {
    const JastrowPade jastrow{0.5, 1.0};
    const PeriodicBoundaryCondition pbc{100.0};
    Particles particles{2U};

    particles.posX()[0] = 0.0;
    particles.posY()[0] = 0.0;
    particles.posZ()[0] = 0.0;

    particles.posX()[1] = 1.0;
    particles.posY()[1] = 0.0;
    particles.posZ()[1] = 0.0;

    const std::size_t stride{particles.paddingStride()};
    std::vector<double> gradX(stride, 0.0);
    std::vector<double> gradY(stride, 0.0);
    std::vector<double> gradZ(stride, 0.0);
    std::vector<double> lap(stride, 0.0);

    jastrow.addDerivatives(particles, pbc, gradX.data(), gradY.data(), gradZ.data(), lap.data());

    requireNearJastrow(gradX[0], -0.125);
    requireNearJastrow(gradX[1], 0.125);
    requireNearJastrow(gradX[0] + gradX[1], 0.0);

    requireNearJastrow(gradY[0], 0.0);
    requireNearJastrow(gradY[1], 0.0);
    requireNearJastrow(gradZ[0], 0.0);
    requireNearJastrow(gradZ[1], 0.0);

    requireNearJastrow(lap[0], 0.125);
    requireNearJastrow(lap[1], 0.125);

    for (std::size_t i = 2; i < stride; ++i) {
        requireNearJastrow(gradX[i], 0.0);
        requireNearJastrow(gradY[i], 0.0);
        requireNearJastrow(gradZ[i], 0.0);
        requireNearJastrow(lap[i], 0.0);
    }
}

TEST_CASE("Jastrow derivatives are unchanged for degenerate pairs", "[jastrow]") {
    const JastrowPade jastrow{0.5, 1.0};
    const PeriodicBoundaryCondition pbc{10.0};
    Particles particles{2U};

    particles.posX()[0] = 4.0;
    particles.posY()[0] = 5.0;
    particles.posZ()[0] = 6.0;

    particles.posX()[1] = 4.0;
    particles.posY()[1] = 5.0;
    particles.posZ()[1] = 6.0;

    const std::size_t stride{particles.paddingStride()};
    std::vector<double> gradX(stride, 3.0);
    std::vector<double> gradY(stride, -2.0);
    std::vector<double> gradZ(stride, 1.5);
    std::vector<double> lap(stride, 7.0);

    jastrow.addDerivatives(particles, pbc, gradX.data(), gradY.data(), gradZ.data(), lap.data());

    for (std::size_t i = 0; i < stride; ++i) {
        requireNearJastrow(gradX[i], 3.0);
        requireNearJastrow(gradY[i], -2.0);
        requireNearJastrow(gradZ[i], 1.5);
        requireNearJastrow(lap[i], 7.0);
    }
}
