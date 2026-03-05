#include <catch2/catch_test_macros.hpp>

#include "pbc/pbc.hpp"

#include <cmath>

namespace {

void require_near_pbc(double actual, double expected, double tolerance = 1e-12) {
    REQUIRE(std::abs(actual - expected) <= tolerance);
}

} // namespace

TEST_CASE("wrap maps values into [0, L)", "[pbc]") {
    const PeriodicBoundaryCondition pbc{10.0};

    require_near_pbc(pbc.wrap(0.0), 0.0);
    require_near_pbc(pbc.wrap(10.0), 0.0);
    require_near_pbc(pbc.wrap(20.0), 0.0);
    require_near_pbc(pbc.wrap(-0.25), 9.75);
    require_near_pbc(pbc.wrap(-10.0), 0.0);
}

TEST_CASE("wrap3 wraps each coordinate independently", "[pbc]") {
    const PeriodicBoundaryCondition pbc{10.0};

    double x{12.5};
    double y{-0.5};
    double z{30.0};

    pbc.wrap3(x, y, z);

    require_near_pbc(x, 2.5);
    require_near_pbc(y, 9.5);
    require_near_pbc(z, 0.0);
}

TEST_CASE("min_image enforces interval (-L/2, L/2]", "[pbc]") {
    const PeriodicBoundaryCondition pbc{10.0};

    require_near_pbc(pbc.min_image(6.0), -4.0);
    require_near_pbc(pbc.min_image(-6.0), 4.0);
    require_near_pbc(pbc.min_image(5.0), 5.0);
    require_near_pbc(pbc.min_image(-5.0), 5.0);
    require_near_pbc(pbc.min_image(15.0), 5.0);
    require_near_pbc(pbc.min_image(-15.0), 5.0);
}

TEST_CASE("displacement and distance use minimum-image coordinates", "[pbc]") {
    const PeriodicBoundaryCondition pbc{10.0};

    double dx{};
    double dy{};
    double dz{};
    pbc.displacement(9.0, 1.0, 1.0, 1.0, 1.0, 1.0, dx, dy, dz);

    require_near_pbc(dx, -2.0);
    require_near_pbc(dy, 0.0);
    require_near_pbc(dz, 0.0);

    require_near_pbc(pbc.distance(9.0, 1.0, 1.0, 1.0, 1.0, 1.0), 2.0);
    require_near_pbc(pbc.distance(1.0, 1.0, 1.0, 9.0, 1.0, 1.0), 2.0);
}

TEST_CASE("wrap and min_image are invariant under full-box translations", "[pbc]") {
    const PeriodicBoundaryCondition pbc{10.0};
    const double L{pbc.L_ptr()};

    const double samples[]{-23.75, -10.1, -0.01, 0.0, 0.01, 4.9, 9.99, 10.01, 31.5};
    for (const double x : samples) {
        const double wrapped{pbc.wrap(x)};
        REQUIRE(wrapped >= 0.0);
        REQUIRE(wrapped < L);
        require_near_pbc(pbc.wrap(x + 2.0 * L), wrapped);
        require_near_pbc(pbc.wrap(x - 3.0 * L), wrapped);

        const double minImg{pbc.min_image(x)};
        require_near_pbc(pbc.min_image(x + 4.0 * L), minImg);
        require_near_pbc(pbc.min_image(x - 5.0 * L), minImg);
    }
}

TEST_CASE("minimum-image displacement is antisymmetric between particle orderings", "[pbc]") {
    const PeriodicBoundaryCondition pbc{10.0};

    double dxAB{};
    double dyAB{};
    double dzAB{};
    pbc.displacement(1.2, 2.7, 4.0, 8.9, 1.3, 6.4, dxAB, dyAB, dzAB);

    double dxBA{};
    double dyBA{};
    double dzBA{};
    pbc.displacement(8.9, 1.3, 6.4, 1.2, 2.7, 4.0, dxBA, dyBA, dzBA);

    require_near_pbc(dxAB, -dxBA);
    require_near_pbc(dyAB, -dyBA);
    require_near_pbc(dzAB, -dzBA);
}
