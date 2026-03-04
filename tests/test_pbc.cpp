#include <catch2/catch_test_macros.hpp>

#include "pbc/pbc.hpp"

#include <cmath>

namespace {

void requireNearPbc(double actual, double expected, double tolerance = 1e-12) {
    REQUIRE(std::abs(actual - expected) <= tolerance);
}

} // namespace

TEST_CASE("wrap maps values into [0, L)", "[pbc]") {
    const PeriodicBoundaryCondition pbc{10.0};

    requireNearPbc(pbc.wrap(0.0), 0.0);
    requireNearPbc(pbc.wrap(10.0), 0.0);
    requireNearPbc(pbc.wrap(20.0), 0.0);
    requireNearPbc(pbc.wrap(-0.25), 9.75);
    requireNearPbc(pbc.wrap(-10.0), 0.0);
}

TEST_CASE("wrap3 wraps each coordinate independently", "[pbc]") {
    const PeriodicBoundaryCondition pbc{10.0};

    double x{12.5};
    double y{-0.5};
    double z{30.0};

    pbc.wrap3(x, y, z);

    requireNearPbc(x, 2.5);
    requireNearPbc(y, 9.5);
    requireNearPbc(z, 0.0);
}

TEST_CASE("minImage enforces interval (-L/2, L/2]", "[pbc]") {
    const PeriodicBoundaryCondition pbc{10.0};

    requireNearPbc(pbc.minImage(6.0), -4.0);
    requireNearPbc(pbc.minImage(-6.0), 4.0);
    requireNearPbc(pbc.minImage(5.0), 5.0);
    requireNearPbc(pbc.minImage(-5.0), 5.0);
    requireNearPbc(pbc.minImage(15.0), 5.0);
    requireNearPbc(pbc.minImage(-15.0), 5.0);
}

TEST_CASE("displacement and distance use minimum-image coordinates", "[pbc]") {
    const PeriodicBoundaryCondition pbc{10.0};

    double dx{};
    double dy{};
    double dz{};
    pbc.displacement(9.0, 1.0, 1.0, 1.0, 1.0, 1.0, dx, dy, dz);

    requireNearPbc(dx, -2.0);
    requireNearPbc(dy, 0.0);
    requireNearPbc(dz, 0.0);

    requireNearPbc(pbc.distance(9.0, 1.0, 1.0, 1.0, 1.0, 1.0), 2.0);
    requireNearPbc(pbc.distance(1.0, 1.0, 1.0, 9.0, 1.0, 1.0), 2.0);
}
