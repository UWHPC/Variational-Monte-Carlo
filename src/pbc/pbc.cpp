#include "pbc.hpp"

#include <cmath>

PeriodicBoundaryCondition::PeriodicBoundaryCondition(double L) noexcept : L_{L}, inv_L_{1.0 / L} {}

double PeriodicBoundaryCondition::wrap(double x) const noexcept {
    const double L{L_get()};
    const double inv_L{inv_L_get()};

    // floor of x / L. done this way to avoid dividing by L
    const double K{std::floor(x * inv_L)};

    x -= K * L;

    // fix for rare float error putting x at L
    // avoids +/- epsilon issues with floating arithmatic
    if (x >= L) {
        x -= L;
    }
    if (x < 0.0) {
        x += L;
    }

    return x;
}

void PeriodicBoundaryCondition::wrap3(double& x, double& y, double& z) const noexcept {
    x = wrap(x);
    y = wrap(y);
    z = wrap(z);
}