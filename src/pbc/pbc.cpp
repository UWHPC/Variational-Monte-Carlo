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

double PeriodicBoundaryCondition::min_image(double dx) const noexcept {
    const double L{L_get()};
    const double inv_L{inv_L_get()};

    // dx -= L * round(dx / L)
    dx -= L * std::round(dx * inv_L);

    // handle edge cases, ensure (-L/2, L/2]
    const double HALF_L{0.5 * L};

    if (dx > HALF_L) {
        dx -= L;
    }
    if (dx <= -HALF_L) {
        dx += L;
    }

    return dx;
}

void PeriodicBoundaryCondition::displacement(double xi, double yi, double zi, double xj, double yj, double zj,
                                             double& dx, double& dy, double& dz) const noexcept {
    dx = xi - xj;
    dy = yi - yj;
    dz = zi - zj;

    dx = min_image(dx);
    dy = min_image(dy);
    dz = min_image(dz);
}

double PeriodicBoundaryCondition::distance(double xi, double yi, double zi, double xj, double yj,
                                           double zj) const noexcept {
    double dx{}, dy{}, dz{};

    displacement(xi, yi, zi, xj, yj, zj, dx, dy, dz);

    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
