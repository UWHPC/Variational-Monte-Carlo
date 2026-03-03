#include "pbc.hpp"

#include <cmath>

// store L and precompute invL
PeriodicBoundaryCondition::PeriodicBoundaryCondition(double L) noexcept
: L_{L}
, invL_{1.0 / L}
{ }

// wrap coordinates into [0, L)
double PeriodicBoundaryCondition::wrap(double x) const noexcept {
    double k{std::floor(x * invL())}; // floor of x / L. done this way to avoid dividing by L
    x = x - k*L();

    // fix for rare float error putting x at L
    // avoids +/- epsilon issues with floating arithmatic
    if (x >= L()) x -= L();
    if (x < 0.0) x += L();

    return x;
}

// wrap x, y, z indep
void PeriodicBoundaryCondition::wrap3(double &x, double &y, double &z) const noexcept {
    x = wrap(x);
    y = wrap(y);
    z = wrap(z);
}

// apply min image mapping to displacement component
double PeriodicBoundaryCondition::minImage(double dx) const noexcept {
    // dx -= L * round(dx / L)
    dx -= L() * std::round(dx * invL());

    // handle edge cases, ensure (-L/2, L/2]
    const double halfL{0.5 * L()};
    if (dx > halfL) dx -= L();
    if (dx <= -halfL) dx += L();

    return dx;
}

// compute min image displacement vector from j to i
void PeriodicBoundaryCondition::displacement(
    double xi, double yi, double zi,
    double xj, double yj, double zj,
    double& dx, double& dy, double& dz
    ) const noexcept {
    dx = xi - xj;
    dy = yi - yj;
    dz = zi - zj;

    dx = minImage(dx);
    dy = minImage(dy);
    dz = minImage(dz);
}

// compute euclidean norm of min image displacement
double PeriodicBoundaryCondition::distance(
    double xi, double yi, double zi,
    double xj, double yj, double zj
    ) const noexcept {

    double dx{};
    double dy{};
    double dz{};

    displacement(
        xi, yi, zi, 
        xj, yj, zj, 
        dx, dy, dz
    );

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}
