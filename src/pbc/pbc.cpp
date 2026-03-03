#include "pbc.hpp"
#include <cmath>


// store L and precompute invL
periodicBoundaryCondition::periodicBoundaryCondition(double L) noexcept
: L_{L}
, invL_{1.0 / L}
{ }

// wrap coordinates into [0, L)
double periodicBoundaryCondition::wrap(double x) const noexcept {
    double k = std::floor(x * invL_); // floor of x / L. done this way to avoid dividing by L
    x = x - k * L_;

    // fix for rare float error putting x at L
    // avoids +/- epsilon issues with floating arithmatic
    if (x >= L_) x -= L_;
    if (x < 0.0) x += L_;

    return x;
}

// wrap x, y, z indep
void periodicBoundaryCondition::wrap3(double &x, double &y, double &z) const noexcept {
    x = wrap(x);
    y = wrap(y);
    z = wrap(z);
}

// apply min image mapping to displacement component
double periodicBoundaryCondition::minImage(double dx) const noexcept {
    // dx -= L * round(dx / L)
    dx -= L_ * std::round(dx * invL_);

    // handle edge cases, ensure (-L/2, L/2]
    const double halfL = 0.5 * L_;
    if (dx > halfL) dx -= L_;
    if (dx <= -halfL) dx += L_;

    return dx;
}

// compute min image displacement vector from j to i
void periodicBoundaryCondition::displacement(
    double xi, double yi, double zi,
    double xj, double yj, double zj,
    double& dx, double& dy, double& dz) const noexcept {
    dx = xi - xj;
    dy = yi - yj;
    dz = zi - zj;

    dx = minImage(dx);
    dy = minImage(dy);
    dz = minImage(dz);
}

// compute euclidean norm of min image displacement
double periodicBoundaryCondition::distance(
    double xi, double yi, double zi,
    double xj, double yj, double zj) const noexcept {

    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;

    displacement(xi, yi, zi, xj, yj, zj, dx, dy, dz);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
