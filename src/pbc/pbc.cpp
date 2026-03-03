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
}

// apply min image mapping to displacement component
double periodicBoundaryCondition::minImage(double dx) const noexcept {
}

// compute min image displacement vector from j to i
void periodicBoundaryCondition::displacement(
    double xi, double yi, double zi,
    double xj, double yj, double zj,
    double& dx, double& dy, double& dz) const noexcept {
}

// compute euclidean norm of min image displacement
double periodicBoundaryCondition::distance(
    double xi, double yi, double zi,
    double xj, double yj, double zj) const noexcept {
}
