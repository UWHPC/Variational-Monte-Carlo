#pragma once

class PeriodicBoundaryCondition {
private:
    double L_{};    // Box length
    double invL_{}; // evaluates to 1.0 / L. avoids repeated division
public:
    // store L and precompute invL
    explicit PeriodicBoundaryCondition(double L) noexcept;

    // Getters for L and invL
    [[nodiscard]] double L() const noexcept { return L_; }
    [[nodiscard]] double invL() const noexcept { return invL_; }

    // wrap coordinates into [0, L)
    [[nodiscard]] double wrap(double x) const noexcept;

    // wrap x, y, z independently
    void wrap3(
        double& x,
        double& y,
        double& z
    ) const noexcept;

    // apply min image mapping to displacement component
    [[nodiscard]] double minImage(double dx) const noexcept;

    // compute min image displacement vector from j to i
    void displacement(
        double xi,
        double yi,
        double zi,
        double xj,
        double yj,
        double zj,
        double& dx,
        double& dy,
        double& dz
    ) const noexcept;

    // compute euclidean norm of min image displacement
    [[nodiscard]] double distance(
        double xi,
        double yi,
        double zi,
        double xj,
        double yj,
        double zj
    ) const noexcept;
};
