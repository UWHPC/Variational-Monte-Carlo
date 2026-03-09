#pragma once

class PeriodicBoundaryCondition {
private:
    double L_{};     // Box length
    double inv_L_{}; // evaluates to 1.0 / L. avoids repeated division
public:
    // store L and precompute invL
    explicit PeriodicBoundaryCondition(double L) noexcept;

    // Getters for L and invL
    [[nodiscard]] double L_get() const noexcept { return L_; }
    [[nodiscard]] double inv_L_get() const noexcept { return inv_L_; }

    // wrap coordinates into [0, L)
    [[nodiscard]] double wrap(double x) const noexcept;

    // wrap x, y, z independently
    void wrap3(double& x, double& y, double& z) const noexcept;
};
