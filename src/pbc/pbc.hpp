#pragma once

class periodicBoundaryCondition {
private:
  double L_{};
  double invL_{}; // avoids repeated division for min image mapping
public:
    explicit periodicBoundaryCondition(double L) noexcept;
    [[nodiscard]] double boxLength() const noexcept { return L_; }
    [[nodiscard]] double wrap(double x) const noexcept;
    void wrap3(double &x, double &y, double &z) const noexcept;

    [[nodiscard]] double minImage(double dx) const noexcept;
    void displacement(double xi, double yi, double zi, double xj, double yj, double zj, double &dx, double &dy, double &dz) const noexcept;
    [[nodiscard]] double distance(double xi, double yi, double zi, double xj, double yj, double zj) const noexcept;
};
