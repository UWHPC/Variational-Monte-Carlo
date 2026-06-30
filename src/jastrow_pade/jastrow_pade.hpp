#pragma once

#include "../particles/particles.hpp"

class JastrowPade {
private:
  double box_length_;
  double a_;
  double b_;

public:
  explicit JastrowPade(double box_length, double a = 0.25, double b = 1) noexcept
  : box_length_{box_length}
  , a_{a}
  , b_{b}
  { }

  [[nodiscard]] double box_length() const noexcept { return box_length_; }
  [[nodiscard]] double a() const noexcept { return a_; }
  [[nodiscard]] double b() const noexcept { return b_; }

  [[nodiscard]] double value(const Particles& particles) const noexcept;

  void add_derivatives(
    const Particles& particles,
    double* RESTRICT grad_x,
    double* RESTRICT grad_y,
    double* RESTRICT grad_z,
    double* RESTRICT laplacian
  ) const noexcept;

  [[nodiscard]] double delta_value(
    const Particles& particles,
    std::size_t moved,
    double old_x,
    double old_y,
    double old_z
  ) const noexcept;

  void update_derivatives_for_move(
    const Particles& particles,
    std::size_t moved,
    double old_x,
    double old_y,
    double old_z,
    double* RESTRICT grad_x,
    double* RESTRICT grad_y,
    double* RESTRICT grad_z,
    double* RESTRICT laplacian
  ) const noexcept;
};
