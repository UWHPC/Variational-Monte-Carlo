#pragma once

#include "../particles/particles.hpp"
#include <xpu/math.hpp>

class JastrowPade {
private:
  real_t box_length_;
  real_t a_;
  real_t b_;

public:
  explicit JastrowPade(real_t box_length, real_t a = 0.25_r, real_t b = 1.0_r) noexcept
  : box_length_{box_length}
  , a_{a}
  , b_{b}
  { }

  [[nodiscard]] real_t box_length() const noexcept { return box_length_; }
  [[nodiscard]] real_t a() const noexcept { return a_; }
  [[nodiscard]] real_t b() const noexcept { return b_; }

  [[nodiscard]] real_t value(const Particles& particles) const noexcept;

  void add_derivatives(
    const Particles& particles,
    real_t* RESTRICT grad_x,
    real_t* RESTRICT grad_y,
    real_t* RESTRICT grad_z,
    real_t* RESTRICT laplacian
  ) const noexcept;

  [[nodiscard]] real_t delta_value(
    const Particles& particles,
    std::size_t moved,
    real_t old_x,
    real_t old_y,
    real_t old_z
  ) const noexcept;

  void update_derivatives_for_move(
    const Particles& particles,
    std::size_t moved,
    real_t old_x,
    real_t old_y,
    real_t old_z,
    real_t* RESTRICT grad_x,
    real_t* RESTRICT grad_y,
    real_t* RESTRICT grad_z,
    real_t* RESTRICT laplacian
  ) const noexcept;
};
