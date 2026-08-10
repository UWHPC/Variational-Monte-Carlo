#pragma once

#include <xpu/math.hpp>
#include "../utilities/macros.hpp"
#include "../utilities/components.hpp"

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

  [[nodiscard]] real_t value(xpu::soa_view<real_t, idx(Axis::NUM)> pos) const noexcept;

  void add_derivatives(
    xpu::soa_view<real_t, idx(Axis::NUM)> pos,
    xpu::soa_view<real_t, idx(Derivatives::NUM)> derivatives
  ) const noexcept;

  [[nodiscard]] real_t delta_value(
    std::size_t moved,
    xpu::array<real_t, idx(Axis::NUM)> old_pos,
    xpu::soa_view<real_t, idx(Axis::NUM)> pos
  ) const noexcept;

  void update_derivatives_for_move(
    std::size_t moved,
    xpu::array<real_t, idx(Axis::NUM)> old_pos,
    xpu::soa_view<real_t, idx(Axis::NUM)> pos,
    xpu::soa_view<real_t, idx(Derivatives::NUM)> derivatives
  ) const noexcept;
};
