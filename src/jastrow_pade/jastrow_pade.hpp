#pragma once

#include <xpu/math.hpp>
#include "../utilities/macros.hpp"
#include "../utilities/components.hpp"

class JastrowPade {
private:
  fp_t box_length_;
  fp_t a_;
  fp_t b_;

public:
  explicit JastrowPade(fp_t box_length, fp_t a = 0.25_fp, fp_t b = 1.0_fp) noexcept
    : box_length_{box_length}
    , a_{a}
    , b_{b}
  { }

  [[nodiscard]] fp_t box_length() const noexcept { return box_length_; }
  [[nodiscard]] fp_t a() const noexcept { return a_; }
  [[nodiscard]] fp_t b() const noexcept { return b_; }

  [[nodiscard]] fp_t value(xpu::soa_view<fp_t, idx(Axis::NUM)> pos) const noexcept;

  void add_derivatives(
    xpu::soa_view<fp_t, idx(Axis::NUM)> pos,
    xpu::soa_view<fp_t, idx(Derivatives::NUM)> derivatives
  ) const noexcept;

  [[nodiscard]] fp_t delta_value(
    std::size_t moved,
    xpu::array<fp_t, idx(Axis::NUM)> old_pos,
    xpu::soa_view<fp_t, idx(Axis::NUM)> pos
  ) const noexcept;

  void update_derivatives_for_move(
    std::size_t moved,
    xpu::array<fp_t, idx(Axis::NUM)> old_pos,
    xpu::soa_view<fp_t, idx(Axis::NUM)> pos,
    xpu::soa_view<fp_t, idx(Derivatives::NUM)> derivatives
  ) const noexcept;
};
