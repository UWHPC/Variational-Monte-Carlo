#pragma once

#include <xpu/math.hpp>
#include "../particles/particles.hpp"
#include "../utilities/macros.hpp"
#include "../utilities/components.hpp"

class JastrowPade {
private:
  fp_t box_length_;
  fp_t a_;
  fp_t b_;

public:
  struct View {
    fp_t box_length;
    fp_t a;
    fp_t b;
  };

  explicit JastrowPade(fp_t box_length, fp_t a = 0.25_fp, fp_t b = 1.0_fp) noexcept
    : box_length_{box_length}
    , a_{a}
    , b_{b}
  { }

  [[nodiscard]] CUDA_CALLABLE fp_t box_length() const noexcept { return box_length_; }
  [[nodiscard]] CUDA_CALLABLE fp_t a() const noexcept { return a_; }
  [[nodiscard]] CUDA_CALLABLE fp_t b() const noexcept { return b_; }

  [[nodiscard]] CUDA_CALLABLE
  View view() const noexcept {
    return {
      this->box_length(),
      this->a(),
      this->b()
    };
  }

  [[nodiscard]] fp_t value(Particles::View particles) const noexcept;

  void add_derivatives(
    Particles::View particles,
    xpu::soa_view<fp_t, idx(Derivatives::NUM)> derivatives
  ) const noexcept;

  [[nodiscard]] fp_t delta_value(
    std::size_t moved,
    xpu::array<fp_t, idx(Axis::NUM)> old_pos,
    Particles::View particles
  ) const noexcept;

  void update_derivatives_for_move(
    std::size_t moved,
    xpu::array<fp_t, idx(Axis::NUM)> old_pos,
    Particles::View particles,
    xpu::soa_view<fp_t, idx(Derivatives::NUM)> derivatives
  ) const noexcept;
};
