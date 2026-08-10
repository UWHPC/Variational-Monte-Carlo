#include <xpu/xpu.hpp>
#include "jastrow_pade.hpp"
#include "jastrow_pade_kernels.hpp"

#include <cstddef>

void JastrowPade::update_derivatives_for_move(
  std::size_t moved,
  xpu::array<real_t, idx(Axis::NUM)> old_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  xpu::soa_view<real_t, idx(Derivatives::NUM)> derivatives
) const noexcept {
  const auto L{box_length_};
  const auto half_L{0.5_r * L};
  const auto neg_two_ab{-2.0_r * this->a() * this->b()};

  kernel::jpade::compute_derivatives(
    moved, L, half_L,
    this->a(), this->b(), neg_two_ab,
    old_pos, pos, derivatives
  );
}

void JastrowPade::add_derivatives(
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  xpu::soa_view<real_t, idx(Derivatives::NUM)> derivatives
) const noexcept {
  const auto L{box_length_};
  const auto half_L{0.5_r * L};
  const auto neg_two_ab{-2.0_r * this->a() * this->b()};

  kernel::jpade::add_derivatives(
    L, half_L,
    this->a(), this->b(), neg_two_ab,
    pos, derivatives
  );
}

real_t JastrowPade::value(
  xpu::soa_view<real_t, idx(Axis::NUM)> pos
) const noexcept {
  const auto L{box_length_};
  const auto half_L{0.5_r * L};

  return kernel::jpade::value(
    L, half_L,
    this->a(), this->b(),
    pos
  );
}

real_t JastrowPade::delta_value(
  std::size_t moved,
  xpu::array<real_t, idx(Axis::NUM)> old_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos
) const noexcept {
  const auto L{box_length_};
  const auto half_L{0.5_r * L};

  return kernel::jpade::delta_value(
    moved,
    L, half_L,
    this->a(), this->b(),
    old_pos, pos
  );
}
