#pragma once

#include "../jastrow_pade/jastrow_pade.hpp"
#include "../particles/particles.hpp"
#include "../slater_plane_wave/slater_plane_wave.hpp"
#include <xpu/soa.hpp>

#include <cstddef>

class WaveFunction {
private:
  JastrowPade jastrow_pade_;
  SlaterPlaneWave slater_plane_wave_;

  bool jastrow_cache_valid_;
  std::size_t steps_since_refresh_;

  enum class ArrayIndex : std::size_t {
    DERIVATIVES,
    NUM_ARRAYS = enum_index(ArrayIndex::DERIVATIVES, Derivatives::NUM)
  };
  xpu::soa<fp_t, idx(ArrayIndex::NUM_ARRAYS)> deriv_;

public:
  explicit WaveFunction(
    const Particles& particles,
    fp_t box_length,
    fp_t a = 0.25_fp,
    fp_t b = 1.0_fp
  )
    : jastrow_pade_{box_length, a, b}
    , slater_plane_wave_{particles, box_length}
    , jastrow_cache_valid_{}
    , steps_since_refresh_{}
    , deriv_{particles.count()}
  { }

  [[nodiscard]]       JastrowPade& jastrow_pade()       noexcept { return jastrow_pade_; }
  [[nodiscard]] const JastrowPade& jastrow_pade() const noexcept { return jastrow_pade_; }

  [[nodiscard]]       SlaterPlaneWave& slater_plane_wave()       noexcept { return slater_plane_wave_; }
  [[nodiscard]] const SlaterPlaneWave& slater_plane_wave() const noexcept { return slater_plane_wave_; }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> j_derivatives() {
    return deriv_.view<idx(Derivatives::NUM), idx(ArrayIndex::DERIVATIVES)>();
  }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<const fp_t, idx(Derivatives::NUM)> j_derivatives() const {
    return deriv_.view<idx(Derivatives::NUM), idx(ArrayIndex::DERIVATIVES)>();
  }

  [[nodiscard]] bool jastrow_cache_valid() const noexcept { return jastrow_cache_valid_; }
  void set_jastrow_cache_valid(bool value) noexcept { jastrow_cache_valid_ = value; }

  [[nodiscard]] std::size_t steps_since_refresh() const noexcept { return steps_since_refresh_; }
  void set_steps_since_refresh(std::size_t value) noexcept { steps_since_refresh_ = value; }

  void evaluate_derivatives(Particles& particles) noexcept;
  void evaluate_derivatives(
    Particles& particles,
    bool move_accepted,
    std::size_t moved,
    xpu::array<fp_t, idx(Axis::NUM)> old_pos
  ) noexcept;

  fp_t evaluate_log_psi(Particles& particles) {
    return this->slater_plane_wave().log_abs_det(particles) + this->jastrow_pade().value(particles.pos());
  }
};
