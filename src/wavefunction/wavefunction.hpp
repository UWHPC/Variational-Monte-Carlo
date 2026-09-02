#pragma once

#include "../jastrow_pade/jastrow_pade.hpp"
#include "../particles/particles.hpp"
#include "../slater_plane_wave/slater_plane_wave.hpp"
#include <xpu/buffer.hpp>
#include <xpu/soa.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>

class WaveFunction {
private:
  JastrowPade jastrow_pade_;
  SlaterPlaneWave slater_plane_wave_;

  xpu::buffer<std::uint8_t> jastrow_cache_valid_;
  xpu::buffer<std::size_t> steps_since_refresh_;

  enum class ArrayIndex : std::size_t {
    DERIVATIVES,
    NUM_ARRAYS = enum_index(ArrayIndex::DERIVATIVES, Derivatives::NUM)
  };
  xpu::soa_batch<fp_t, idx(ArrayIndex::NUM_ARRAYS)> deriv_;

public:
  struct View {
    JastrowPade::View jastrow{};
    SlaterPlaneWave::View slater{};
    xpu::soa_view<fp_t, idx(Derivatives::NUM)> jastrow_derivatives{nullptr, 0uz};
    std::uint8_t* jastrow_cache_valid{};
    std::size_t* steps_since_refresh{};
  };

  explicit WaveFunction(
    const Particles& particles,
    fp_t box_length,
    fp_t a = 0.25_fp,
    fp_t b = 1.0_fp
  )
    : jastrow_pade_{box_length, a, b}
    , slater_plane_wave_{particles, box_length}
    , jastrow_cache_valid_{particles.walker_count()}
    , steps_since_refresh_{particles.walker_count()}
    , deriv_{particles.walker_count(), particles.count()}
  {
    xpu::zero_n(jastrow_cache_valid_.data(), particles.walker_count());
    xpu::zero_n(steps_since_refresh_.data(), particles.walker_count());
  }

  [[nodiscard]]       JastrowPade& jastrow_pade()       noexcept { return jastrow_pade_; }
  [[nodiscard]] const JastrowPade& jastrow_pade() const noexcept { return jastrow_pade_; }

  [[nodiscard]]       SlaterPlaneWave& slater_plane_wave()       noexcept { return slater_plane_wave_; }
  [[nodiscard]] const SlaterPlaneWave& slater_plane_wave() const noexcept { return slater_plane_wave_; }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> j_derivatives(std::size_t walker = 0uz) {
    return deriv_.view<idx(Derivatives::NUM), idx(ArrayIndex::DERIVATIVES)>(walker);
  }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<const fp_t, idx(Derivatives::NUM)> j_derivatives(std::size_t walker = 0uz) const {
    return deriv_.view<idx(Derivatives::NUM), idx(ArrayIndex::DERIVATIVES)>(walker);
  }

  [[nodiscard]] CUDA_CALLABLE
  View view(std::size_t walker = 0uz) noexcept {
    assert(walker < deriv_.batch_count());
    return {
      this->jastrow_pade().view(),
      this->slater_plane_wave().view(walker),
      this->j_derivatives(walker),
      jastrow_cache_valid_.data() + walker,
      steps_since_refresh_.data() + walker
    };
  }

  [[nodiscard]] bool jastrow_cache_valid(std::size_t walker = 0uz) const noexcept {
    assert(walker < deriv_.batch_count());
    return jastrow_cache_valid_.data()[walker] != 0u;
  }
  void set_jastrow_cache_valid(bool value, std::size_t walker = 0uz) noexcept {
    assert(walker < deriv_.batch_count());
    jastrow_cache_valid_.data()[walker] = static_cast<std::uint8_t>(value);
  }

  [[nodiscard]] std::size_t steps_since_refresh(std::size_t walker = 0uz) const noexcept {
    assert(walker < deriv_.batch_count());
    return steps_since_refresh_.data()[walker];
  }
  void set_steps_since_refresh(std::size_t value, std::size_t walker = 0uz) noexcept {
    assert(walker < deriv_.batch_count());
    steps_since_refresh_.data()[walker] = value;
  }

  void evaluate_derivatives(Particles::View particles, std::size_t walker = 0uz) noexcept;
  void evaluate_derivatives(
    Particles::View particles,
    bool move_accepted,
    std::size_t moved,
    xpu::array<fp_t, idx(Axis::NUM)> old_pos,
    std::size_t walker = 0uz
  ) noexcept;

  fp_t evaluate_log_psi(Particles::View particles, std::size_t walker = 0uz) {
    return this->slater_plane_wave().log_abs_det(particles, walker)
      + this->jastrow_pade().value(particles);
  }
};
