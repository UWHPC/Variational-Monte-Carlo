#pragma once

#include "../utilities/components.hpp"
#include "../particles/particles.hpp"
#include <xpu/buffer.hpp>
#include <xpu/soa.hpp>
#include <xpu/math.hpp>
#include <cstddef>

class EnergyTracker {
private:
  real_t box_length_;

  static constexpr real_t EWALD_RECIPROCAL_TOLERANCE{1.0e-6_r};
  real_t ewald_alpha_;
  real_t ewald_correction_;
  real_t ewald_background_;

  real_t V_recip_;
  real_t V_real_;

  std::size_t num_g_vectors_;

  enum class ArrayIndex : std::size_t {
    G_X,
    G_Y,
    G_Z,
    G_WEIGHTS,
    S_REAL,
    S_IMAG,
    NUM_ARRAYS
  };
  xpu::soa<real_t, idx(ArrayIndex::NUM_ARRAYS)> data_;
  mutable xpu::buffer<real_t> reduction_scratch_;

public:
  struct View {
    std::size_t num_g_vectors;
    real_t ewald_alpha;
    xpu::soa_view<real_t, idx(Axis::NUM)> g_vector;
    const real_t* g_weights;
    real_t* sum_real;
    real_t* sum_imag;
  };

  explicit EnergyTracker(real_t box_length, std::size_t num_particles);

  [[nodiscard]] std::size_t num_g_vectors() const noexcept {
    return num_g_vectors_;
  }

  void initialize_reciprocal_energy() noexcept;
  void initialize_real_energy(const Particles& particles) noexcept;

  void initialize_structure_factors(const Particles& particles) noexcept;

  void update_structure_factors(
    xpu::array<real_t, idx(Axis::NUM)> old_pos,
    xpu::array<real_t, idx(Axis::NUM)> new_pos
  ) noexcept;

  void update_real_energy(
    std::size_t moved,
    xpu::array<real_t, idx(Axis::NUM)> old_pos,
    const Particles& particles
  ) noexcept;

  real_t eval_total_energy(const Particles& particles) const noexcept {
    return kinetic_energy(particles) + potential_energy();
  }

  [[nodiscard]] CUDA_CALLABLE
  View view() noexcept {
    return {
      this->num_g_vectors(),
      ewald_alpha_,
      this->g_vector(),
      this->g_weights(),
      this->sum_real(),
      this->sum_imag()
    };
  }

  void accept_move(
    real_t real_energy_delta,
    real_t reciprocal_energy
  ) noexcept {
    V_real_ += real_energy_delta;
    V_recip_ = reciprocal_energy;
  }

private:
  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<real_t, idx(Axis::NUM)> g_vector() {
    return data_.view<idx(Axis::NUM), idx(ArrayIndex::G_X)>();
  }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<const real_t, idx(Axis::NUM)> g_vector() const {
    return data_.view<idx(Axis::NUM), idx(ArrayIndex::G_X)>();
  }

  [[nodiscard]] real_t*       g_weights()       noexcept { return data_[idx(ArrayIndex::G_WEIGHTS)]; }
  [[nodiscard]] real_t const* g_weights() const noexcept { return data_[idx(ArrayIndex::G_WEIGHTS)]; }

  [[nodiscard]] real_t*       sum_real()       noexcept { return data_[idx(ArrayIndex::S_REAL)]; }
  [[nodiscard]] real_t const* sum_real() const noexcept { return data_[idx(ArrayIndex::S_REAL)]; }

  [[nodiscard]] real_t*       sum_imag()       noexcept { return data_[idx(ArrayIndex::S_IMAG)]; }
  [[nodiscard]] real_t const* sum_imag() const noexcept { return data_[idx(ArrayIndex::S_IMAG)]; }

  [[nodiscard]] real_t* reduction_scratch() const noexcept {
    return reduction_scratch_.data();
  }

  real_t kinetic_energy(const Particles& particles) const noexcept;
  inline real_t potential_energy() const noexcept {
    return V_real_ + V_recip_ + ewald_correction_ + ewald_background_;
  }
};
