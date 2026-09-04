#pragma once

#include "../utilities/components.hpp"
#include "../particles/particles.hpp"
#include <xpu/buffer.hpp>
#include <xpu/soa.hpp>
#include <xpu/math.hpp>
#include <cstddef>

class EnergyTracker {
private:
  fp_t box_length_;

  static constexpr fp_t EWALD_RECIPROCAL_TOLERANCE{1.0e-6_fp};
  fp_t ewald_alpha_;
  fp_t ewald_correction_;
  fp_t ewald_background_;

  std::size_t num_g_vectors_;
  std::size_t num_walkers_;

  enum class SharedArray : std::size_t {
    G_X,
    G_Y,
    G_Z,
    G_WEIGHTS,
    NUM_ARRAYS
  };
  enum class WalkerArray : std::size_t {
    S_REAL,
    S_IMAG,
    NUM_ARRAYS
  };
  enum class WalkerScalar : std::size_t {
    V_REAL,
    V_RECIP,
    NUM_ARRAYS
  };

  xpu::soa<fp_t, idx(SharedArray::NUM_ARRAYS)> shared_data_;
  xpu::soa_batch<fp_t, idx(WalkerArray::NUM_ARRAYS)> walker_data_;
  xpu::soa<fp_t, idx(WalkerScalar::NUM_ARRAYS)> walker_scalars_;
  mutable xpu::buffer<fp_t> reduction_scratch_;

public:
  struct View {
    fp_t box_length{};
    std::size_t num_g_vectors{};
    fp_t ewald_alpha{};
    fp_t ewald_correction{};
    fp_t ewald_background{};
    xpu::soa_view<fp_t, idx(Axis::NUM)> g_vector{nullptr, 0uz};
    const fp_t* g_weights{};
    fp_t* sum_real{};
    fp_t* sum_imag{};
    fp_t* real_energy{};
    fp_t* reciprocal_energy{};
    fp_t* reduction_scratch{};
  };

  explicit EnergyTracker(
    fp_t box_length,
    std::size_t num_particles,
    std::size_t num_walkers = 1uz
  );
  explicit EnergyTracker(fp_t box_length, const Particles& particles)
    : EnergyTracker{
        box_length,
        particles.count(),
        particles.walker_count()
      }
  { }

  [[nodiscard]] std::size_t num_g_vectors() const noexcept {
    return num_g_vectors_;
  }
  [[nodiscard]] std::size_t walker_count() const noexcept {
    return num_walkers_;
  }

  void initialize_reciprocal_energy(std::size_t walker = 0uz) noexcept;
  void initialize_real_energy(
    Particles::View particles,
    std::size_t walker = 0uz
  ) noexcept;

  void initialize_structure_factors(
    Particles::View particles,
    std::size_t walker = 0uz
  ) noexcept;

  void update_structure_factors(
    xpu::array<fp_t, idx(Axis::NUM)> old_pos,
    xpu::array<fp_t, idx(Axis::NUM)> new_pos,
    std::size_t walker = 0uz
  ) noexcept;

  void update_real_energy(
    std::size_t moved,
    xpu::array<fp_t, idx(Axis::NUM)> old_pos,
    Particles::View particles,
    std::size_t walker = 0uz
  ) noexcept;

  fp_t eval_total_energy(
    Particles::View particles,
    std::size_t walker = 0uz
  ) noexcept {
    return kinetic_energy(particles, walker) + potential_energy(walker);
  }

  [[nodiscard]]
  View view(std::size_t walker = 0uz) noexcept {
    return {
      box_length_,
      this->num_g_vectors(),
      ewald_alpha_,
      ewald_correction_,
      ewald_background_,
      this->g_vector(),
      this->g_weights(),
      this->sum_real(walker),
      this->sum_imag(walker),
      &this->real_energy(walker),
      &this->reciprocal_energy_value(walker),
      this->reduction_scratch(walker)
    };
  }

  void accept_move(
    fp_t real_energy_delta,
    fp_t reciprocal_energy,
    std::size_t walker = 0uz
  ) noexcept {
    real_energy(walker) += real_energy_delta;
    reciprocal_energy_value(walker) = reciprocal_energy;
  }

private:
  [[nodiscard]]
  xpu::soa_view<fp_t, idx(Axis::NUM)> g_vector() {
    return shared_data_.view<idx(Axis::NUM), idx(SharedArray::G_X)>();
  }

  [[nodiscard]]
  xpu::soa_view<const fp_t, idx(Axis::NUM)> g_vector() const {
    return shared_data_.view<idx(Axis::NUM), idx(SharedArray::G_X)>();
  }

  [[nodiscard]] fp_t* g_weights() noexcept {
    return shared_data_[idx(SharedArray::G_WEIGHTS)];
  }
  [[nodiscard]] fp_t const* g_weights() const noexcept {
    return shared_data_[idx(SharedArray::G_WEIGHTS)];
  }

  [[nodiscard]] fp_t* sum_real(std::size_t walker) noexcept {
    return walker_data_.view<1uz, idx(WalkerArray::S_REAL)>(walker)[0uz];
  }
  [[nodiscard]] fp_t const* sum_real(std::size_t walker) const noexcept {
    return walker_data_.view<1uz, idx(WalkerArray::S_REAL)>(walker)[0uz];
  }

  [[nodiscard]] fp_t* sum_imag(std::size_t walker) noexcept {
    return walker_data_.view<1uz, idx(WalkerArray::S_IMAG)>(walker)[0uz];
  }
  [[nodiscard]] fp_t const* sum_imag(std::size_t walker) const noexcept {
    return walker_data_.view<1uz, idx(WalkerArray::S_IMAG)>(walker)[0uz];
  }

  [[nodiscard]] fp_t& real_energy(std::size_t walker) noexcept {
    return walker_scalars_[idx(WalkerScalar::V_REAL)][walker];
  }
  [[nodiscard]] fp_t real_energy(std::size_t walker) const noexcept {
    return walker_scalars_[idx(WalkerScalar::V_REAL)][walker];
  }
  [[nodiscard]] fp_t& reciprocal_energy_value(std::size_t walker) noexcept {
    return walker_scalars_[idx(WalkerScalar::V_RECIP)][walker];
  }
  [[nodiscard]] fp_t reciprocal_energy_value(std::size_t walker) const noexcept {
    return walker_scalars_[idx(WalkerScalar::V_RECIP)][walker];
  }

  [[nodiscard]] fp_t* reduction_scratch(std::size_t walker) const noexcept {
    return reduction_scratch_.data() + walker;
  }

  fp_t kinetic_energy(
    Particles::View particles,
    std::size_t walker
  ) noexcept;
  inline fp_t potential_energy(std::size_t walker) const noexcept {
    return
      real_energy(walker) +
      reciprocal_energy_value(walker) +
      ewald_correction_ +
      ewald_background_;
  }
};
