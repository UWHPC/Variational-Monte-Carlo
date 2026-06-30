#pragma once

#include "../particles/particles.hpp"
#include "../utilities/aligned_soa.hpp"
#include "../utilities/ptr3d.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

class EnergyTracker {
private:
  double box_length_;

  static constexpr double EWALD_RECIPROCAL_TOLERANCE{1.0e-6};
  double ewald_alpha_;
  double ewald_correction_;
  double ewald_background_;

  double V_recip_;
  double V_real_;

  std::size_t num_g_vectors_;

  enum ArrayIndex : std::size_t {
    G_X_,
    G_Y_,
    G_Z_,
    G_WEIGHTS_,
    S_REAL_,
    S_IMAG_,
    D_REAL_TEMP_,
    D_IMAG_TEMP_,
    NUM_ARRAYS
  };
  AlignedSoA<double> data_;

public:
  explicit EnergyTracker(double box_length, double num_particles);

  void initialize_reciprocal_energy() noexcept;
  void initialize_real_energy(const Particles& particles) noexcept;

  void initialize_structure_factors(const Particles& particles) noexcept;

  void update_structure_factors(
    double old_x,
    double old_y,
    double old_z,
    double new_x,
    double new_y,
    double new_z
  ) noexcept;

  void update_real_energy(
    std::size_t moved_idx,
    double old_x,
    double old_y,
    double old_z,
    const Particles& particles
  ) noexcept;

  double eval_total_energy(const Particles& particles) const noexcept;

private:
  Ptr3D<      double> g_vector()       noexcept { return {data_[G_X_], data_[G_Y_], data_[G_Z_]}; }
  Ptr3D<const double> g_vector() const noexcept { return {data_[G_X_], data_[G_Y_], data_[G_Z_]}; }

  [[nodiscard]] double*       g_weights()       noexcept { return data_[G_WEIGHTS_]; }
  [[nodiscard]] double const* g_weights() const noexcept { return data_[G_WEIGHTS_]; }

  [[nodiscard]] double*       sum_real()       noexcept { return data_[S_REAL_]; }
  [[nodiscard]] double const* sum_real() const noexcept { return data_[S_REAL_]; }

  [[nodiscard]] double*       sum_imag()       noexcept { return data_[S_IMAG_]; }
  [[nodiscard]] double const* sum_imag() const noexcept { return data_[S_IMAG_]; }

  [[nodiscard]] double*       d_real_temp()       noexcept { return data_[D_REAL_TEMP_]; }
  [[nodiscard]] double const* d_real_temp() const noexcept { return data_[D_REAL_TEMP_]; }

  [[nodiscard]] double*       d_imag_temp()       noexcept { return data_[D_IMAG_TEMP_]; }
  [[nodiscard]] double const* d_imag_temp() const noexcept { return data_[D_IMAG_TEMP_]; }

  double kinetic_energy(const Particles& particles) const noexcept;
  double potential_energy() const noexcept;
};