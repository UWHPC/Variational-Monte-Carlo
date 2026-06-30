#pragma once

#include "../jastrow_pade/jastrow_pade.hpp"
#include "../particles/particles.hpp"
#include "../slater_plane_wave/slater_plane_wave.hpp"
#include "../utilities/aligned_soa.hpp"

#include <cstddef>
#include <vector>

class WaveFunction {
private:
  JastrowPade jastrow_pade_;
  SlaterPlaneWave slater_plane_wave_;

  bool jastrow_cache_valid_;
  std::size_t steps_since_refresh_;

  enum ArrayIndex : std::size_t { GRAD_X, GRAD_Y, GRAD_Z, LAP, NUM_ARRAYS };
  AlignedSoA<double> deriv_;

public:
  explicit WaveFunction(
    const Particles& particles,
    double box_length,
    double a = 0.25,
    double b = 1.0
  )
  : jastrow_pade_{box_length, a, b}
  , slater_plane_wave_{particles, box_length}
  , jastrow_cache_valid_{}
  , steps_since_refresh_{}
  , deriv_{particles.size(), NUM_ARRAYS}
  { }

  [[nodiscard]]       JastrowPade&     jastrow_pade()             noexcept { return jastrow_pade_; }
  [[nodiscard]] const JastrowPade&     jastrow_pade()       const noexcept { return jastrow_pade_; }

  [[nodiscard]]       SlaterPlaneWave& slater_plane_wave()        noexcept { return slater_plane_wave_; }
  [[nodiscard]] const SlaterPlaneWave& slater_plane_wave()  const noexcept { return slater_plane_wave_; }

  Ptr3D<      double> j_grad()       noexcept { return {deriv_[GRAD_X], deriv_[GRAD_Y], deriv_[GRAD_Z]}; }
  Ptr3D<const double> j_grad() const noexcept { return {deriv_[GRAD_X], deriv_[GRAD_Y], deriv_[GRAD_Z]}; }

  [[nodiscard]] double*       j_lap()       noexcept { return deriv_[LAP]; }
  [[nodiscard]] double const* j_lap() const noexcept { return deriv_[LAP]; }

  [[nodiscard]] bool jastrow_cache_valid() const noexcept { return jastrow_cache_valid_; }
  void set_jastrow_cache_valid(bool value) noexcept { jastrow_cache_valid_ = value; }

  [[nodiscard]] std::size_t steps_since_refresh() const noexcept { return steps_since_refresh_; }
  void set_steps_since_refresh(std::size_t value) noexcept { steps_since_refresh_ = value; }

  void evaluate_derivatives(Particles& particles) noexcept;
  void evaluate_derivatives(
    Particles& particles,
    bool move_accepted,
    std::size_t moved,
    double old_x, double old_y, double old_z
  ) noexcept;

  double evaluate_log_psi(const Particles& particles);
};
