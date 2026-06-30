#pragma once

#include "../utilities/aligned_soa.hpp"
#include "../utilities/macros.hpp"
#include "../utilities/ptr3d.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>

class Particles {
private:
  std::size_t size_;

  enum ArrayIndex : std::size_t {
    POS_X,
    POS_Y,
    POS_Z,
    GRAD_X,
    GRAD_Y,
    GRAD_Z,
    LAP_LOG_PSI,
    NUM_SUB_ARRAYS
  };
  AlignedSoA<double> data_;

public:
  explicit Particles(std::size_t num_particles)
  : size_{num_particles}
  , data_{num_particles, NUM_SUB_ARRAYS}
  { }

  [[nodiscard]] std::size_t size() const { return size_; }
  [[nodiscard]] std::size_t p_stride() const { return data_.stride(); }

  [[nodiscard]] double*       lap_log_psi()       noexcept { return data_[LAP_LOG_PSI]; }
  [[nodiscard]] double const* lap_log_psi() const noexcept { return data_[LAP_LOG_PSI]; }

  Ptr3D<      double> pos()       noexcept { return {data_[POS_X], data_[POS_Y], data_[POS_Z]}; }
  Ptr3D<const double> pos() const noexcept { return {data_[POS_X], data_[POS_Y], data_[POS_Z]}; }

  Ptr3D<      double> grad_log_psi()       noexcept { return {data_[GRAD_X], data_[GRAD_Y], data_[GRAD_Z]}; }
  Ptr3D<const double> grad_log_psi() const noexcept { return {data_[GRAD_X], data_[GRAD_Y], data_[GRAD_Z]}; }
};
