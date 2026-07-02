#pragma once

#include "../utilities/aligned_soa.cuh"
#include "../utilities/macros.cuh"
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
  AlignedSoA<real_t> data_;

public:
  explicit Particles(std::size_t num_particles)
  : size_{num_particles}
  , data_{num_particles, NUM_SUB_ARRAYS}
  { }

  [[nodiscard]] std::size_t size() const { return size_; }
  [[nodiscard]] std::size_t p_stride() const { return data_.stride(); }

  [[nodiscard]] real_t*       lap_log_psi()       noexcept { return data_[LAP_LOG_PSI]; }
  [[nodiscard]] real_t const* lap_log_psi() const noexcept { return data_[LAP_LOG_PSI]; }

  Ptr3D<      real_t> pos()       noexcept { return {data_[POS_X], data_[POS_Y], data_[POS_Z]}; }
  Ptr3D<const real_t> pos() const noexcept { return {data_[POS_X], data_[POS_Y], data_[POS_Z]}; }

  Ptr3D<      real_t> grad_log_psi()       noexcept { return {data_[GRAD_X], data_[GRAD_Y], data_[GRAD_Z]}; }
  Ptr3D<const real_t> grad_log_psi() const noexcept { return {data_[GRAD_X], data_[GRAD_Y], data_[GRAD_Z]}; }
};
