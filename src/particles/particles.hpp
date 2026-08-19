#pragma once

#include <xpu/soa.hpp>
#include "../utilities/macros.hpp"
#include "../utilities/components.hpp"

#include <cstdlib>
#include <cstring>

class Particles {
private:
  enum class ArrayIndex {
    POS,
    DERIVATIVES = enum_index(Axis::NUM, POS),
    NUM_SUB_ARRAYS = enum_index(Derivatives::NUM, DERIVATIVES)
  };
  xpu::soa<fp_t, idx(ArrayIndex::NUM_SUB_ARRAYS)> data_;

public:
  explicit Particles(std::size_t num_particles)
    : data_{num_particles}
  { }

  [[nodiscard]]
  std::size_t count() const {
    return data_.view().count();
  }

  [[nodiscard]]
  std::size_t stride() const {
    return data_.view().stride();
  }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos() {
    return data_.view<idx(Axis::NUM), idx(ArrayIndex::POS)>();
  }
  
  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<const fp_t, idx(Axis::NUM)> pos() const {
    return data_.view<idx(Axis::NUM), idx(ArrayIndex::POS)>();
  }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> derivatives() {
    return data_.view<idx(Derivatives::NUM), idx(ArrayIndex::DERIVATIVES)>();
  }
  
  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<const fp_t, idx(Derivatives::NUM)> derivatives() const {
    return data_.view<idx(Derivatives::NUM), idx(ArrayIndex::DERIVATIVES)>();
  }
};
