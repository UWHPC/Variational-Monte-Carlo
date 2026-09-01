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
  xpu::soa_batch<fp_t, idx(ArrayIndex::NUM_SUB_ARRAYS)> data_;

public:
  explicit Particles(
    std::size_t num_particles,
    std::size_t num_walkers = 1uz
  )
    : data_{num_walkers, num_particles}
  { }

  [[nodiscard]]
  std::size_t count() const {
    return data_.element_count();
  }

  [[nodiscard]]
  std::size_t walker_count() const {
    return data_.batch_count();
  }

  [[nodiscard]]
  std::size_t stride() const {
    return data_.array_stride();
  }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos(std::size_t walker = 0uz) {
    return data_.view<idx(Axis::NUM), idx(ArrayIndex::POS)>(walker);
  }
  
  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<const fp_t, idx(Axis::NUM)> pos(std::size_t walker = 0uz) const {
    return data_.view<idx(Axis::NUM), idx(ArrayIndex::POS)>(walker);
  }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> derivatives(std::size_t walker = 0uz) {
    return data_.view<idx(Derivatives::NUM), idx(ArrayIndex::DERIVATIVES)>(walker);
  }
  
  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<const fp_t, idx(Derivatives::NUM)> derivatives(std::size_t walker = 0uz) const {
    return data_.view<idx(Derivatives::NUM), idx(ArrayIndex::DERIVATIVES)>(walker);
  }
};
