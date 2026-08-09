#pragma once

#include "../utilities/macros.hpp"
#include <xpu/math.hpp>

#include <cstddef>
#include <utility>

class BlockingAnalysis {
private:
  std::size_t num_blocks_;
  real_t running_mean_;
  real_t running_m2_;

  std::size_t block_size_;
  std::size_t in_block_;
  real_t block_sum_;

public:
  explicit BlockingAnalysis(std::size_t block_size);

  [[nodiscard]] std::pair<real_t, real_t> mean_and_standard_error() const;
  void add(real_t local_energy);
  [[nodiscard]] bool ready() const noexcept;
};