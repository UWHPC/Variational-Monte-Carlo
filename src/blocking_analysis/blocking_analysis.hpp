#pragma once

#include "../utilities/macros.hpp"
#include <xpu/math.hpp>

#include <cstddef>
#include <utility>

class BlockingAnalysis {
private:
  std::size_t num_blocks_;
  fp_t running_mean_;
  fp_t running_m2_;

  std::size_t block_size_;
  std::size_t in_block_;
  fp_t block_sum_;

public:
  explicit BlockingAnalysis(std::size_t block_size);

  [[nodiscard]] std::pair<fp_t, fp_t> mean_and_standard_error() const;
  void add(fp_t local_energy);
  [[nodiscard]] bool ready() const noexcept;
};