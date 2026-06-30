#pragma once

#include <cstddef>
#include <utility>

class BlockingAnalysis {
private:
  std::size_t num_blocks_;
  double running_mean_;
  double running_m2_;

  std::size_t block_size_;
  std::size_t in_block_;
  double block_sum_;

public:
  explicit BlockingAnalysis(std::size_t block_size);

  [[nodiscard]] std::pair<double, double> mean_and_standard_error() const;
  void add(double local_energy);
  [[nodiscard]] bool ready() const noexcept;
};