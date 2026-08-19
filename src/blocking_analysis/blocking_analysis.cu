#include "blocking_analysis.hpp"

#include <stdexcept>

BlockingAnalysis::BlockingAnalysis(std::size_t block_size)
: num_blocks_{}
, running_mean_{}
, running_m2_{}
, block_size_{block_size}
, in_block_{}
, block_sum_{}
{ }

std::pair<fp_t, fp_t> BlockingAnalysis::mean_and_standard_error() const {
  if (num_blocks_ < 2) {
    throw std::runtime_error("Not enough blocks");
  }
  const fp_t variance{running_m2_ / static_cast<fp_t>(num_blocks_ - 1)};
  const fp_t standard_error{xpu::sqrt(variance / static_cast<fp_t>(num_blocks_))};

  return {running_mean_, standard_error};
}

void BlockingAnalysis::add(fp_t local_energy) {
  block_sum_ += local_energy;
  ++in_block_;

  if (in_block_ == block_size_) {
    const fp_t block_mean{block_sum_ / static_cast<fp_t>(block_size_)};

    ++num_blocks_;

    const fp_t delta{block_mean - running_mean_};
    running_mean_ += delta / static_cast<fp_t>(num_blocks_);
    running_m2_ += delta * (block_mean - running_mean_);

    block_sum_ = 0.0_fp;
    in_block_ = 0U;
  }
}

bool BlockingAnalysis::ready() const noexcept { return num_blocks_ >= 2; }