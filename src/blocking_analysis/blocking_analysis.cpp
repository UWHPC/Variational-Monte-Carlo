#include "blocking_analysis.hpp"

BlockingAnalysis::BlockingAnalysis(std::size_t block_size)
    : num_blocks_{0}, running_mean_{0.0}, running_m2_{0.0}, block_size_{block_size}, in_block_{0}, block_sum_{0.0} {}

std::pair<double, double> BlockingAnalysis::mean_and_standard_error() const {
    if (num_blocks_ < 2) {
        throw std::runtime_error("Not enough blocks");
    }
    const double variance{running_m2_ / static_cast<double>(num_blocks_ - 1)};
    const double standard_error{std::sqrt(variance / static_cast<double>(num_blocks_))};

    return {running_mean_, standard_error};
}

void BlockingAnalysis::add(double local_energy) {
    block_sum_ += local_energy;
    ++in_block_;

    // When the block is full, calculate its average and save it
    if (in_block_ == block_size_) {
        // Welford's online algorithm:
        const double block_mean{block_sum_ / static_cast<double>(block_size_)};
        ++num_blocks_;
        const double delta{block_mean - running_mean_};

        running_mean_ += delta / static_cast<double>(num_blocks_);
        running_m2_ += delta * (block_mean - running_mean_);

        // Reset for the next block
        block_sum_ = 0;
        in_block_ = 0;
    }
}

bool BlockingAnalysis::ready() const noexcept { return num_blocks_ >= 2; }