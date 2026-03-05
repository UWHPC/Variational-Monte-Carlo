#include "blocking_analysis.hpp"

BlockingAnalysis::BlockingAnalysis(std::size_t block_size)
    : block_size_{block_size}, in_block_{}, block_sum_{}, block_means_{} {}

std::pair<double, double> BlockingAnalysis::mean_and_standard_error() const {
    if (!ready()) {
        // The document says if K < 2, report mean only or flag it [cite: 35]
        throw std::runtime_error("Not enough blocks to calculate Standard Error.");
    }
    const std::size_t block_size{block_means_get().size()};
    const double K{static_cast<double>(block_size)};

    double overall_mean{};
    auto& block_means{block_means_get()};

    // calculating average
    const double inv_K{1.0 / K};
    overall_mean = std::reduce(block_means.begin(), block_means.end(), 0.0) * inv_K;

    // 2. Calculate the variance between the blocks (Eq. 39)
    double variance_sum{};

    for (std::size_t i = 0; i < block_size; ++i) {
        variance_sum += (block_means[i] - overall_mean) * (block_means[i] - overall_mean);
    }
    const double denom{1.0 / (K - 1.0)};
    const double s_sq{variance_sum * denom};

    // 3. Calculate the Standard Error (Eq. 40)
    const double standard_error{std::sqrt(s_sq * inv_K)};

    return {overall_mean, standard_error};
}

void BlockingAnalysis::add(double x) {
    block_sum_ += x;
    in_block_++;

    // When the block is full, calculate its average and save it
    if (in_block_ == block_size_) {
        block_means_.push_back(block_sum_ / block_size_); // Eq. 37

        // Reset for the next block
        block_sum_ = 0;
        in_block_ = 0;
    }
}

bool BlockingAnalysis::ready() const noexcept { return block_means_.size() >= 2; }
