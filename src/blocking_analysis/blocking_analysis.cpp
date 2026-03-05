#include "blocking_analysis.hpp"

BlockingAnalysis::BlockingAnalysis(std::size_t block_size)
    : block_size_{block_size}, in_block_{}, block_sum_{}, block_means_{} {}

std::pair<double, double> BlockingAnalysis::mean_and_standard_error() const {
    if (!ready()) {
        // The document says if K < 2, report mean only or flag it [cite: 35]
        throw std::runtime_error("Not enough blocks to calculate Standard Error.");
    }

    const double K{static_cast<double>(block_means_get().size())};
    double overall_mean{};

    // calculating average
    overall_mean = std::accumulate(block_means_get().begin(), block_means_get().end(), 0.0) / K;

    // 2. Calculate the variance between the blocks (Eq. 39)
    double variance_sum{};
    for (std::size_t i = 0; i < block_means_get().size(); ++i) {
        variance_sum += (block_means_get()[i] - overall_mean) * (block_means_get()[i] - overall_mean);
    }
    double s_sq = variance_sum / (K - 1.0);

    // 3. Calculate the Standard Error (Eq. 40)
    double standard_error{std::sqrt(s_sq / K)};

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
