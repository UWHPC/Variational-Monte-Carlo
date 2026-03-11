#pragma once

#include "../config/config.hpp"

#include <cstddef>

class BlockingAnalysis {
private:
    // For Welford's online algorithm:
    std::size_t num_blocks_; // completed block count
    double running_mean_;    // Welford running mean
    double running_m2_;      // Welford sum of squared deviations

    std::size_t block_size_;
    std::size_t in_block_;
    double block_sum_;

public:
    explicit BlockingAnalysis(std::size_t block_size);

    std::pair<double, double> mean_and_standard_error() const;
    void add(double local_energy);
    bool ready() const noexcept;
};