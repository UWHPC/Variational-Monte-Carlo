#pragma once

#include "../config/config.hpp"

#include <numeric>
#include <vector>

class BlockingAnalysis {
private:
    std::size_t block_size_;
    std::size_t in_block_;

    double block_sum_;
    std::vector<double> block_means_;

    [[nodiscard]] std::size_t in_block_get() const { return in_block_; }
    [[nodiscard]] double block_sum_get() const { return block_sum_; }
    [[nodiscard]] std::size_t block_size_get() const { return block_size_; }

    [[nodiscard]] std::vector<double>& block_means_get() { return block_means_; }
    [[nodiscard]] const std::vector<double>& block_means_get() const { return block_means_; }

    [[nodiscard]] double& block_sum_set() { return block_sum_; }
    [[nodiscard]] std::size_t& in_block_set() { return in_block_; }

public:
    explicit BlockingAnalysis(std::size_t block_size);

    std::pair<double, double> mean_and_standard_error() const;
    void add(double local_energy);
    bool ready() const noexcept;
};