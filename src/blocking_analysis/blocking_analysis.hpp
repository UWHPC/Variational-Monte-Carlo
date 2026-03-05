#pragma once

#include "../config/config.hpp"

class BlockingAnalysis {
private:
    std::size_t block_size_;
    std::size_t K_;

public:
    explicit BlockingAnalysis(std::size_t block_size);
    void add(double x);
    bool ready() const noexcept;
};
