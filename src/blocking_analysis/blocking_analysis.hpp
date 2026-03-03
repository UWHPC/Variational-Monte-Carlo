#pragma once

#include "../config/config.hpp"

class BlockingAnalysis {
private:
    std::size_t blockSize_;
    std::size_t K_;
public:
    explicit BlockingAnalysis(std::size_t blockSize);
    void add(double x);
    bool ready() const noexcept;
};
