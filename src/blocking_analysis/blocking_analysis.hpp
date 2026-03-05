#pragma once
#include <vector>
#include <numeric>
#include "../config/config.hpp"

class BlockingAnalysis {
private:
    std::size_t blockSize;
    std::size_t K;
    double blockSum{};
    std::size_t inBlock{};
    std::vector<double> blockMeans;
public:
    explicit BlockingAnalysis(std::size_t blockSize);
    std::pair<double, double> BlockingAnalysis::meanAndStandardError() const;
    void add(double x);
    bool ready() const noexcept;

    


};