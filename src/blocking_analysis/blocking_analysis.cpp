#include "blocking_analysis.hpp"


std::pair<double, double> BlockingAnalysis::meanAndStandardError() const {
        if (!ready()) {
            // The document says if K < 2, report mean only or flag it [cite: 35]
            throw std::runtime_error("Not enough blocks to calculate Standard Error.");
        }

        const double K{static_cast<double>(blockMeans.size())};
        double overallMean{};

        // calculating average
        overallMean = std::accumulate(blockMeans.begin(), blockMeans.end(), 0.0) / K;

        // 2. Calculate the variance between the blocks (Eq. 39) 
        double varianceSum{};
        for (double i = 0; i < K; i++ ) {
            varianceSum += (blockMeans[i] - overallMean) * (blockMeans[i] - overallMean);
        }
        double s_sq = varianceSum / (K - 1.0);

        // 3. Calculate the Standard Error (Eq. 40) 
        double standardError{std::sqrt(s_sq / K)};

        return {overallMean, standardError};
    }

    void BlockingAnalysis::add(double x){
        blockSum += x;
        inBlock++;

        // When the block is full, calculate its average and save it
        if (inBlock == blockSize) {
            blockMeans.push_back(blockSum / blockSize); // Eq. 37 
            
            // Reset for the next block
            blockSum = 0;
            inBlock = 0;
        }
    }

    bool BlockingAnalysis::ready() const noexcept {
        return blockMeans.size() >= 2; 
    }
