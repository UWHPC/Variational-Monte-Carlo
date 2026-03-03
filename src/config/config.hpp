#pragma once

#include <cstdint>
#include <cstddef>

struct Config {
    std::size_t numParticles;  // Number of particles
    double boxLength;          // Length of box (grid)
    std::size_t warmupSteps;   // Warm up 
    std::size_t measureSteps;  // M
    double stepSize;           // proposal half-width s
    uint64_t seed;             // Random seed
    std::size_t blockSize;     // Size of block
};
