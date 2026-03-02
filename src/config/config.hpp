#include <cstdint>

struct Config {
  std::size_t numParticles; // Number of particles
  double boxLength;         // Length of box (grid)
  std::size_t warmupSteps;  // 
  std::size_t measureSteps; //
  double stepSize;          // 
  uint64_t seed;            // Seed
  std::size_t blockSize;    // Size of block
};