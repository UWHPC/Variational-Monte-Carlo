#include <exception>
#include <iostream>
#include <stdexcept>

#include "particles/particles.hpp"

struct Config {
  int N;
  double L;
  int warmup_steps;
  int measure_steps;
  double step_size;
  uint64_t seed;
  int block_size;
};

int main(int argc, char* argv[]) {
  return 0;
}
