#pragma once

#if defined(__CUDACC__)
#include <curand_kernel.h>
#endif

#include <random>
#include <memory>

#include "../config/config.hpp"
#include "../utilities/macros.cuh"

class WalkerRNG {
private:
#if defined(__CUDA_ARCH__)
  curandStatePhilox4_32_10_t rng_;
#else
  std::mt19937_64 rng_;
  std::uniform_real_distribution<real_t> uniform01_{0.0_r, 1.0_r};
  std::uniform_real_distribution<real_t> proposal_;
  std::uniform_int_distribution<std::size_t> pick_particle_;
#endif

  std::size_t N_;
  real_t step_size_;

public:
  #if defined(__CUDA_ARCH__)
  __device__ 
  WalkerRNG() = default;

  __device__ 
  void init(const Config& config, uint64_t walker_id, uint64_t offset = 0) {
    step_size_ = config.step_size;
    curand_init(config.master_seed, walker_id, offset, &rng_);
    N_ = config.num_particles;

  }
  #else
  WalkerRNG() = default;

  void init(const Config& config, uint64_t walker_id, uint64_t offset = 0) {
    step_size_ = config.step_size;
    rng_.seed(config.master_seed + walker_id);
    proposal_ = std::uniform_real_distribution<real_t>(-config.step_size, config.step_size);
    pick_particle_ = std::uniform_int_distribution<std::size_t>(0, config.num_particles - 1);
  }
  #endif

  #if defined(__CUDA_ARCH__) 
  [[nodiscard]] __device__ 
  real_t rand_uniform() {
    #ifdef FP_64
      return 1.0_r - curand_uniform_double(&rng_);
    #else
      return 1.0_r - curand_uniform(&rng_);
    #endif
  }

  [[nodiscard]] __device__ 
  real_t rand_proposal() {
    #ifdef FP_64
      return (curand_uniform_double(&rng_) * 2.0_r - 1.0_r) * step_size_;
    #else 
      return (curand_uniform(&rng_) * 2.0_r - 1.0_r) * step_size_;
    #endif
  }

  [[nodiscard]] __device__ 
  std::size_t rand_particle() {
    return curand(&rng_) % N_;
  }

  __device__ 
  void change_step_size(real_t step_size) {
    step_size_ = step_size;
  }
  #else
  [[nodiscard]] real_t rand_uniform() {
    return uniform01_(rng_);
  }

  [[nodiscard]] real_t rand_proposal() {
    return proposal_(rng_);
  }

  [[nodiscard]] std::size_t rand_particle() {
    return pick_particle_(rng_);
  }

  void change_step_size(real_t step_size) {
    step_size_ = step_size;
    proposal_.param(std::uniform_real_distribution<real_t>::param_type(-step_size, step_size));
  }
  #endif
};

