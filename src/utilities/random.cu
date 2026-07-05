#pragma once

#if defined(__CUDACC__)
#include <curand_kernel.h>
#else
#include <random>
#endif

#include <memory>

#include "../config/config.hpp"
#include "../utilities/macros.cuh"

class WalkerRNG {
private:
#if defined(__CUDACC__)
  curandStatePhilox4_32_10_t rng_;
#else
  std::mt19937_64 rng_;
  std::uniform_real_distribution<real_t> uniform01_{0.0_r, 1.0_r};
  std::unique_ptr<std::uniform_real_distribution<real_t>> proposal_;
  std::unique_ptr<std::uniform_int_distribution<std::size_t>> pick_particle_;
#endif

  std::size_t N_;
  real_t step_size_;

public:
  CUDA_CALLABLE WalkerRNG() {}

  CUDA_CALLABLE void init(const Config& config, uint64_t walker_id, uint64_t offset = 0) {
    step_size_ = config.step_size;
  #if defined(__CUDACC__)
    curand_init(config.master_seed, walker_id, offset, &rng_);
    N_ = config.num_particles;
  #else 
    rng_.seed(config.master_seed + walker_id);
    *proposal_ = std::uniform_real_distribution<real_t>(-config.step_size, config.step_size);
    *pick_particle_ = std::uniform_int_distribution<std::size_t>(0, config.num_particles - 1);
  #endif
  }

  [[nodiscard]] CUDA_CALLABLE real_t rand_uniform() { 
  #if defined(__CUDACC__)
    #ifdef FP_64
      return 1.0_r - curand_uniform_double(&rng_);
    #else
      return 1.0_r - curand_uniform(&rng_);
    #endif
  #else
    return uniform01_(rng_); 
  #endif
  }

  [[nodiscard]] CUDA_CALLABLE real_t rand_proposal() { 
  #if defined(__CUDACC__)
    #ifdef FP_64
      return (curand_uniform_double(&rng_) * 2.0_r - 1.0_r) * step_size_; 
    #else
      return (curand_uniform(&rng_) * 2.0_r - 1.0_r) * step_size_; 
    #endif
  #else    
    return (*proposal_)(rng_); 
  #endif
  }

  [[nodiscard]] CUDA_CALLABLE std::size_t rand_particle() { 
  #if defined(__CUDACC__)
    return curand(&rng_) % N_;
  #else
    return (*pick_particle_)(rng_); 
  #endif
  }

  CUDA_CALLABLE void change_step_size(real_t step_size) {
    step_size_ = step_size;
  #if defined(__CUDACC__)
    return;
  #else
    (*proposal_).param(std::uniform_real_distribution<real_t>::param_type(-step_size, step_size));
  #endif
  }
};

