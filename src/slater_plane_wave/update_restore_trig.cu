#include "slater_plane_wave.cuh"

#if defined(__CUDACC__)
#include <cstddef>
namespace {
  __global__ void cudaUpdateTrigRow(
    std::size_t particle, std::size_t num_k, std::size_t row_stride,
    const real_t* RESTRICT p_x, const real_t* RESTRICT p_y, const real_t* RESTRICT p_z,
    const real_t* RESTRICT k_x, const real_t* RESTRICT k_y, const real_t* RESTRICT k_z,
    real_t* RESTRICT s_cache, real_t* RESTRICT c_cache
  ) {
    const std::size_t k{blockIdx.x * blockDim.x + threadIdx.x};
    if (k >= num_k) { return; }

    const std::size_t offset{particle * row_stride};

    const real_t dot{
      k_x[k] * p_x[particle] +
      k_y[k] * p_y[particle] +
      k_z[k] * p_z[particle]
    };

    vmc::sincos(dot, &s_cache[offset + k], &c_cache[offset + k]);
  }
}
#else
#include "particles/particles.cuh"
#include "utilities/aligned_soa.cuh"
#include <cstring>
#endif


void SlaterPlaneWave::restore_trig_row(std::size_t particle) {
  const std::size_t num_k{num_unique_k()};
  const std::size_t ROW_STRIDE{trig_row_stride()};
  const std::size_t offset{particle * ROW_STRIDE};

#if defined(__CUDACC__)
  CUDA_CHECK(cudaMemcpy(sin_cache() + offset, trig_scratch_[SIN_SAVED], num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy(cos_cache() + offset, trig_scratch_[COS_SAVED], num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
#else
  std::memcpy(sin_cache() + offset, trig_scratch_[SIN_SAVED], num_k * sizeof(real_t));
  std::memcpy(cos_cache() + offset, trig_scratch_[COS_SAVED], num_k * sizeof(real_t));
#endif
}

void SlaterPlaneWave::update_trig_cache(
  std::size_t particle, const Particles& particles) {
  const std::size_t num_k{num_unique_k()};
  const std::size_t ROW_STRIDE{trig_row_stride()};
  const std::size_t offset{particle * ROW_STRIDE};
#if defined(__CUDACC__)
  // save
  CUDA_CHECK(cudaMemcpy(trig_scratch_[SIN_SAVED], sin_cache() + offset, num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy(trig_scratch_[COS_SAVED], cos_cache() + offset, num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
  
  // update
  dim3 updateTrigRowThreads(256);
  dim3 updateTrigRowBlocks(
    vmc::cudaNumBlocks(num_k, threads.x)
  );
  cudaUpdateTrigRow<<<updateTrigRowBlocks, updateTrigRowThreads>>>(
    particle, num_k, ROW_STRIDE,
    particles.pos().x_, particles.pos().y_, particles.pos().z_,
    k_vector().x_, k_vector().y_, k_vector().z_,
    sin_cache(), cos_cache()
  );
  CUDA_CHECK(cudaGetLastError());

  CUDA_CHECK(cudaDeviceSynchronize());
#else
  // save
  std::memcpy(trig_scratch_[SIN_SAVED], sin_cache() + offset, num_k * sizeof(real_t));
  std::memcpy(trig_scratch_[COS_SAVED], cos_cache() + offset, num_k * sizeof(real_t));

  // update
  const auto pos{particles.pos().align()};
  const auto kv{k_vector().align()};

  real_t* RESTRICT c_row{cos_cache() + particle * ROW_STRIDE};
  real_t* RESTRICT s_row{sin_cache() + particle * ROW_STRIDE};

  ASSUME_ALIGNED(c_row, SIMD_BYTES);
  ASSUME_ALIGNED(s_row, SIMD_BYTES);

  const real_t px{pos.x_[particle]};
  const real_t py{pos.y_[particle]};
  const real_t pz{pos.z_[particle]};

  // Not vectorized: loop-carried data dependency
  #pragma omp simd
  for (std::size_t k = 0; k < num_k; ++k) {
    const real_t dot{
      kv.x_[k] * px +
      kv.y_[k] * py +
      kv.z_[k] * pz
    };

    vmc::sincos(dot, &s_row[k], &c_row[k]);
  }

#endif
}
