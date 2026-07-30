#include "slater_plane_wave.cuh"
#include "slater_kernels.cuh"

#ifdef VMC_CUDA_BACKEND
#include <cstddef>
namespace {

__global__ void cudaUpdateTrigRow(
  std::size_t num_k, std::size_t offset,
  real_t p_x, real_t p_y, real_t p_z,
  const real_t* RESTRICT k_x, const real_t* RESTRICT k_y, const real_t* RESTRICT k_z,
  real_t* RESTRICT sin_cache, real_t* RESTRICT cos_cache
) {
  const auto [i]{vmc::cudaThreadIdx<1>()};
  if (i >= num_k) { return; }

  trig_cell(i, offset, p_x, p_y, p_z, k_x, k_y, k_z, sin_cache, cos_cache);
}

} // namespace
#else
#include "particles/particles.cuh"
#include "utilities/aligned_soa.cuh"
#include <cstring>
#endif


void SlaterPlaneWave::restore_trig_row(std::size_t particle) {
  const std::size_t num_k{this->num_unique_k()};
  const std::size_t ROW_STRIDE{this->trig_row_stride()};
  const std::size_t offset{particle * ROW_STRIDE};

#ifdef VMC_CUDA_BACKEND
  CUDA_CHECK(cudaMemcpy(this->sin_cache() + offset, trig_scratch_[SIN_SAVED], num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy(this->cos_cache() + offset, trig_scratch_[COS_SAVED], num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
#else
  std::memcpy(this->sin_cache() + offset, trig_scratch_[SIN_SAVED], num_k * sizeof(real_t));
  std::memcpy(this->cos_cache() + offset, trig_scratch_[COS_SAVED], num_k * sizeof(real_t));
#endif
}

void SlaterPlaneWave::update_trig_cache(
  std::size_t particle, const Particles& particles) {
  const std::size_t num_k{this->num_unique_k()};
  const std::size_t ROW_STRIDE{this->trig_row_stride()};
  const std::size_t offset{particle * ROW_STRIDE};
#ifdef VMC_CUDA_BACKEND
  // save
  CUDA_CHECK(cudaMemcpy(trig_scratch_[SIN_SAVED], this->sin_cache() + offset, num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy(trig_scratch_[COS_SAVED], this->cos_cache() + offset, num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
  
  // update
  dim3 updateTrigRowThreads(256);
  dim3 updateTrigRowBlocks(
    vmc::cudaNumBlocks(num_k, updateTrigRowThreads.x)
  );
  cudaUpdateTrigRow<<<updateTrigRowBlocks, updateTrigRowThreads>>>(
    num_k, particle * ROW_STRIDE,
    particles.pos().x_[particle], particles.pos().y_[particle], particles.pos().z_[particle],
    this->k_vector().x_, this->k_vector().y_, this->k_vector().z_,
    this->sin_cache(), this->cos_cache()
  );
  CUDA_CHECK(cudaGetLastError());

  CUDA_CHECK(cudaDeviceSynchronize());
#else
  // save
  std::memcpy(trig_scratch_[SIN_SAVED], this->sin_cache() + offset, num_k * sizeof(real_t));
  std::memcpy(trig_scratch_[COS_SAVED], this->cos_cache() + offset, num_k * sizeof(real_t));

  // update
  const auto pos{particles.pos().align()};
  const auto kv{this->k_vector().align()};

  real_t* RESTRICT cos_cache{this->cos_cache()}; ASSUME_ALIGNED(cos_cache, SIMD_BYTES);
  real_t* RESTRICT sin_cache{this->sin_cache()}; ASSUME_ALIGNED(sin_cache, SIMD_BYTES);

  const real_t px{pos.x_[particle]};
  const real_t py{pos.y_[particle]};
  const real_t pz{pos.z_[particle]};

  // Not vectorized: loop-carried data dependency
  #pragma omp simd
  for (std::size_t k = 0; k < num_k; ++k) {
    trig_cell(k, offset, px, py, pz, kv.x_, kv.y_, kv.z_, sin_cache, cos_cache);
  }
#endif
}
