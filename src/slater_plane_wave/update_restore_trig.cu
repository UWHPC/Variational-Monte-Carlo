#include <xpu/xpu.hpp>
#include "slater_plane_wave.hpp"

#ifdef XPU_CUDA
#include <cstddef>
namespace {

__global__ void cudaUpdateTrigRow(
  std::size_t num_k, std::size_t offset,
  real_t p_x, real_t p_y, real_t p_z,
  xpu::soa_view<const real_t, 3> k_vector,
  real_t* RESTRICT s_cache, real_t* RESTRICT c_cache
) {
  const auto [i]{xpu::global_index<1>()};
  if (i >= num_k) { return; }

  const real_t dot{
    k_vector[idx(Axis::X)][i] * p_x +
    k_vector[idx(Axis::Y)][i] * p_y +
    k_vector[idx(Axis::Z)][i] * p_z
  };

  const auto idx{offset + i};
  xpu::sincos(dot, &s_cache[idx], &c_cache[idx]);
}

} // namespace
#else
#include "particles/particles.hpp"
#include "utilities/aligned_soa.hpp"
#include <cstring>
#endif


void SlaterPlaneWave::restore_trig_row(std::size_t particle) {
  const std::size_t num_k{this->num_unique_k()};
  const std::size_t ROW_STRIDE{this->trig_row_stride()};
  const std::size_t offset{particle * ROW_STRIDE};

#ifdef XPU_CUDA
  xpu::cu_check(cudaMemcpy(this->sin_cache() + offset, trig_scratch_[SIN_SAVED], num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
  xpu::cu_check(cudaMemcpy(this->cos_cache() + offset, trig_scratch_[COS_SAVED], num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
#else
  std::memcpy(this->sin_cache() + offset, trig_scratch_[SIN_SAVED], num_k * sizeof(real_t));
  std::memcpy(this->cos_cache() + offset, trig_scratch_[COS_SAVED], num_k * sizeof(real_t));
#endif
}

void SlaterPlaneWave::update_trig_cache(
  std::size_t particle, const Particles& particles
) {
  const std::size_t num_k{this->num_unique_k()};
  const std::size_t ROW_STRIDE{this->trig_row_stride()};
  const std::size_t offset{particle * ROW_STRIDE};
#ifdef XPU_CUDA
  // save
  xpu::cu_check(cudaMemcpy(trig_scratch_[SIN_SAVED], this->sin_cache() + offset, num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
  xpu::cu_check(cudaMemcpy(trig_scratch_[COS_SAVED], this->cos_cache() + offset, num_k * sizeof(real_t), cudaMemcpyDeviceToDevice));
  
  // update
  dim3 updateTrigRowThreads{256};
  dim3 updateTrigRowBlocks{
    xpu::block_per_dim(k_vector().count(), updateTrigRowThreads.x)
  };
  cudaUpdateTrigRow<<<updateTrigRowBlocks, updateTrigRowThreads>>>(
    num_k, particle * ROW_STRIDE,
    particles.pos().x_[particle], particles.pos().y_[particle], particles.pos().z_[particle],
    this->k_vector(),
    this->sin_cache(), this->cos_cache()
  );
  xpu::cu_check(cudaGetLastError());

  xpu::cu_check(cudaDeviceSynchronize());
#else
  // save
  std::memcpy(trig_scratch_[SIN_SAVED], this->sin_cache() + offset, num_k * sizeof(real_t));
  std::memcpy(trig_scratch_[COS_SAVED], this->cos_cache() + offset, num_k * sizeof(real_t));

  // update
  const auto pos{particles.pos()};
  const auto kv{this->k_vector()};

  real_t* RESTRICT c_row{this->cos_cache() + particle * ROW_STRIDE};
  real_t* RESTRICT s_row{this->sin_cache() + particle * ROW_STRIDE};

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

    xpu::sincos(dot, &s_row[k], &c_row[k]);
  }
#endif
}
