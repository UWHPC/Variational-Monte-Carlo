#include "slater_plane_wave.cuh"
#include "slater_kernels.cuh"

#ifdef VMC_CUDA_BACKEND
#include <cublas_v2.h>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace {

inline void cuSolverGetrfBufferSize(
  cusolverDnHandle_t handle,
  int rows,
  int cols,
  real_t* matrix,
  int leading_dim,
  int* work_size
) {
#ifdef FP_64
  CUSOLVER_CHECK(cusolverDnDgetrf_bufferSize(
    handle, rows, cols, matrix, leading_dim, work_size
  ));
#else
  CUSOLVER_CHECK(cusolverDnSgetrf_bufferSize(
    handle, rows, cols, matrix, leading_dim, work_size
  ));
#endif
}

inline void cusolverGetrf(
  cusolverDnHandle_t handle,
  int rows,
  int cols,
  real_t* matrix,
  int leading_dim,
  real_t* work,
  int* pivot,
  int* info
) {
#ifdef FP_64
  CUSOLVER_CHECK(cusolverDnDgetrf(
    handle, rows, cols, matrix, leading_dim, work, pivot, info
  ));
#else
  CUSOLVER_CHECK(cusolverDnSgetrf(
    handle, rows, cols, matrix, leading_dim, work, pivot, info
  ));
#endif
}

inline void cusolverGetrs(
  cusolverDnHandle_t handle,
  cublasOperation_t transpose,
  int n,
  int num_rhs,
  const real_t* lower_upper,
  int leading_dim,
  const int* pivot,
  real_t* rhs,
  int rhs_leading_dim,
  int* info
) {
#ifdef FP_64
  CUSOLVER_CHECK(cusolverDnDgetrs(
    handle, transpose, n, num_rhs, lower_upper, leading_dim,
    pivot, rhs, rhs_leading_dim, info
  ));
#else
  CUSOLVER_CHECK(cusolverDnSgetrs(
    handle, transpose, n, num_rhs, lower_upper, leading_dim,
    pivot, rhs, rhs_leading_dim, info
  ));
#endif
}

void checkCusolverInfo(int info, const char* name) {
  if (info < 0) {
    std::fprintf(stderr, "cuSOLVER %s invalid argument: %d\n", name, -info);
    std::abort();
  }
}

__global__
void cudaBuildTrigCache(
  std::size_t N, std::size_t trig_S, std::size_t num_unique_k,
  const real_t* RESTRICT p_x, const real_t* RESTRICT p_y, const real_t* RESTRICT p_z,
  const real_t* RESTRICT k_x, const real_t* RESTRICT k_y, const real_t* RESTRICT k_z,
  real_t* RESTRICT sin_cache, real_t* RESTRICT cos_cache
) {
  const auto [i, j]{vmc::cudaThreadIdx<2>()};
  if (i >= num_unique_k || j >= N) {return; }

  const std::size_t offset{j * trig_S};

  trig_cell(i, offset, p_x[j], p_y[j], p_z[j], k_x, k_y, k_z, sin_cache, cos_cache);
}

__global__
void cudaBuildDetFromCache(
  std::size_t N, std::size_t trig_S, std::size_t mat_S,
  const real_t* RESTRICT sin_cache, const real_t* RESTRICT cos_cache,
  const std::size_t* RESTRICT k_idx, const std::uint8_t* RESTRICT orb_t,
  real_t* RESTRICT det_mat
) {
  const auto [i, j]{vmc::cudaThreadIdx<2>()};
  if (i >= N || j >= N) { return; }

  build_row_cell(i, j, trig_S, sin_cache, cos_cache, k_idx[i], orb_t[i], det_mat + j * mat_S);
}

__global__
void cudaComputeLogAbsDet(
  std::size_t N, std::size_t mat_S,
  const real_t* lower_upper,
  real_t* log_abs_det
) {
  const auto [i]{vmc::cudaThreadIdx<1>()};
  if (i >= N) { return; }

  const real_t U_diag{lower_upper[i * mat_S + i]};

  atomicAdd(log_abs_det, vmc::log(vmc::abs(U_diag)));
}

__global__
void cudaBuildIdentity(
  std::size_t N,
  std::size_t mat_S,
  real_t* matrix
) {
  const auto [i]{vmc::cudaThreadIdx<1>()};
  if (i >= N * mat_S) { return; }

  const std::size_t row{i / mat_S};
  const std::size_t col{i - row * mat_S};

  matrix[i] = (row == col) ? static_cast<real_t>(1) : static_cast<real_t>(0);
}

} // namespace
#else
#include "slater_plane_wave.cuh"
#include "../utilities/matrix.hpp"
#include "particles/particles.cuh"
#include "utilities/aligned_soa.cuh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <vector>
#endif


real_t SlaterPlaneWave::log_abs_det(const Particles& particles) {
#ifdef VMC_CUDA_BACKEND
  AlignedSoA<real_t> log_abs_det{1, 1};

  const int N{static_cast<int>(particles.size())};
  const int mat_S{static_cast<int>(this->matrix_row_stride())};

  cusolverDnHandle_t cusolver_handle{};
  CUSOLVER_CHECK(cusolverDnCreate(&cusolver_handle));

  int work_size{};
  cuSolverGetrfBufferSize(
    cusolver_handle,
    N, N,
    this->lower_upper(), mat_S,
    &work_size
  );

  AlignedSoA<real_t> work{static_cast<std::size_t>(work_size), 1};
  AlignedSoA<int> info{1, 1};

  dim3 buildTrigCacheThreads(16, 16);
  dim3 buildTrigCacheBlocks(
    vmc::cudaNumBlocks(num_unique_k(), buildTrigCacheThreads.x),
    vmc::cudaNumBlocks(particles.size(), buildTrigCacheThreads.y)
  );
  cudaBuildTrigCache<<<buildTrigCacheBlocks, buildTrigCacheThreads>>>(
    particles.size(), this->trig_row_stride(), this->num_unique_k(),
    particles.pos().x_, particles.pos().y_, particles.pos().z_,
    this->k_vector().x_, this->k_vector().y_, this->k_vector().z_,
    this->sin_cache(), this->cos_cache()
  );
  CUDA_CHECK(cudaGetLastError());

  dim3 buildDetFromCacheThreads(16, 16);
  dim3 buildDetFromCacheBlocks(
    vmc::cudaNumBlocks(particles.size(), buildDetFromCacheThreads.x),
    vmc::cudaNumBlocks(particles.size(), buildDetFromCacheThreads.y)
  );
  cudaBuildDetFromCache<<<buildDetFromCacheBlocks, buildDetFromCacheThreads>>>(
    particles.size(),
    this->trig_row_stride(), this->matrix_row_stride(),
    this->sin_cache(), this->cos_cache(),
    this->orbital_k_index(), this->orbital_type(),
    this->determinant()
  );
  CUDA_CHECK(cudaGetLastError());

  CUDA_CHECK(cudaMemcpyAsync(
    this->lower_upper(),
    this->determinant(),
    this->matrix_row_stride() *
    particles.size()          *
    sizeof(real_t),
    cudaMemcpyDeviceToDevice
  ));

  cusolverGetrf(
    cusolver_handle,
    N, N,
    this->lower_upper(), mat_S,
    work[0], this->pivot(), info[0]
  );

  CUDA_CHECK(cudaDeviceSynchronize());
  checkCusolverInfo(*info[0], "getrf");
  if (*info[0] > 0) {
    CUSOLVER_CHECK(cusolverDnDestroy(cusolver_handle));
    return -std::numeric_limits<real_t>::infinity();
  }

  dim3 computeLogAbsDetThreads(256);
  dim3 computeLogAbsDetBlocks(
    vmc::cudaNumBlocks(particles.size(), computeLogAbsDetThreads.x)
  );
  cudaComputeLogAbsDet<<<computeLogAbsDetBlocks, computeLogAbsDetThreads>>>(
    particles.size(), this->matrix_row_stride(),
    this->lower_upper(),
    log_abs_det[0]
  );
  CUDA_CHECK(cudaGetLastError());

  dim3 buildIdentityThreads(256);
  dim3 buildIdentityBlocks(
    vmc::cudaNumBlocks(
      particles.size() * this->matrix_row_stride(),
      buildIdentityThreads.x
    )
  );
  cudaBuildIdentity<<<buildIdentityBlocks, buildIdentityThreads>>>(
    particles.size(),
    this->matrix_row_stride(),
    this->inv_determinant()
  );
  CUDA_CHECK(cudaGetLastError());

  cusolverGetrs(
    cusolver_handle,
    CUBLAS_OP_T,
    N, N,
    this->lower_upper(), mat_S,
    this->pivot(),
    this->inv_determinant(), mat_S,
    info[0]
  );

  CUDA_CHECK(cudaDeviceSynchronize());
  checkCusolverInfo(*info[0], "getrs");
  CUSOLVER_CHECK(cusolverDnDestroy(cusolver_handle));

  if (!std::isfinite(*log_abs_det[0])) {
    return -std::numeric_limits<real_t>::infinity();
  }

  return *log_abs_det[0];
#else
  const std::size_t N{this->num_orbitals()};
  const std::size_t S{this->matrix_row_stride()};
  const std::size_t padded_N{particles.p_stride()};

  const auto pos{particles.pos().align()};
  const auto kv{this->k_vector().align()};

  const auto* k_index{this->orbital_k_index()};
  const auto* orb_type{this->orbital_type()};

  real_t* RESTRICT det_matrix{this->determinant()}; ASSUME_ALIGNED(det_matrix, SIMD_BYTES);
  real_t* RESTRICT lower_upper_matrix{this->lower_upper()}; ASSUME_ALIGNED(lower_upper_matrix, SIMD_BYTES);
  real_t* RESTRICT inv_det_matrix{this->inv_determinant()}; ASSUME_ALIGNED(inv_det_matrix, SIMD_BYTES);

  int* RESTRICT pivot_vector{this->pivot()}; ASSUME_ALIGNED(pivot_vector, SIMD_BYTES);

  real_t* RESTRICT rhs{this->rhs()}; ASSUME_ALIGNED(rhs, SIMD_BYTES);
  real_t* RESTRICT solution{this->solution()}; ASSUME_ALIGNED(solution, SIMD_BYTES);

  const std::size_t num_k{this->num_unique_k()};
  real_t* RESTRICT cos_cache{this->cos_cache()}; ASSUME_ALIGNED(cos_cache, SIMD_BYTES);
  real_t* RESTRICT sin_cache{this->sin_cache()}; ASSUME_ALIGNED(sin_cache, SIMD_BYTES);

  const std::size_t ROW_STRIDE{this->trig_row_stride()};
  for (std::size_t particle = 0; particle < N; ++particle) {
    const std::size_t offset{particle * ROW_STRIDE};
    const real_t px{pos.x_[particle]};
    const real_t py{pos.y_[particle]};
    const real_t pz{pos.z_[particle]};

    // Not vectorized: data dependency
    #pragma omp simd
    for (std::size_t k = 0; k < num_k; ++k) {
      trig_cell(k, offset, px, py, pz, kv.x_, kv.y_, kv.z_, sin_cache, cos_cache);
    }

    real_t* RESTRICT det_row{det_matrix + particle * S};

    // Not vectorized: non-contiguous memory accesses
    #pragma omp simd
    for (std::size_t orbital = 0; orbital < N; ++orbital) {
      build_row_cell(orbital, particle, ROW_STRIDE, sin_cache, cos_cache, k_index[orbital], orb_type[orbital], det_row);
    }
  }

  std::copy_n(det_matrix, S * N, lower_upper_matrix);

  static_cast<void>(
    lower_upper_decomp(lower_upper_matrix, pivot_vector, N, S)
  );

  real_t log_abs_det{};

  // Not vectorzied: loop-carried data dependency
  #pragma omp simd reduction(+ : log_abs_det)
  for (std::size_t diag = 0; diag < N; ++diag) {
    const real_t U_ii{lower_upper_matrix[diag * S + diag]};
    const real_t abs_U_ii{vmc::abs(U_ii)};

    log_abs_det += vmc::log(abs_U_ii);
  }

  if (!std::isfinite(log_abs_det)) {
    return -std::numeric_limits<real_t>::infinity();
  }

  for (std::size_t column = 0; column < N; ++column) {
    std::fill(rhs, rhs + padded_N, 0.0_r);
    rhs[column] = 1.0_r;

    solve_lower_upper(lower_upper_matrix, pivot_vector, rhs, solution, N, S);

    for (std::size_t row = 0; row < N; ++row) {
      inv_det_matrix[column * S + row] = solution[row];
    }
  }

  return log_abs_det;
#endif
}
