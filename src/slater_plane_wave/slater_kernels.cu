#include "slater_plane_wave.cuh"

#if defined(__CUDACC__)
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
  real_t* RESTRICT s_cache, real_t* RESTRICT c_cache
) {
  const std::size_t j{blockIdx.x * blockDim.x + threadIdx.x};
  if (j >= num_unique_k) { return; }

  const std::size_t i{blockIdx.y * blockDim.y + threadIdx.y};
  if (i >= N) { return; }

  const real_t dot{
    k_x[j] * p_x[i] +
    k_y[j] * p_y[i] +
    k_z[j] * p_z[i]
  };

  const std::size_t offset{i * trig_S};
  const std::size_t sc_idx{offset + j};

  vmc::sincos(dot, &s_cache[sc_idx], &c_cache[sc_idx]);
}

__global__
void cudaBuildDetFromCache(
  std::size_t N, std::size_t trig_S, std::size_t mat_S,
  const real_t* RESTRICT s_cache, const real_t* RESTRICT c_cache,
  const std::size_t* RESTRICT k_idx, const std::uint8_t* RESTRICT orb_t,
  real_t* RESTRICT det_mat
) {
  const std::size_t j{blockIdx.x * blockDim.x + threadIdx.x};
  if (j >= N) { return; }

  const std::size_t i{blockIdx.y * blockDim.y + threadIdx.y};
  if (i >= N) { return; }

  const std::size_t offset{i * trig_S};
  const std::size_t trig_idx{offset + k_idx[j]};
  const std::size_t mat_idx{i * mat_S + j};

  const real_t type{static_cast<real_t>(orb_t[j])};
  const real_t cos_term{c_cache[trig_idx]};
  const real_t sin_term{s_cache[trig_idx]};

  det_mat[mat_idx] = cos_term + type * (sin_term - cos_term);
}

__global__
void cudaComputeLogAbsDet(
  std::size_t N, std::size_t mat_S,
  const real_t* lower_upper,
  real_t* log_abs_det
) {
  const std::size_t i{blockIdx.x * blockDim.x + threadIdx.x};
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
  const std::size_t i{blockIdx.x * blockDim.x + threadIdx.x};
  const std::size_t size{N * mat_S};
  if (i >= size) { return; }

  const std::size_t row{i / mat_S};
  const std::size_t col{i - row * mat_S};

  matrix[i] = (row == col) ? static_cast<real_t>(1) : static_cast<real_t>(0);
}

} // namespace

real_t SlaterPlaneWave::cudaLogAbsDet(const Particles& particles) {
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
}

#endif