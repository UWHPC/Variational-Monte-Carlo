#include <xpu/xpu.hpp>
#include "slater_plane_wave.hpp"
#include "slater_plane_wave_kernels.hpp"

#if defined(XPU_CUDA)
  #include <cublas_v2.h>
  #include <cstdio>
  #include <cstdlib>
#else
  #include "../utilities/matrix.hpp"
#endif

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

#if defined(XPU_CUDA)

namespace {

inline void cuSolverGetrfBufferSize(
  cusolverDnHandle_t handle,
  int rows,
  int columns,
  real_t* matrix,
  int leading_dim,
  int* work_size
) {
#if defined(FP_64)
  xpu::cu_check(cusolverDnDgetrf_bufferSize(
    handle, rows, columns, matrix, leading_dim, work_size
  ));
#else
  xpu::cu_check(cusolverDnSgetrf_bufferSize(
    handle, rows, columns, matrix, leading_dim, work_size
  ));
#endif
}

inline void cusolverGetrf(
  cusolverDnHandle_t handle,
  int rows,
  int columns,
  real_t* matrix,
  int leading_dim,
  real_t* work,
  int* pivot,
  int* info
) {
#if defined(FP_64)
  xpu::cu_check(cusolverDnDgetrf(
    handle, rows, columns, matrix, leading_dim, work, pivot, info
  ));
#else
  xpu::cu_check(cusolverDnSgetrf(
    handle, rows, columns, matrix, leading_dim, work, pivot, info
  ));
#endif
}

inline void cusolverGetrs(
  cusolverDnHandle_t handle,
  cublasOperation_t transpose,
  int num_orbitals,
  int num_rhs,
  const real_t* lower_upper,
  int leading_dim,
  const int* pivot,
  real_t* rhs,
  int rhs_leading_dim,
  int* info
) {
#if defined(FP_64)
  xpu::cu_check(cusolverDnDgetrs(
    handle, transpose, num_orbitals, num_rhs,
    lower_upper, leading_dim,
    pivot, rhs, rhs_leading_dim,
    info
  ));
#else
  xpu::cu_check(cusolverDnSgetrs(
    handle, transpose, num_orbitals, num_rhs,
    lower_upper, leading_dim,
    pivot, rhs, rhs_leading_dim,
    info
  ));
#endif
}

void checkCusolverInfo(int info, const char* name) {
  if (info < 0) {
    std::fprintf(
      stderr,
      "cuSOLVER %s invalid argument: %d\n",
      name, -info
    );
    std::abort();
  }
}

} // namespace

void SlaterPlaneWave::init_cuda_scratch() {
  const auto num_orbitals{static_cast<int>(this->num_orbitals())};
  const auto matrix_row_stride{
    static_cast<int>(this->matrix_row_stride())
  };

  xpu::cu_check(cusolverDnCreate(&cuda_scratch_.handle));

  auto work_size{0};
  cuSolverGetrfBufferSize(
    cuda_scratch_.handle,
    num_orbitals, num_orbitals,
    this->lower_upper(), matrix_row_stride,
    &work_size
  );

  cuda_scratch_.work = xpu::unique_ptr<real_t>{
    static_cast<std::size_t>(work_size)
  };
  cuda_scratch_.info = xpu::unique_ptr<int>{1uz};
  cuda_scratch_.log_abs_det = xpu::unique_ptr<real_t>{1uz};
}

#endif

SlaterPlaneWave::~SlaterPlaneWave() {
#if defined(XPU_CUDA)
  if (cuda_scratch_.handle) {
    xpu::cu_check(cusolverDnDestroy(cuda_scratch_.handle));
  }
#endif
}

real_t SlaterPlaneWave::log_abs_det(const Particles& particles) {
  kernel::slater::build_trig_cache(
    this->num_unique_k(),
    this->trig_row_stride(),
    particles.pos(), k_vector(),
    this->sin_cache(), this->cos_cache()
  );

  kernel::slater::build_determinant(
    this->num_orbitals(),
    this->trig_row_stride(), this->matrix_row_stride(),
    this->sin_cache(), this->cos_cache(),
    this->orbital_k_index(), this->orbital_type(),
    this->determinant()
  );

  xpu::copy_n(
    this->lower_upper(),
    this->determinant(),
    this->matrix_size()
  );

#if defined(XPU_CUDA)
  const auto num_orbitals{static_cast<int>(this->num_orbitals())};
  const auto matrix_row_stride{
    static_cast<int>(this->matrix_row_stride())
  };
  auto& scratch{cuda_scratch_};

  cusolverGetrf(
    scratch.handle,
    num_orbitals, num_orbitals,
    this->lower_upper(), matrix_row_stride,
    scratch.work.get(),
    this->pivot(),
    scratch.info.get()
  );

  auto info{0};
  xpu::copy_n(&info, scratch.info.get(), 1uz);
  checkCusolverInfo(info, "getrf");
  if (info > 0) {
    return -std::numeric_limits<real_t>::infinity();
  }

  const auto log_abs_det{
    kernel::slater::compute_log_abs_det(
      this->num_orbitals(),
      this->matrix_row_stride(),
      this->lower_upper(),
      scratch.log_abs_det.get()
    )
  };
  if (!std::isfinite(log_abs_det)) {
    return -std::numeric_limits<real_t>::infinity();
  }

  kernel::slater::build_identity(
    this->matrix_size(),
    this->matrix_row_stride(),
    this->inv_determinant()
  );

  cusolverGetrs(
    scratch.handle,
    CUBLAS_OP_T,
    num_orbitals, num_orbitals,
    this->lower_upper(), matrix_row_stride,
    this->pivot(),
    this->inv_determinant(), matrix_row_stride,
    scratch.info.get()
  );

  xpu::copy_n(&info, scratch.info.get(), 1uz);
  checkCusolverInfo(info, "getrs");
#else
  static_cast<void>(lower_upper_decomp(
    this->lower_upper(),
    this->pivot(),
    this->num_orbitals(),
    this->matrix_row_stride()
  ));

  const auto log_abs_det{
    kernel::slater::compute_log_abs_det(
      this->num_orbitals(),
      this->matrix_row_stride(),
      this->lower_upper(),
      nullptr
    )
  };
  if (!std::isfinite(log_abs_det)) {
    return -std::numeric_limits<real_t>::infinity();
  }

  for (auto column = 0uz; column < this->num_orbitals(); ++column) {
    xpu::zero_n(this->rhs(), this->num_orbitals());
    this->rhs()[column] = 1.0_r;

    solve_lower_upper(
      this->lower_upper(),
      this->pivot(),
      this->rhs(), this->solution(),
      this->num_orbitals(),
      this->matrix_row_stride()
    );

    xpu::copy_n(
      &this->inv_determinant()[column * this->matrix_row_stride()],
      this->solution(),
      this->num_orbitals()
    );
  }
#endif

  return log_abs_det;
}
