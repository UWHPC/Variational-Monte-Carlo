#pragma once

#include "../particles/particles.hpp"
#include <xpu/buffer.hpp>
#include <xpu/soa.hpp>
#include "../utilities/macros.hpp"
#include <xpu/math.hpp>

#include <cstddef>
#include <cstdlib>
#include <cstring>

#ifdef XPU_CUDA
  #include <cusolverDn.h>

struct CudaScratch {
  cusolverDnHandle_t handle{};
  xpu::unique_ptr<real_t> work{};
  xpu::unique_ptr<int> info{};
  xpu::unique_ptr<real_t> log_abs_det{};
};
#endif

class SlaterPlaneWave {
private:
  std::size_t num_orbitals_;
  std::size_t num_unique_k_;
  std::size_t trig_row_stride_;
  std::size_t matrix_row_stride_;
  std::size_t matrix_size_;
  real_t box_length_;

  enum OrbKIndex : std::size_t { K, NUM_ORB_K };
  xpu::soa<std::size_t, NUM_ORB_K> orbital_k_index_;
  
  enum OrbType : std::size_t { O, NUM_ORB_TYPE };
  xpu::soa<std::uint8_t, NUM_ORB_TYPE> orbital_type_;

  enum PivotIndex : std::size_t { N_X, N_Y, N_Z, PIVOT, NUM_INT_VECTORS };
  xpu::soa<int, NUM_INT_VECTORS> int_vec_;

  enum VectorIndex : std::size_t {
    K_X,
    K_Y,
    K_Z,
    RHS,
    SOLUTION,
    NEW_ROW,
    INV_D_COL,
    NUM_DOUBLE_VECTORS
  };
  xpu::soa<real_t, NUM_DOUBLE_VECTORS> fp_vec_;

  enum TrigIndex : std::size_t { SIN_CACHE, COS_CACHE, NUM_TRIG_ARRAYS };
  xpu::soa<real_t, NUM_TRIG_ARRAYS> trig_cache_;

  enum ScratchTrigIndex : std::size_t { SIN_SAVED, COS_SAVED, NUM_SCRATCH_TRIG };
  xpu::soa<real_t, NUM_SCRATCH_TRIG> trig_scratch_;

  enum MatrixIndex : std::size_t { D, INV_D, LU, NUM_MATRIX };
  xpu::soa<real_t, NUM_MATRIX> matrices_;

#ifdef XPU_CUDA
  CudaScratch cuda_scratch_;
#endif

public:
  explicit SlaterPlaneWave(const Particles& particles, real_t box_length);
  void initialize(const Particles& particles);

  ~SlaterPlaneWave();

  SlaterPlaneWave(const SlaterPlaneWave&) = delete;
  SlaterPlaneWave& operator=(const SlaterPlaneWave&) = delete;

  SlaterPlaneWave(SlaterPlaneWave&&) = delete;
  SlaterPlaneWave& operator=(SlaterPlaneWave&&) = delete;

  [[nodiscard]] std::size_t num_orbitals() const noexcept { return num_orbitals_; }
  [[nodiscard]] std::size_t num_unique_k() const noexcept { return num_unique_k_; }
  [[nodiscard]] std::size_t trig_row_stride() const noexcept { return trig_row_stride_; }
  [[nodiscard]] std::size_t matrix_row_stride() const noexcept { return matrix_row_stride_; }
  [[nodiscard]] std::size_t matrix_size() const noexcept { return matrix_size_; }
  [[nodiscard]] real_t box_length() const noexcept { return box_length_; }

  [[nodiscard]]       std::size_t* orbital_k_index()       noexcept { return orbital_k_index_[K]; }
  [[nodiscard]] const std::size_t* orbital_k_index() const noexcept { return orbital_k_index_[K]; }

  [[nodiscard]]       std::uint8_t* orbital_type()       noexcept { return orbital_type_[O]; }
  [[nodiscard]] const std::uint8_t* orbital_type() const noexcept { return orbital_type_[O]; }

  [[nodiscard]] real_t*       determinant()       noexcept { return matrices_[D]; }
  [[nodiscard]] real_t const* determinant() const noexcept { return matrices_[D]; }

  [[nodiscard]] real_t*       inv_determinant()       noexcept { return matrices_[INV_D]; }
  [[nodiscard]] real_t const* inv_determinant() const noexcept { return matrices_[INV_D]; }

  [[nodiscard]] real_t*       lower_upper()       noexcept { return matrices_[LU]; }
  [[nodiscard]] real_t const* lower_upper() const noexcept { return matrices_[LU]; }

  [[nodiscard]] int*       pivot()       noexcept { return int_vec_[PIVOT]; }
  [[nodiscard]] int const* pivot() const noexcept { return int_vec_[PIVOT]; }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<int, idx(Axis::NUM)> n_vector() {
    return int_vec_.view<idx(Axis::NUM), N_X>();
  }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<const int, idx(Axis::NUM)> n_vector() const {
    return int_vec_.view<idx(Axis::NUM), N_X>();
  }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<real_t, idx(Axis::NUM)> k_vector() {
    return fp_vec_.view<idx(Axis::NUM), K_X>();
  }

  [[nodiscard]] CUDA_CALLABLE
  xpu::soa_view<const real_t, idx(Axis::NUM)> k_vector() const {
    return fp_vec_.view<idx(Axis::NUM), K_X>();
  }

  [[nodiscard]] real_t*       solution()       noexcept { return fp_vec_[SOLUTION]; }
  [[nodiscard]] real_t const* solution() const noexcept { return fp_vec_[SOLUTION]; }

  [[nodiscard]] real_t*       rhs()       noexcept { return fp_vec_[RHS]; }
  [[nodiscard]] real_t const* rhs() const noexcept { return fp_vec_[RHS]; }

  [[nodiscard]] real_t*       sin_cache()       noexcept { return trig_cache_[SIN_CACHE]; }
  [[nodiscard]] real_t const* sin_cache() const noexcept { return trig_cache_[SIN_CACHE]; }

  [[nodiscard]] real_t*       cos_cache()       noexcept { return trig_cache_[COS_CACHE]; }
  [[nodiscard]] real_t const* cos_cache() const noexcept { return trig_cache_[COS_CACHE]; }

  void restore_trig_row(std::size_t particle);
  void update_trig_cache(std::size_t particle, Particles& particles);

  real_t log_abs_det(const Particles& particles);

  real_t* build_row(std::size_t particle) noexcept;

  [[nodiscard]] real_t determinant_ratio(
    std::size_t particle,
    const real_t* new_row
  ) const noexcept;

  void accept_move(std::size_t particle, const real_t* new_row, real_t ratio) noexcept;

  void add_derivatives(
    xpu::soa_view<real_t, idx(Derivatives::NUM)> derivatives
  ) const noexcept;

private:
  void save_trig_row(std::size_t particle);
#ifdef XPU_CUDA
  void init_cuda_scratch();
#endif

  [[nodiscard]] real_t* new_row() noexcept { return fp_vec_[NEW_ROW]; }
  [[nodiscard]] real_t* inv_d_col() noexcept { return fp_vec_[INV_D_COL]; }
};
