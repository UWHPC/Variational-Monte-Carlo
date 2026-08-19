#pragma once

#include "../particles/particles.hpp"
#include <xpu/buffer.hpp>
#include <xpu/linear_algebra.hpp>
#include <xpu/soa.hpp>
#include "../utilities/macros.hpp"
#include <xpu/math.hpp>

#include <cstddef>
#include <cstdlib>
#include <cstring>

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

  enum IntVectorIndex : std::size_t { N_X, N_Y, N_Z, NUM_INT_VECTORS };
  xpu::soa<int, NUM_INT_VECTORS> int_vec_;

  enum VectorIndex : std::size_t {
    K_X,
    K_Y,
    K_Z,
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

  xpu::buffer<real_t> reduction_scratch_;
  xpu::linalg::lu_factorization<real_t> lu_factorization_;

public:
  struct View {
    std::size_t num_orbitals;
    std::size_t num_unique_k;
    std::size_t trig_row_stride;
    std::size_t matrix_row_stride;

    xpu::soa_view<real_t, idx(Axis::NUM)> k_vector;
    const std::size_t* orbital_k_index;
    const std::uint8_t* orbital_type;

    real_t* determinant;
    real_t* inv_determinant;
    real_t* solution;
    real_t* new_row;
    real_t* inv_d_col;

    real_t* sin_cache;
    real_t* cos_cache;
    real_t* sin_saved;
    real_t* cos_saved;
  };

  explicit SlaterPlaneWave(const Particles& particles, real_t box_length);
  void initialize(const Particles& particles);

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

  [[nodiscard]] CUDA_CALLABLE
  View view() noexcept {
    return {
      this->num_orbitals(),
      this->num_unique_k(),
      this->trig_row_stride(),
      this->matrix_row_stride(),
      this->k_vector(),
      this->orbital_k_index(),
      this->orbital_type(),
      this->determinant(),
      this->inv_determinant(),
      this->solution(),
      this->new_row(),
      this->inv_d_col(),
      this->sin_cache(),
      this->cos_cache(),
      trig_scratch_[SIN_SAVED],
      trig_scratch_[COS_SAVED]
    };
  }

private:
  void save_trig_row(std::size_t particle);

  [[nodiscard]] real_t* new_row() noexcept { return fp_vec_[NEW_ROW]; }
  [[nodiscard]] real_t* inv_d_col() noexcept { return fp_vec_[INV_D_COL]; }
  [[nodiscard]] real_t* reduction_scratch() noexcept { return reduction_scratch_.data(); }
};
