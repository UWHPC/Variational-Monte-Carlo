#pragma once

#include "../particles/particles.cuh"
#include "../utilities/aligned_soa.cuh"
#include "../utilities/macros.cuh"
#include "../utilities/math.cuh"
#include "../utilities/ptr3d.hpp"

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <numbers>
#include <vector>

class SlaterPlaneWave {
private:
  std::size_t num_orbitals_;
  std::size_t num_unique_k_;
  std::size_t trig_row_stride_;
  std::size_t matrix_row_stride_;
  std::size_t matrix_size_;
  real_t box_length_;

  enum OrbKIndex : std::size_t { K, NUM_ORB_K };
  AlignedSoA<std::size_t> orbital_k_index_;
  
  enum OrbType : std::size_t { O, NUM_ORB_TYPE };
  AlignedSoA<std::uint8_t> orbital_type_;

  enum PivotIndex : std::size_t { N_X, N_Y, N_Z, PIVOT, NUM_INT_VECTORS };
  AlignedSoA<int> int_vec_;

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
  AlignedSoA<real_t> fp_vec_;

  enum TrigIndex : std::size_t { SIN_CACHE, COS_CACHE, NUM_TRIG_ARRAYS };
  AlignedSoA<real_t> trig_cache_;

  enum ScratchTrigIndex : std::size_t { SIN_SAVED, COS_SAVED, NUM_SCRATCH_TRIG };
  AlignedSoA<real_t> trig_scratch_;

  enum MatrixIndex : std::size_t { D, INV_D, LU, NUM_MATRIX };
  AlignedSoA<real_t> matrices_;

public:
  explicit SlaterPlaneWave(const Particles& particles, real_t box_length);

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

  Ptr3D<      int> n_vector()       noexcept { return {int_vec_[N_X], int_vec_[N_Y], int_vec_[N_Z]}; }
  Ptr3D<const int> n_vector() const noexcept { return {int_vec_[N_X], int_vec_[N_Y], int_vec_[N_Z]}; }

  Ptr3D<      real_t> k_vector()       noexcept { return {fp_vec_[K_X], fp_vec_[K_Y], fp_vec_[K_Z]}; }
  Ptr3D<const real_t> k_vector() const noexcept { return {fp_vec_[K_X], fp_vec_[K_Y], fp_vec_[K_Z]}; }

  [[nodiscard]] real_t*       solution()       noexcept { return fp_vec_[SOLUTION]; }
  [[nodiscard]] real_t const* solution() const noexcept { return fp_vec_[SOLUTION]; }

  [[nodiscard]] real_t*       rhs()       noexcept { return fp_vec_[RHS]; }
  [[nodiscard]] real_t const* rhs() const noexcept { return fp_vec_[RHS]; }

  [[nodiscard]] real_t*       sin_cache()       noexcept { return trig_cache_[SIN_CACHE]; }
  [[nodiscard]] real_t const* sin_cache() const noexcept { return trig_cache_[SIN_CACHE]; }

  [[nodiscard]] real_t*       cos_cache()       noexcept { return trig_cache_[COS_CACHE]; }
  [[nodiscard]] real_t const* cos_cache() const noexcept { return trig_cache_[COS_CACHE]; }

  void save_trig_row(std::size_t particle) noexcept;
  void restore_trig_row(std::size_t particle) noexcept;

  void update_trig_cache(std::size_t particle, const Particles& particles) noexcept;

  #if defined(__CUDACC__)
  real_t cudaLogAbsDet(const Particles& particles);
  #endif
  real_t cpu_log_abs_det(const Particles& particles);

  real_t log_abs_det(const Particles& particles) {
  #if defined(__CUDACC__)
    return cudaLogAbsDet(particles);
  #else
    return cpu_log_abs_det(particles);
  #endif
  }

  real_t* build_row(std::size_t particle) noexcept;

  [[nodiscard]] real_t determinant_ratio(
    std::size_t particle,
    const real_t* new_row
  ) const noexcept;

  void accept_move(std::size_t particle, const real_t* new_row, real_t ratio) noexcept;

  void add_derivatives(
    real_t* RESTRICT grad_x,
    real_t* RESTRICT grad_y,
    real_t* RESTRICT grad_z,
    real_t* RESTRICT laplacian
  ) const noexcept;

private:
  [[nodiscard]] real_t* new_row() noexcept { return fp_vec_[NEW_ROW]; }
  [[nodiscard]] real_t* inv_d_col() noexcept { return fp_vec_[INV_D_COL]; }
};