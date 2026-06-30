#pragma once

#include "../particles/particles.hpp"
#include "../utilities/aligned_soa.hpp"
#include "../utilities/macros.hpp"
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
  double box_length_;

  std::vector<std::size_t> orbital_k_index_;
  std::vector<std::uint8_t> orbital_type_;

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
  AlignedSoA<double> double_vec_;

  enum TrigIndex : std::size_t { SIN_CACHE, COS_CACHE, NUM_TRIG_ARRAYS };
  AlignedSoA<double> trig_cache_;

  enum ScratchTrigIndex : std::size_t { SIN_SAVED, COS_SAVED, NUM_SCRATCH_TRIG };
  AlignedSoA<double> trig_scratch_;

  enum MatrixIndex : std::size_t { D, INV_D, LU, NUM_MATRIX };
  AlignedSoA<double> matrices_;

public:
  explicit SlaterPlaneWave(const Particles& particles, double box_length);

  [[nodiscard]] std::size_t num_orbitals() const noexcept { return num_orbitals_; }
  [[nodiscard]] std::size_t num_unique_k() const noexcept { return num_unique_k_; }
  [[nodiscard]] std::size_t trig_row_stride() const noexcept { return trig_row_stride_; }
  [[nodiscard]] std::size_t matrix_row_stride() const noexcept { return matrix_row_stride_; }
  [[nodiscard]] std::size_t matrix_size() const noexcept { return matrix_size_; }
  [[nodiscard]] double box_length() const noexcept { return box_length_; }

  [[nodiscard]]       std::vector<std::size_t>& orbital_k_index()       noexcept { return orbital_k_index_; }
  [[nodiscard]] const std::vector<std::size_t>& orbital_k_index() const noexcept { return orbital_k_index_; }

  [[nodiscard]]       std::vector<std::uint8_t>& orbital_type()       noexcept { return orbital_type_; }
  [[nodiscard]] const std::vector<std::uint8_t>& orbital_type() const noexcept { return orbital_type_; }

  [[nodiscard]] double*       determinant()       noexcept { return matrices_[D]; }
  [[nodiscard]] double const* determinant() const noexcept { return matrices_[D]; }

  [[nodiscard]] double*       inv_determinant()       noexcept { return matrices_[INV_D]; }
  [[nodiscard]] double const* inv_determinant() const noexcept { return matrices_[INV_D]; }

  [[nodiscard]] double*       lower_upper()       noexcept { return matrices_[LU]; }
  [[nodiscard]] double const* lower_upper() const noexcept { return matrices_[LU]; }

  [[nodiscard]] int*       pivot()       noexcept { return int_vec_[PIVOT]; }
  [[nodiscard]] int const* pivot() const noexcept { return int_vec_[PIVOT]; }

  Ptr3D<      int> n_vector()       noexcept { return {int_vec_[N_X], int_vec_[N_Y], int_vec_[N_Z]}; }
  Ptr3D<const int> n_vector() const noexcept { return {int_vec_[N_X], int_vec_[N_Y], int_vec_[N_Z]}; }

  Ptr3D<      double> k_vector()       noexcept { return {double_vec_[K_X], double_vec_[K_Y], double_vec_[K_Z]}; }
  Ptr3D<const double> k_vector() const noexcept { return {double_vec_[K_X], double_vec_[K_Y], double_vec_[K_Z]}; }

  [[nodiscard]] double*       solution()       noexcept { return double_vec_[SOLUTION]; }
  [[nodiscard]] double const* solution() const noexcept { return double_vec_[SOLUTION]; }

  [[nodiscard]] double*       rhs()       noexcept { return double_vec_[RHS]; }
  [[nodiscard]] double const* rhs() const noexcept { return double_vec_[RHS]; }

  [[nodiscard]] double*       sin_cache()       noexcept { return trig_cache_[SIN_CACHE]; }
  [[nodiscard]] double const* sin_cache() const noexcept { return trig_cache_[SIN_CACHE]; }

  [[nodiscard]] double*       cos_cache()       noexcept { return trig_cache_[COS_CACHE]; }
  [[nodiscard]] double const* cos_cache() const noexcept { return trig_cache_[COS_CACHE]; }

  void save_trig_row(std::size_t particle) noexcept;
  void restore_trig_row(std::size_t particle) noexcept;

  void update_trig_cache(std::size_t particle, const Particles& particles) noexcept;

  double log_abs_det(const Particles& particles);

  double* build_row(std::size_t particle) noexcept;

  [[nodiscard]] double determinant_ratio(
    std::size_t particle,
    const double* new_row
  ) const noexcept;

  void accept_move(std::size_t particle, const double* new_row, double ratio) noexcept;

  void add_derivatives(
    double* RESTRICT grad_x,
    double* RESTRICT grad_y,
    double* RESTRICT grad_z,
    double* RESTRICT laplacian
  ) const noexcept;

private:
  [[nodiscard]] double* new_row() noexcept { return double_vec_[NEW_ROW]; }
  [[nodiscard]] double* inv_d_col() noexcept { return double_vec_[INV_D_COL]; }
};