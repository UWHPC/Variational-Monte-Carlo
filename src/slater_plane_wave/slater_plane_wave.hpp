#pragma once

#include "../particles/particles.hpp"
#include "../utilities/aligned_soa.hpp"
#include "../utilities/memory.hpp"

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <numbers>
#include <vector>

class SlaterPlaneWave {
private:
    std::size_t num_orbitals_;
    std::size_t num_unique_k_;
    std::size_t matrix_size_;
    double box_length_;

    // Per-orbital mapping: which unique k-vector does orbital j use
    std::vector<std::size_t> orbital_k_index_;

    // Per-orbital type: 0 = cos, 1 = sin
    std::vector<std::uint8_t> orbital_type_;

    // Int vectors (sized to num_orbitals_ for pivot, num_unique_k_ for n-vectors):
    enum PivotIndex : std::size_t { N_X_, N_Y_, N_Z_, PIVOT_, NUM_INT_VECTORS_ };
    AlignedSoA<int> int_vectors_;

    // Double vectors (k-vectors sized to num_unique_k_; rest to num_orbitals_):
    enum VectorIndex : std::size_t { K_X_, K_Y_, K_Z_, RHS_, SOLUTION_, NEW_ROW_, INV_D_COL_, NUM_DOUBLE_VECTORS_ };
    AlignedSoA<double> double_vectors_;

    // Trig cache array - size of N * num_k-vectors
    enum TrigIndex : std::size_t { SIN_CACHE_, COS_CACHE_, NUM_TRIG_ARRAYS_ };
    AlignedSoA<double> trig_cache_;

    enum ScratchTrigIndex : std::size_t { SIN_SAVED_, COS_SAVED_, NUM_SCRATCH_TRIG_ };
    AlignedSoA<double> trig_scratch_;

    // All matrices (sized to num_orbitals_^2):
    enum MatrixIndex : std::size_t { D_, INV_D_, LU_, NUM_MATRIX_ };
    AlignedSoA<double> matrices_;

public:
    explicit SlaterPlaneWave(const Particles& particles, double box_length);

    // Getters:
    // Num orbitals - N (num_particles)
    [[nodiscard]] std::size_t num_orbitals_get() const noexcept { return num_orbitals_; }

    // Number of unique k-vectors (after +-n deduplication)
    [[nodiscard]] std::size_t num_unique_k_get() const noexcept { return num_unique_k_; }
    [[nodiscard]] std::size_t& num_unique_k_set() noexcept { return num_unique_k_; }

    // Matrix size - N^2
    [[nodiscard]] std::size_t matrix_size_get() const noexcept { return matrix_size_; }

    // Box length - L
    [[nodiscard]] double box_length_get() const noexcept { return box_length_; }

    // Per-orbital k-vector index
    [[nodiscard]] std::vector<std::size_t>& orbital_k_index_get() noexcept { return orbital_k_index_; }
    [[nodiscard]] const std::vector<std::size_t>& orbital_k_index_get() const noexcept { return orbital_k_index_; }

    // Per-orbital type (0=cos, 1=sin)
    [[nodiscard]] std::vector<std::uint8_t>& orbital_type_get() noexcept { return orbital_type_; }
    [[nodiscard]] const std::vector<std::uint8_t>& orbital_type_get() const noexcept { return orbital_type_; }

    // Det. matrix
    [[nodiscard]] double* determinant_get() noexcept { return matrices_[D_]; }
    [[nodiscard]] double const* determinant_get() const noexcept { return matrices_[D_]; }

    // Inv. det. matrix
    [[nodiscard]] double* inv_determinant_get() noexcept { return matrices_[INV_D_]; }
    [[nodiscard]] double const* inv_determinant_get() const noexcept { return matrices_[INV_D_]; }

    // Lower upper matrix
    [[nodiscard]] double* lower_upper_get() noexcept { return matrices_[LU_]; }
    [[nodiscard]] double const* lower_upper_get() const noexcept { return matrices_[LU_]; }

    // Pivot vector
    [[nodiscard]] int* pivot_get() noexcept { return int_vectors_[PIVOT_]; }
    [[nodiscard]] int const* pivot_get() const noexcept { return int_vectors_[PIVOT_]; }

    // X component of n (length = num_unique_k_)
    [[nodiscard]] int* n_vector_x_get() noexcept { return int_vectors_[N_X_]; }
    [[nodiscard]] int const* n_vector_x_get() const noexcept { return int_vectors_[N_X_]; }

    // Y component of n
    [[nodiscard]] int* n_vector_y_get() noexcept { return int_vectors_[N_Y_]; }
    [[nodiscard]] int const* n_vector_y_get() const noexcept { return int_vectors_[N_Y_]; }

    // Z component of n
    [[nodiscard]] int* n_vector_z_get() noexcept { return int_vectors_[N_Z_]; }
    [[nodiscard]] int const* n_vector_z_get() const noexcept { return int_vectors_[N_Z_]; }

    // Solution vector
    [[nodiscard]] double* solution_get() noexcept { return double_vectors_[SOLUTION_]; }
    [[nodiscard]] double const* solution_get() const noexcept { return double_vectors_[SOLUTION_]; }

    // RHS vector
    [[nodiscard]] double* rhs_get() noexcept { return double_vectors_[RHS_]; }
    [[nodiscard]] double const* rhs_get() const noexcept { return double_vectors_[RHS_]; }

    // X component of k (length = num_unique_k_)
    [[nodiscard]] double* k_vector_x_get() noexcept { return double_vectors_[K_X_]; }
    [[nodiscard]] double const* k_vector_x_get() const noexcept { return double_vectors_[K_X_]; }

    // Y component of k
    [[nodiscard]] double* k_vector_y_get() noexcept { return double_vectors_[K_Y_]; }
    [[nodiscard]] double const* k_vector_y_get() const noexcept { return double_vectors_[K_Y_]; }

    // Z component of k
    [[nodiscard]] double* k_vector_z_get() noexcept { return double_vectors_[K_Z_]; }
    [[nodiscard]] double const* k_vector_z_get() const noexcept { return double_vectors_[K_Z_]; }

    [[nodiscard]] double* sin_cache_get() noexcept { return trig_cache_[SIN_CACHE_]; }
    [[nodiscard]] double const* sin_cache_get() const noexcept { return trig_cache_[SIN_CACHE_]; }

    [[nodiscard]] double* cos_cache_get() noexcept { return trig_cache_[COS_CACHE_]; }
    [[nodiscard]] double const* cos_cache_get() const noexcept { return trig_cache_[COS_CACHE_]; }

    void save_trig_row(std::size_t particle) noexcept;
    void restore_trig_row(std::size_t particle) noexcept;

    void update_trig_cache(std::size_t particle, const Particles& particles) noexcept;

    // Computes log|det(D)| via full LU - use for initialization only.
    double log_abs_det(const Particles& particles);

    // Builds the new Slater row for a moved particle into the internal NEW_ROW_ buffer.
    // Returns a pointer to the internal buffer (valid until next call).
    double* build_row(std::size_t particle) noexcept;

    // Call after build_row. O(N).
    [[nodiscard]] double determinant_ratio(std::size_t particle, const double* new_row) const noexcept;

    // Applies the Sherman-Morrison inverse update and patches D.
    // Call only after accepting a move. O(N^2).
    void accept_move(std::size_t particle, const double* new_row, double ratio) noexcept;

    // Accumulates Slater contributions into grad/lap (length = stride/at least N).
    void add_derivatives(double* RESTRICT grad_x, double* RESTRICT grad_y, double* RESTRICT grad_z,
                         double* RESTRICT laplacian) const noexcept;

private:
    // Scratch buffer: new Slater row for moved particle
    [[nodiscard]] double* new_row_get() noexcept { return double_vectors_[NEW_ROW_]; }

    // Scratch buffer: column of D_inv used by Sherman-Morrison
    [[nodiscard]] double* inv_d_col_get() noexcept { return double_vectors_[INV_D_COL_]; }
};
