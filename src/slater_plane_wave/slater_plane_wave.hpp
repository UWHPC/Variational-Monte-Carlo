#pragma once

#include "../particles/particles.hpp"
#include "../utilities/aligned_soa.hpp"
#include "../utilities/memory.hpp"

#include <cstddef>
#include <cstdlib>

class SlaterPlaneWave {
private:
    std::size_t num_orbitals_;
    std::size_t matrix_size_;
    double box_length_;

    // All vectors:
    enum VectorIndex : std::size_t { K_X_, K_Y_, K_Z_, NUM_K_VECTORS_ };
    AlignedSoA<double> k_vector_;

    // Pivot:
    enum PivotIndex : std::size_t { N_X_, N_Y_, N_Z_, PIVOT_, NUM_INT_VECTORS_ };
    AlignedSoA<int> int_vectors_;

    // All matrices:
    enum MatrixIndex : std::size_t { D_, INV_D_, LU_, NUM_MATRIX_ };
    AlignedSoA<double> matrices_;

public:
    explicit SlaterPlaneWave(std::size_t N, double L)
        : num_orbitals_{N}, matrix_size_{N * N}, box_length_{L}, k_vector_{N, NUM_K_VECTORS_},
          int_vectors_{N, NUM_INT_VECTORS_}, matrices_{N * N, NUM_MATRIX_} {};

    // Getters:
    // Num orbitals - N
    [[nodiscard]] std::size_t num_orbitals_get() const noexcept { return num_orbitals_; }

    // Matrix size - N^2
    [[nodiscard]] std::size_t matrix_size_get() const noexcept { return matrix_size_; }

    // Box length - L
    [[nodiscard]] double box_length_get() const noexcept { return box_length_; }

    // Det. matrix
    [[nodiscard]] double* determinant_get() noexcept { return matrices_[D_]; }
    [[nodiscard]] double const* determinant_get() const noexcept { return matrices_[D_]; }

    // Inv. det. matrix
    [[nodiscard]] double* inv_determinant_get() noexcept { return matrices_[INV_D_]; }
    [[nodiscard]] double const* inv_determinant_get() const noexcept { return matrices_[INV_D_]; }

    // Lower upper matrix
    [[nodiscard]] double* lower_upper_get() noexcept { return matrices_[LU_]; }
    [[nodiscard]] double const* lower_upper_get() const noexcept { return matrices_[LU_]; }

    // Pivot matrix
    [[nodiscard]] int* pivot_get() noexcept { return int_vectors_[PIVOT_]; }
    [[nodiscard]] int const* pivot_get() const noexcept { return int_vectors_[PIVOT_]; }

    // X component of n
    [[nodiscard]] int* n_vector_x_get() noexcept { return int_vectors_[N_X_]; }
    [[nodiscard]] int const* n_vector_x_get() const noexcept { return int_vectors_[N_X_]; }

    // Y component of n
    [[nodiscard]] int* n_vector_y_get() noexcept { return int_vectors_[N_Y_]; }
    [[nodiscard]] int const* n_vector_y_get() const noexcept { return int_vectors_[N_Y_]; }

    // Z component of n
    [[nodiscard]] int* n_vector_z_get() noexcept { return int_vectors_[N_Z_]; }
    [[nodiscard]] int const* n_vector_z_get() const noexcept { return int_vectors_[N_Z_]; }

    // X component of k
    [[nodiscard]] double* k_vector_x_get() noexcept { return k_vector_[K_X_]; }
    [[nodiscard]] double const* k_vector_x_get() const noexcept { return k_vector_[K_X_]; }

    // Y component of k
    [[nodiscard]] double* k_vector_y_get() noexcept { return k_vector_[K_Y_]; }
    [[nodiscard]] double const* k_vector_y_get() const noexcept { return k_vector_[K_Y_]; }

    // Z component of k
    [[nodiscard]] double* k_vector_z_get() noexcept { return k_vector_[K_Z_]; }
    [[nodiscard]] double const* k_vector_z_get() const noexcept { return k_vector_[K_Z_]; }

    // Computes log|det(D)| and updates internal cached inverse/LU.
    [[nodiscard]] double log_abs_det(const Particles& particles);

    // Accumulates Slater contributions into grad/lap (length = stride/at least N).
    void add_derivatives(const Particles& particles, double* RESTRICT grad_x, double* RESTRICT grad_y,
                         double* RESTRICT grad_z, double* RESTRICT laplacian) const noexcept;
};
