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
    enum VectorIndex : std::size_t { K_X_, K_Y_, K_Z_, NUM_VECTORS_ };
    AlignedSoA<double> vectors_;

    // Pivot:
    enum PivotIndex : std::size_t { PIVOT_, NUM_PIVOT_ };
    AlignedSoA<std::size_t> pivot_;

    // All matrices:
    enum MatrixIndex : std::size_t { D_, INV_D_, LU_, NUM_MATRIX_ };
    AlignedSoA<double> matrices_;

public:
    explicit SlaterPlaneWave(std::size_t N, double L)
        : num_orbitals_{N}, matrix_size_{N * N}, box_length_{L}, vectors_{N, NUM_VECTORS_}, pivot_{N, NUM_PIVOT_},
          matrices_{N * N, NUM_MATRIX_} {};

    // Getters:
    // Num orbitals - N
    [[nodiscard]] std::size_t num_orbitals_ptr() const noexcept { return num_orbitals_; }

    // Matrix size - N^2
    [[nodiscard]] std::size_t matrix_size_ptr() const noexcept { return matrix_size_; }

    // Box length - L
    [[nodiscard]] double box_length_ptr() const noexcept { return box_length_; }

    // Det. matrix
    [[nodiscard]] double* determinant_ptr() noexcept { return matrices_[D_]; }
    [[nodiscard]] double const* determinant_ptr() const noexcept { return matrices_[D_]; }

    // Inv. det. matrix
    [[nodiscard]] double* inv_determinant_ptr() noexcept { return matrices_[INV_D_]; }
    [[nodiscard]] double const* inv_determinant_ptr() const noexcept { return matrices_[INV_D_]; }

    // Lower upper matrix
    [[nodiscard]] double* lower_upper_ptr() noexcept { return matrices_[LU_]; }
    [[nodiscard]] double const* lower_ipper_ptr() const noexcept { return matrices_[LU_]; }

    // Pivot matrix
    [[nodiscard]] std::size_t* pivot_ptr() noexcept { return pivot_[PIVOT_]; }
    [[nodiscard]] std::size_t const* pivot_ptr() const noexcept { return pivot_[PIVOT_]; }

    // X component of k
    [[nodiscard]] double* k_vector_x_ptr() noexcept { return vectors_[K_X_]; }
    [[nodiscard]] double const* k_vector_x_ptr() const noexcept { return vectors_[K_X_]; }

    // Y component of k
    [[nodiscard]] double* k_vector_y_ptr() noexcept { return vectors_[K_Y_]; }
    [[nodiscard]] double const* k_vector_y_ptr() const noexcept { return vectors_[K_Y_]; }

    // Z component of k
    [[nodiscard]] double* k_vector_z_ptr() noexcept { return vectors_[K_Z_]; }
    [[nodiscard]] double const* k_vector_z_ptr() const noexcept { return vectors_[K_Z_]; }

    // Computes log|det(D)| and updates internal cached inverse/LU.
    [[nodiscard]] double log_abs_det(const Particles& particles);

    // Accumulates Slater contributions into grad/lap (length = stride/at least N).
    void add_derivatives(const Particles& particles, double* RESTRICT grad_x, double* RESTRICT grad_y,
                         double* RESTRICT grad_z, double* RESTRICT laplacian) const noexcept;
};
