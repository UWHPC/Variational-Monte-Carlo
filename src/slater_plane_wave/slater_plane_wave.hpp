#pragma once

#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "../utilities/aligned_soa.hpp"
#include "../utilities/memory.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <stdexcept>

class SlaterPlaneWave {
private:
    std::size_t numOrbitals_;
    std::size_t matrixSize_;
    double boxLength_;

    // All vectors:
    enum VectorIndex : std::size_t { K_X_, K_Y_, K_Z_, NUM_VECTORS_ };
    AlignedSoA<double> vectors_;

    // Pivot:
    enum PivotIndex : std::size_t { PIVOT_, NUM_PIVOT_ };
    AlignedSoA<int> pivot_;

    // All matrices:
    enum MatrixIndex : std::size_t { D_, INV_D_, LU_, NUM_MATRIX_ };
    AlignedSoA<double> matrices_;

public:
    explicit SlaterPlaneWave(std::size_t N, double L)
        : numOrbitals_{N}, matrixSize_{N * N}, boxLength_{L}, vectors_{N, NUM_VECTORS_}, pivot_{N, NUM_PIVOT_},
          matrices_{N * N, NUM_MATRIX_} {};

    // Getters:
    // Num orbitals - N
    [[nodiscard]] std::size_t numOrbitals() const noexcept { return numOrbitals_; }

    // Matrix size - N^2
    [[nodiscard]] std::size_t matrixSize() const noexcept { return matrixSize_; }

    // Box length - L
    [[nodiscard]] double boxLength() const noexcept { return boxLength_; }

    // Det. matrix
    [[nodiscard]] double* determinant() noexcept { return matrices_[D_]; }
    [[nodiscard]] double const* determinant() const noexcept { return matrices_[D_]; }

    // Inv. det. matrix
    [[nodiscard]] double* invDeterminant() noexcept { return matrices_[INV_D_]; }
    [[nodiscard]] double const* invDeterminant() const noexcept { return matrices_[INV_D_]; }

    // Lower upper matrix
    [[nodiscard]] double* lowerUpper() noexcept { return matrices_[LU_]; }
    [[nodiscard]] double const* lowerUpper() const noexcept { return matrices_[LU_]; }

    // Pivot matrix
    [[nodiscard]] int* pivot() noexcept { return pivot_[PIVOT_]; }
    [[nodiscard]] int const* pivot() const noexcept { return pivot_[PIVOT_]; }

    // X component of k
    [[nodiscard]] double* kVectorX() noexcept { return vectors_[K_X_]; }
    [[nodiscard]] double const* kVectorX() const noexcept { return vectors_[K_X_]; }

    // Y component of k
    [[nodiscard]] double* kVectorY() noexcept { return vectors_[K_Y_]; }
    [[nodiscard]] double const* kVectorY() const noexcept { return vectors_[K_Y_]; }

    // Z component of k
    [[nodiscard]] double* kVectorZ() noexcept { return vectors_[K_Z_]; }
    [[nodiscard]] double const* kVectorZ() const noexcept { return vectors_[K_Z_]; }

    // Computes log|det(D)| and updates internal cached inverse/LU.
    [[nodiscard]] double logAbsDet(const Particles& particles, const PeriodicBoundaryCondition& pbc);

    // Accumulates Slater contributions into grad/lap (length = stride/at least N).
    void addDerivatives(const Particles& particles, const PeriodicBoundaryCondition& pbc, double* RESTRICT gradX,
                        double* RESTRICT gradY, double* RESTRICT gradZ, double* RESTRICT lap) const noexcept;
};
