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
    double boxLength_;

    // All sub-arrays:
    static constexpr std::size_t D_{0}, INV_D_{1};
    static constexpr std::size_t LU_{2};
    static constexpr std::size_t PIVOT_{3};
    static constexpr std::size_t K_X_{4}, K_Y_{5}, K_Z_{6};

    // Number of sub-arrays:
    static constexpr std::size_t NUM_SUB_ARRAYS{7};

    // Memory data:
    AlignedSoA slaterPlaneWaveData_;

    // static constexpr std::size_t numVectorComponents_{8};  // Number of components
    static constexpr std::size_t alignmentBytes_{SIMD_BYTES}; // SIMD byte alignment

public:
    explicit SlaterPlaneWave(std::size_t N, double L)
        : numOrbitals_{N}, boxLength_{L}, slaterPlaneWaveData_{N, NUM_SUB_ARRAYS} {};

    [[nodiscard]] std::size_t numOrbitals() const noexcept { return numOrbitals_; }
    [[nodiscard]] double boxLength() const noexcept { return boxLength_; }

    // --- Matrix buffers (row-major N x N stored in padded block) ---
    [[nodiscard]] double* determinant() noexcept { return slaterPlaneWaveData_[D_]; }
    [[nodiscard]] double* invDeterminant() noexcept { return slaterPlaneWaveData_[INV_D_]; }
    [[nodiscard]] double* lowerUpper() noexcept { return slaterPlaneWaveData_[LU_]; }

    [[nodiscard]] double const* determinant() const noexcept { return slaterPlaneWaveData_[D_]; }
    [[nodiscard]] double const* invDeterminant() const noexcept { return slaterPlaneWaveData_[INV_D_]; }
    [[nodiscard]] double const* lowerUpper() const noexcept { return slaterPlaneWaveData_[LU_]; }

    // --- Pivot buffer ---

    [[nodiscard]] double* pivot() noexcept { return slaterPlaneWaveData_[PIVOT_]; }             // MUT - pivot
    [[nodiscard]] double const* pivot() const noexcept { return slaterPlaneWaveData_[PIVOT_]; } // IMMUT - pivot

    // --- k-vector buffers (length N, padded to vecStride_) ---

    [[nodiscard]] double* kVectorX() noexcept { return slaterPlaneWaveData_[K_X_]; } // MUT - X component of k
    [[nodiscard]] double* kVectorY() noexcept { return slaterPlaneWaveData_[K_Y_]; } // MUT - Y component of k
    [[nodiscard]] double* kVectorZ() noexcept { return slaterPlaneWaveData_[K_Z_]; } // MUT - Z component of k

    // IMMUT - X component of k
    [[nodiscard]] double const* kVectorX() const noexcept { return slaterPlaneWaveData_[K_X_]; }

    // IMMUT - Y component of k
    [[nodiscard]] double const* kVectorY() const noexcept { return slaterPlaneWaveData_[K_Y_]; }

    // IMMUT - Z component of k
    [[nodiscard]] double const* kVectorZ() const noexcept { return slaterPlaneWaveData_[K_Z_]; }

    // Computes log|det(D)| and updates internal cached inverse/LU.
    [[nodiscard]] double logAbsDet(const Particles& particles, const PeriodicBoundaryCondition& pbc);

    // Accumulates Slater contributions into grad/lap (length = stride/at least N).
    void addDerivatives(const Particles& particles, const PeriodicBoundaryCondition& pbc, double* gradX, double* gradY,
                        double* gradZ, double* la) const noexcept;
};
