#pragma once

#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "../memory/memory.hpp"

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <cstddef>
#include <stdexcept>

class SlaterPlaneWave {
private:
    std::size_t N_{}; // number of orbitals (initialized to 0)
    double L_{};      // box length (initialized to 0.0)

    // static constexpr std::size_t numVectorComponents_{8};                         // Number of components
    static constexpr std::size_t alignmentBytes_{SIMD_BYTES};                     // SIMD byte alignment

    std::unique_ptr<double[], AlignedDeleter> memoryBlock_; // Memory block size

    std::size_t matStride_{}; // padded length for N*N arrays
    std::size_t vecStride_{}; // padded length for N arrays

    double* D_{nullptr};
    double* invD_{nullptr};
    double* LU_{nullptr};
    double* piv_{nullptr}; //pivot indices stored as double
    double* kx_{nullptr};
    double* ky_{nullptr};
    double* kz_{nullptr};


public:
    explicit SlaterPlaneWave(std::size_t N, double L);
    [[nodiscard]] std::size_t N() const noexcept { return N_; }
    [[nodiscard]] double L() const noexcept { return L_; }
    [[nodiscard]] std::size_t matStride() const noexcept { return matStride_; }
    [[nodiscard]] std::size_t vecStride() const noexcept { return vecStride_; }

    // --- Matrix buffers (row-major N x N stored in padded block) ---
    [[nodiscard]] double* D() noexcept { return D_; }
    [[nodiscard]] double* invD() noexcept { return invD_; }
    [[nodiscard]] double* LU() noexcept { return LU_; }

    [[nodiscard]] double const* D() const noexcept { return D_; }
    [[nodiscard]] double const* invD() const noexcept { return invD_; }
    [[nodiscard]] double const* LU() const noexcept { return LU_; }

    // --- Pivot buffer ---

    [[nodiscard]] double* piv() noexcept { return piv_; }             // MUT - pivot
    [[nodiscard]] double const* piv() const noexcept { return piv_; } // IMMUT - pivot

    // --- k-vector buffers (length N, padded to vecStride_) ---

    [[nodiscard]] double* kx() noexcept { return kx_; } // MUT - X component of k
    [[nodiscard]] double* ky() noexcept { return ky_; } // MUT - Y component of k
    [[nodiscard]] double* kz() noexcept { return kz_; } // MUT - Z component of k

    [[nodiscard]] double const* kx() const noexcept { return kx_; } // IMMUT - X component of k
    [[nodiscard]] double const* ky() const noexcept { return ky_; } // IMMUT - Y component of k
    [[nodiscard]] double const* kz() const noexcept { return kz_; } // IMMUT - Z component of k

    // Computes log|det(D)| and updates internal cached inverse/LU.
    [[nodiscard]] double logAbsDet(const Particles& particles, const periodicBoundaryCondition& pbc);

    // Accumulates Slater contributions into grad/lap (length = stride/at least N).
    void addDerivatives(
        const Particles& particles,
        const PeriodicBoundaryCondition& pbc,
        double* gradX,
        double* gradY,
        double* gradZ,
        double* la
    ) const noexcept;
};