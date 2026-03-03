#include "slater_plane_wave.hpp"

#include <cmath>
#include <algorithm>
#include <cstddef>
#include <new>      // std::bad_alloc

namespace {

inline std::size_t roundUpToSimd(std::size_t n) noexcept {
    constexpr std::size_t doublesPerAlignment = SIMD_BYTES / sizeof(double);
    return (n + doublesPerAlignment - 1) & ~(doublesPerAlignment - 1);
}

// matrix is stored as a single array to access the ij-th value u do i * N + j.
inline std::size_t id(std::size_t i, std::size_t j, std::size_t N) noexcept {
    return i * N + j;
}

/*
see https://www.geeksforgeeks.org/dsa/doolittle-algorithm-lu-decomposition/
In-Place LU with partial pivoting in LU (size N*N, row-major)
pivot is length >= N, which stores row permutation info as doubles,
return numbers of row swaps (parity info, if you need det sign).
*/
int lu_decompose(double* LU, double* piv, std::size_t N) {
    // Track row swaps
    int swapCount{};

    for (std::size_t row = 0; row < N; ++row) piv[row] = static_cast<double>(row);
    for (std::size_t col = 0; col < N; ++col) {
        // Pivot selection
        // Find row >= col maximizing |LU(row, col)|
        std::size_t pivotRow{ col };
        double maxAbs{ std::abs(LU[id(col, col, N)])};

        for(std::size_t row = col + 1; row < N; ++row) {
            const double value = std::abs(LU[id(row, col, N)]);
            if (value > maxAbs) { 
                maxAbs = value; 
                pivotRow = row;
            };
        }

        // max abs = 0.0 implies the pivot column is 0 & det = 0.
        if (maxAbs == 0.0) continue;

        if (pivotRow != col) {
            for (std::size_t col2 = 0; col2 < N; ++col2) {
                std::swap(LU[id(col, col2, N)], LU[id(pivotRow, col2, N)]);
            }
            std::swap(piv[col], piv[pivotRow]);
            ++swapCount;
        }

        // eliminate
        const double pivotValue = LU[id(col, col, N)];
        for (std::size_t row = col + 1; row < N; ++row) {
            LU[id(row, col, N)] /= pivotValue; // L (i,k)
            const double multiplier = LU[id(row, col, N)];
            for (std::size_t col2 = col + 1; col2 < N; ++col2) {   
                LU[id(row, col2, N)] -= multiplier * LU[id(col, col2, N)];
            }
        }
    }

    return swapCount;
}

/*
see https://www.geeksforgeeks.org/dsa/doolittle-algorithm-lu-decomposition/
Solving LU*x = b using LU, piv permutation.
We use x as an output buffer (length N), b as input (length N).
*/
void lu_solve(const double* LU, const double* piv,
              const double* b, double* x, std::size_t N)
{   
    // Apply permutation: x = Pb
    for (std::size_t i = 0; i < N; ++i) {
        const std::size_t pi{ static_cast<std::size_t>(piv[i]) };
        x[i] = b[pi];
    }

    // forward solve: ly = Pb (L has implicit on diagonal)
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            x[i] -= LU[id(i, j, N)] * x[j];
        }
    }

    // backward solve: Ux = y
    for (std::size_t ii = 0; ii < N; ++ii) {
        const std::size_t i = N - 1 - ii;
        for (std::size_t j = i + 1; j < N; ++j) {
            x[i] -= LU[id(i, j, N)] * x[j];
        }
        x[i] /= LU[id(i, i, N)];
    }
}

}

SlaterPlaneWave::SlaterPlaneWave(std::size_t N, double L)
: N_{N}
, L_{L}
{
    vecStride_ = { roundUpToSimd(N_) };
    matStride_ = { roundUpToSimd(N_ * N_) };

    // Layout: [ D | invD | LU | piv | kx | ky | kz ]
    const std::size_t totalDoubles {
        3 * matStride_ +      // D, invD, LU
        4 * vecStride_        // piv, kx, ky, kz
    };      
    const std::size_t totalBytes = totalDoubles * sizeof(double);

    double* ptr = static_cast<double*>(alignedAlloc(alignmentBytes_, totalBytes));
    if (!ptr) throw std::bad_alloc();

    std::fill_n(ptr, totalDoubles, 0.0);
    memoryBlock_.reset(ptr);

    double* cur = memoryBlock_.get();

    D_    = cur; cur += matStride_;
    invD_ = cur; cur += matStride_;
    LU_   = cur; cur += matStride_;

    piv_  = cur; cur += vecStride_;
    kx_   = cur; cur += vecStride_;
    ky_   = cur; cur += vecStride_;
    kz_   = cur; cur += vecStride_;

    if (cur != memoryBlock_.get() + totalDoubles) throw std::runtime_error("Bad slice");
}

double SlaterPlaneWave::logAbsDet(const Particles& particles, const PeriodicBoundaryCondition& pbc) 
{
    const std::size_t N = N_;
    const double* x = particles.posX();
    const double* y = particles.posY();
    const double* z = particles.posX();

    // 1.
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            const double dot = kx_[j] * x[i] + ky_[j] * y[i] + kz_[j] * z[i];
            D_[id(i, j, N)] = std::cos(dot);
        }
    }

    // 2. copy D_ -> LU
    std::copy_n(D_, N * N, LU_);
}