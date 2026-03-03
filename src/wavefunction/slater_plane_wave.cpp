#include "slater_plane_wave.hpp"

#include <algorithm>
#include <cstddef>
#include <new>      // std::bad_alloc

namespace {

inline std::size_t roundUpToSimd(std::size_t n) noexcept {
    constexpr std::size_t doublesPerAlignment = SIMD_BYTES / sizeof(double);
    return (n + doublesPerAlignment - 1) & ~(doublesPerAlignment - 1);
}

inline std::size_t id(std::size_t i, std::size_t j, std::size_t N) noexcept {
    return i * N + j;
}

/*
In-Place LU with partial pivoting in LU (size N*N, row-major)
pivot is length >= N, which stores row permutation info as doubles,
return numbers of row swaps (parity info, if you need det sign).
*/
int lu_decompose(double* LU, double* piv, std::size_t N) {
    for(std::size_t i = 0; i < N; ++i) piv[i] = static_cast<double>(i);
    int swapCount{};
    for (std::size_t k = 0; k < N; ++k) {
        // we chose a pivot row p
        std::size_t p{ k };
        double maxAbs{ std::abs(LU[id(k,k, N)])};
        for(std::size_t i = k + 1; i < N; ++i) {
            const double v = std::abs(LU[id(k, k, N)]);
            if (v > maxAbs) { maxAbs = v; p = i; };
        }
        if (maxAbs == 0.0) continue; // singular column det will be 0 so just continue
        if (p != k) {
            // we swap rows p and k in LU
            for (std::size_t j = 0; j < N; ++j)
            {
                std::swap(LU[id(k, j, N)], LU[id(p, j, N)]);
            }
            std::swap(piv[k], piv[p]);
            ++swapCount;
        }

        const double Akk = LU[id(k, k, N)];
        for (std::size_t i = k + 1; i < N; ++i) {
            LU[id(i, k, N)] /= Akk; // L (i,k)
            const double Lik = LU[id(i, k, N)];
            for (std::size_t j = k + 1; j < N; ++j) {
                LU[id(i, j, N)] -= Lik * LU[id(k, j, N)];
            }
        }
    }

    return swapCount;
}

/*
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
}