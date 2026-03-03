#include "slater_plane_wave.hpp"

#include <algorithm>
#include <cstddef>
#include <new>      // std::bad_alloc

namespace {
    inline std::size_t roundUpToSimd(std::size_t n) noexcept {
        constexpr std::size_t doublesPerAlignment = SIMD_BYTES / sizeof(double);
        return (n + doublesPerAlignment - 1) & ~(doublesPerAlignment - 1);
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