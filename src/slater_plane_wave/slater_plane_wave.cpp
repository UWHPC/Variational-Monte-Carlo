#include "slater_plane_wave.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <new>      // std::bad_alloc

namespace {

inline std::size_t roundUpToSimd(std::size_t n) noexcept {
    constexpr std::size_t doublesPerAlignment = SIMD_BYTES / sizeof(double);
    return (n + doublesPerAlignment - 1) & ~(doublesPerAlignment - 1);
}

// @brief helper function to convert i-jth indices -> n
// @param stride is the difference between the row and the column.
// @return the appropriate i-th row j-th column as a size_t.
inline std::size_t index(std::size_t row, std::size_t col, std::size_t stride) noexcept {
    return row * stride + col;
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
        double maxAbs{ std::abs(LU[index(col, col, N)])};

        for(std::size_t row = col + 1; row < N; ++row) {
            const double value = std::abs(LU[index(row, col, N)]);
            if (value > maxAbs) { 
                maxAbs = value; 
                pivotRow = row;
            };
        }

        // max abs = 0.0 implies the pivot column is 0 & det = 0.
        if (maxAbs == 0.0) continue;

        if (pivotRow != col) {
            for (std::size_t col2 = 0; col2 < N; ++col2) {
                std::swap(LU[index(col, col2, N)], LU[index(pivotRow, col2, N)]);
            }
            std::swap(piv[col], piv[pivotRow]);
            ++swapCount;
        }

        // eliminate
        const double pivotValue = LU[index(col, col, N)];
        for (std::size_t row = col + 1; row < N; ++row) {
            LU[index(row, col, N)] /= pivotValue; // L (i,k)
            const double multiplier = LU[index(row, col, N)];
            for (std::size_t col2 = col + 1; col2 < N; ++col2) {   
                LU[index(row, col2, N)] -= multiplier * LU[index(col, col2, N)];
            }
        }
    }

    return swapCount;
}

/*
see https://www.geeksforgeeks.org/dsa/doolittle-algorithm-lu-decomposition/
solve (P^-1)LU x = b. given combined LU and pivot permutation piv.
piv encodes the row permutation applied during LU so that
we first permute b: y = P b, then solve L z = y, then U x = z.
*/
void lu_solve(const double* LU, const double* piv,
              const double* b, double* x, std::size_t N)
{   
    // Apply permutation: x = Pb
    // store y in x temporarily
    for (std::size_t row = 0; row < N; ++row) {
        const std::size_t permRow{ static_cast<std::size_t>(piv[row]) };
        x[row] = b[permRow];
    }

    // forward solve: ly = Pb (L has implicit on diagonal)
    for (std::size_t row = 0; row < N; ++row) {
        double sum = x[row];
        for (std::size_t col = 0; col < row; ++col) {
            sum -= LU[index(row, col, N)] * x[col];
        }
        x[row] = sum;
    }

    // backward solve: Ux = y
    for (std::size_t rev = 0; rev < N; ++rev) {
        const std::size_t row = N - 1 - rev;
        double sum = x[row];
        for (std::size_t col = row + 1; col < N; ++col) {
            sum -= LU[index(row, col, N)] * x[col];
        }
        x[row] = sum / LU[index(row, row, N)];
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
    const double* posX = particles.posX();
    const double* posY = particles.posY();
    const double* posZ = particles.posZ();

    // 1.
    for (std::size_t particle = 0; particle < N; ++particle) {
        const double x = posX[particle];
        const double y = posY[particle];
        const double z = posZ[particle];

        for (std::size_t orbital = 0; orbital < N; ++orbital) {
            const double kdotr = kx_[orbital] * x + ky_[orbital] * y + kz_[orbital] * z;
            D_[index(particle, orbital, N)] = std::cos(kdotr);
        }
    }

    // copy D_ -> LU
    std::copy_n(D_, N * N, LU_);
    // cast void because we dont need the # of swaps.
    (void)lu_decompose(LU_, piv_, N);

    // log|det(D)| = sum log|U(ii)|
    double logAbsDet = 0.0;
    for (std::size_t diag = 0; diag < N; ++diag) {
        const double Uii = LU_[index(diag, diag, N)];
        const double absUii  = std::abs(Uii);
        if (absUii == 0.0 || !std::isfinite(absUii)) {
            return -std::numeric_limits<double>::infinity();
        }
        logAbsDet += std::log(absUii);
    }
    // TO DO CHANGE THESE VECTORS 
    std::vector<double> rhs(N, 0.0);
    std::vector<double> sol(N, 0.0);

    for (std::size_t col = 0; col < N; ++col) {
        std::fill(rhs.begin(), rhs.end(), 0.0);
        rhs[col] = 1.0;
        lu_solve(LU_, piv_, rhs.data(), sol.data(), N);

        for (std::size_t row = 0; row < N; ++row) {
            invD_[index(row, col, N)] = sol[row];
        }
    }

    return logAbsDet;
}