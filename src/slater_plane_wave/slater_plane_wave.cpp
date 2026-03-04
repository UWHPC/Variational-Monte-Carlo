#include "slater_plane_wave.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <new>
#include <vector>

namespace {

// @brief helper function to convert i-jth indices -> n
// @param stride is the difference between the row and the column.
// @return the appropriate i-th row j-th column as a size_t.
inline std::size_t index(std::size_t row, std::size_t col, std::size_t stride) noexcept { return row * stride + col; }

/*
see https://www.geeksforgeeks.org/dsa/doolittle-algorithm-lu-decomposition/
In-Place LU with partial pivoting in LU (size N*N, row-major)
pivot is length >= N, which stores row permutation info as doubles,
return numbers of row swaps (parity info, if you need det sign).
*/
int lu_decompose(double* LU, double* piv, std::size_t N) {
    // Track row swaps
    int swapCount{};

    for (std::size_t row = 0; row < N; ++row)
        piv[row] = static_cast<double>(row);

    for (std::size_t col = 0; col < N; ++col) {
        // Pivot selection
        // Find row >= col maximizing |LU(row, col)|
        std::size_t pivotRow{col};
        double maxAbs{std::abs(LU[index(col, col, N)])};

        for (std::size_t row = col + 1; row < N; ++row) {
            const double value = std::abs(LU[index(row, col, N)]);
            if (value > maxAbs) {
                maxAbs = value;
                pivotRow = row;
            };
        }

        // max abs = 0.0 implies the pivot column is 0 & det = 0.
        if (maxAbs == 0.0)
            continue;

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
void lu_solve(const double* LU, const double* piv, const double* b, double* x, std::size_t N) {
    // Apply permutation: x = Pb
    // store y in x temporarily
    for (std::size_t row = 0; row < N; ++row) {
        const std::size_t permRow{static_cast<std::size_t>(piv[row])};
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

} // namespace

double SlaterPlaneWave::logAbsDet(const Particles& particles, const PeriodicBoundaryCondition& pbc) {
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
        const double absUii = std::abs(Uii);
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

void SlaterPlaneWave::addDerivatives(const Particles& particles, const PeriodicBoundaryCondition& pbc, double* gradX,
                                     double* gradY, double* gradZ, double* lap) const noexcept {
    const std::size_t N{numOrbitals()};

    const double* RESTRICT posX{particles.posX()};
    const double* RESTRICT posY{particles.posY()};
    const double* RESTRICT posZ{particles.posZ()};

    const double* RESTRICT k_x{kVectorX()};
    const double* RESTRICT k_y{kVectorY()};
    const double* RESTRICT k_z{kVectorZ()};

    const double* RESTRICT inv_det{invDeterminant()};

    for (std::size_t particle = 0; particle < N; ++particle) {
        const double p_x{posX[particle]};
        const double p_y{posY[particle]};
        const double p_z{posZ[particle]};

        double d_LogDet_dx{}, d_LogDet_dy{}, d_LogDet_dz{};

        // Σ_j (D^{-1})_{j,particle} * (∇^2 D_{particle,j})
        double laplace_det_term{};

        for (std::size_t orbital = 0; orbital < N; ++orbital) {
            const double kx_orbital{k_x[orbital]};
            const double ky_orbital{k_y[orbital]};
            const double kz_orbital{k_z[orbital]};

            const double k_dot_r{kx_orbital * p_x + ky_orbital * p_y + kz_orbital * p_z};
            const double cosine_term{std::cos(k_dot_r)};
            const double sine_term{std::sin(k_dot_r)};
            // D(particle,orbital) = cos(k·r)
            // ∇_particle D = -sin(k·r) * k
            const double dD_dx{-sine_term * kx_orbital};
            const double dD_dy{-sine_term * ky_orbital};
            const double dD_dz{-sine_term * kz_orbital};

            // ∇^2 D = -cos(k·r) * |k|^2
            const double kSquared{kx_orbital * kx_orbital + ky_orbital * ky_orbital + kz_orbital * kz_orbital};
            const double lapD{-cosine_term * kSquared};
            // Weight = (D^{-1})_{orbital,particle}
            // invD_ is row-major, so entry (row=orbital, col=particle):
            const double weight{inv_det[index(orbital, particle, N)]};

            d_LogDet_dx += weight * dD_dx;
            d_LogDet_dy += weight * dD_dy;
            d_LogDet_dz += weight * dD_dz;

            laplace_det_term += weight * lapD;
        }
        // ∇^2 log det D = Σ_j (D^{-1})_{j,i} ∇^2 D_{i,j}  -  ||∇ log det D||^2
        const double gradSquared{d_LogDet_dx * d_LogDet_dx + d_LogDet_dy * d_LogDet_dy + d_LogDet_dz * d_LogDet_dz};
        const double lapLogDet{laplace_det_term - gradSquared};

        gradX[particle] += d_LogDet_dx;
        gradY[particle] += d_LogDet_dy;
        gradZ[particle] += d_LogDet_dz;
        lap[particle] += lapLogDet;
    }
}