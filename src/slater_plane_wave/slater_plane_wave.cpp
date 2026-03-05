#include "slater_plane_wave.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace {

// @brief helper function to convert i-jth indices -> n
// @param stride is the difference between the row and the column.
// @return the appropriate i-th row j-th column as a size_t.
inline std::size_t index(std::size_t row, std::size_t col, std::size_t stride) noexcept { return row * stride + col; }

/*
see https://www.geeksforgeeks.org/dsa/doolittle-algorithm-lu-decomposition/
In-Place LU with partial pivoting in LU (size N*N, row-major)
pivot is length >= N storing row permutation indices
return numbers of row swaps (parity info, if you need det sign).
*/
int lowerUpperDecompose(double* lowerUpper, std::size_t* pivot, std::size_t N) {
    // Track row swaps
    int swapCount{};

    for (std::size_t row = 0; row < N; ++row)
        pivot[row] = row;

    for (std::size_t col = 0; col < N; ++col) {
        // Pivot selection
        // Find row >= col maximizing |LU(row, col)|
        std::size_t pivotRow{col};
        double maxAbs{std::abs(lowerUpper[index(col, col, N)])};

        for (std::size_t row = col + 1; row < N; ++row) {
            const double value = std::abs(lowerUpper[index(row, col, N)]);
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
                std::swap(lowerUpper[index(col, col2, N)], lowerUpper[index(pivotRow, col2, N)]);
            }
            std::swap(pivot[col], pivot[pivotRow]);
            ++swapCount;
        }

        // eliminate
        const double pivotValue{lowerUpper[index(col, col, N)]};
        for (std::size_t row = col + 1; row < N; ++row) {
            lowerUpper[index(row, col, N)] /= pivotValue; // L (i,k)
            const double multiplier{lowerUpper[index(row, col, N)]};
            for (std::size_t col2 = col + 1; col2 < N; ++col2) {
                lowerUpper[index(row, col2, N)] -= multiplier * lowerUpper[index(col, col2, N)];
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
void lowerUpperSolve(const double* LU, const std::size_t* pivot, const double* b, double* x, std::size_t N) {
    // Apply permutation: x = Pb
    // store y in x temporarily
    for (std::size_t row = 0; row < N; ++row) {
        const std::size_t permRow{pivot[row]};
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

double SlaterPlaneWave::log_abs_det(const Particles& particles) {
    const std::size_t N{num_orbitals_get()};

    const double* pos_x{particles.pos_x_get()};
    const double* pos_y{particles.pos_y_get()};
    const double* pos_z{particles.pos_z_get()};

    double* RESTRICT det_matrix{determinant_get()};
    double* RESTRICT lower_upper_matrix{lower_upper_get()};
    double* RESTRICT inv_det_matrix{inv_determinant_get()};
    std::size_t* RESTRICT pivot_vector{pivot_get()};

    const double* k_x_comp{k_vector_x_get()};
    const double* k_y_comp{k_vector_y_get()};
    const double* k_z_comp{k_vector_z_get()};

    // Build determinant matrix D
    for (std::size_t particle = 0; particle < N; ++particle) {
        const double P_X{pos_x[particle]};
        const double P_Y{pos_y[particle]};
        const double P_Z{pos_z[particle]};

        for (std::size_t orbital = 0; orbital < N; ++orbital) {
            const double K_DOT_R = k_x_comp[orbital] * P_X + k_y_comp[orbital] * P_Y + k_z_comp[orbital] * P_Z;

            det_matrix[index(particle, orbital, N)] = std::cos(K_DOT_R);
        }
    }

    // Copy determinant matrix into LU storage
    std::copy_n(det_matrix, N * N, lower_upper_matrix);

    // Perform LU decomposition
    (void)lowerUpperDecompose(lower_upper_matrix, pivot_vector, N);

    // Compute log|det(D)| = Σ log|U_ii|
    double logAbsoluteDeterminant = 0.0;

    for (std::size_t diag = 0; diag < N; ++diag) {
        const double Uii = lower_upper_matrix[index(diag, diag, N)];

        const double absUii = std::abs(Uii);

        if (absUii == 0.0 || !std::isfinite(absUii)) {
            return -std::numeric_limits<double>::infinity();
        }

        logAbsoluteDeterminant += std::log(absUii);
    }

    std::vector<double> rhs(N);
    std::vector<double> solution(N);

    for (std::size_t column = 0; column < N; ++column) {
        std::fill(rhs.begin(), rhs.end(), 0.0);
        rhs[column] = 1.0;

        lowerUpperSolve(lower_upper_matrix, pivot_vector, rhs.data(), solution.data(), N);

        for (std::size_t row = 0; row < N; ++row) {
            inv_det_matrix[index(row, column, N)] = solution[row];
        }
    }

    return logAbsoluteDeterminant;
}

void SlaterPlaneWave::add_derivatives(const Particles& particles, double* RESTRICT grad_x, double* RESTRICT grad_y,
                                      double* RESTRICT grad_z, double* RESTRICT laplacian) const noexcept {
    const std::size_t N{num_orbitals_get()};

    const double* RESTRICT pos_x{particles.pos_x_get()};
    const double* RESTRICT pos_y{particles.pos_y_get()};
    const double* RESTRICT pos_z{particles.pos_z_get()};

    const double* RESTRICT K_X{k_vector_x_get()};
    const double* RESTRICT K_Y{k_vector_y_get()};
    const double* RESTRICT K_Z{k_vector_z_get()};

    const double* RESTRICT INV_DET{inv_determinant_get()};

    for (std::size_t particle = 0; particle < N; ++particle) {
        const double P_X{pos_x[particle]};
        const double P_Y{pos_y[particle]};
        const double P_Z{pos_z[particle]};

        double d_Log_det_dx{}, d_Log_det_dy{}, d_Log_det_dz{};

        // Σ_j (D^{-1})_{j,particle} * (∇^2 D_{particle,j})
        double laplace_det_term{};

        for (std::size_t orbital = 0; orbital < N; ++orbital) {
            const double K_X_ORBITAL{K_X[orbital]};
            const double K_Y_ORBITAL{K_Y[orbital]};
            const double K_Z_ORBITAL{K_Z[orbital]};

            const double K_DOT_R{K_X_ORBITAL * P_X + K_Y_ORBITAL * P_Y + K_Z_ORBITAL * P_Z};
            const double COS_TERM{std::cos(K_DOT_R)};
            const double SIN_TERM{std::sin(K_DOT_R)};
            // D(particle,orbital) = cos(k·r)
            // ∇_particle D = -sin(k·r) * k
            const double dD_dx{-SIN_TERM * K_X_ORBITAL};
            const double dD_dy{-SIN_TERM * K_Y_ORBITAL};
            const double dD_dz{-SIN_TERM * K_Z_ORBITAL};

            // ∇^2 D = -cos(k·r) * |k|^2
            const double K_SQ{K_X_ORBITAL * K_X_ORBITAL + K_Y_ORBITAL * K_Y_ORBITAL + K_Z_ORBITAL * K_Z_ORBITAL};
            const double LAP_D{-COS_TERM * K_SQ};
            // WEIGHT = (D^{-1})_{orbital,particle}
            // invD_ is row-major, so entry (row=orbital, col=particle):
            const double WEIGHT{INV_DET[index(orbital, particle, N)]};

            d_Log_det_dx += WEIGHT * dD_dx;
            d_Log_det_dy += WEIGHT * dD_dy;
            d_Log_det_dz += WEIGHT * dD_dz;

            laplace_det_term += WEIGHT * LAP_D;
        }
        // ∇^2 log det D = Σ_j (D^{-1})_{j,i} ∇^2 D_{i,j}  -  ||∇ log det D||^2
        const double GRAD_SQ{d_Log_det_dx * d_Log_det_dx + d_Log_det_dy * d_Log_det_dy + d_Log_det_dz * d_Log_det_dz};
        const double LAP_LOG_DET{laplace_det_term - GRAD_SQ};

        grad_x[particle] += d_Log_det_dx;
        grad_y[particle] += d_Log_det_dy;
        grad_z[particle] += d_Log_det_dz;
        laplacian[particle] += LAP_LOG_DET;
    }
}