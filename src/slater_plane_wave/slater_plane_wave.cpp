#include "slater_plane_wave.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

SlaterPlaneWave::SlaterPlaneWave(std::size_t N, double L)
    : num_orbitals_{N}, matrix_size_{N * N}, box_length_{L}, int_vectors_{N, NUM_INT_VECTORS_},
      k_vector_{N, NUM_K_VECTORS_}, matrices_{N * N, NUM_MATRIX_} {

    const std::size_t N{num_orbitals_get()};

    int* RESTRICT n_x{n_vector_x_get()};
    int* RESTRICT n_y{n_vector_y_get()};
    int* RESTRICT n_z{n_vector_z_get()};

    // Start at 1 since n[0] = (0,0,0):
    // TODO: implement the algorithm to determine n
    for (std::size_t i = 1; i < N; ++i) {
        n_x[i] = 0;
        n_y[i] = 0;
        n_z[i] = 0;
    }

    double* RESTRICT k_x{k_vector_x_get()};
    double* RESTRICT k_y{k_vector_y_get()};
    double* RESTRICT k_z{k_vector_z_get()};

    const double inv_L{1 / box_length_get()};

    // Follows the calculation: K = (2pi/L) * n;
    for (std::size_t i = 0; i < N; ++i) {
        k_x[i] = 2 * std::numbers::pi * inv_L * static_cast<double>(n_x[i]);
        k_y[i] = 2 * std::numbers::pi * inv_L * static_cast<double>(n_y[i]);
        k_z[i] = 2 * std::numbers::pi * inv_L * static_cast<double>(n_z[i]);
    }
};

double SlaterPlaneWave::log_abs_det(const Particles& particles) {
    const std::size_t N{num_orbitals_get()};

    const double* pos_x{particles.pos_x_get()};
    const double* pos_y{particles.pos_y_get()};
    const double* pos_z{particles.pos_z_get()};

    double* RESTRICT det_matrix{determinant_get()};
    double* RESTRICT lower_upper_matrix{lower_upper_get()};
    double* RESTRICT inv_det_matrix{inv_determinant_get()};
    int* RESTRICT pivot_vector{pivot_get()};

    const double* k_x_comp{k_vector_x_get()};
    const double* k_y_comp{k_vector_y_get()};
    const double* k_z_comp{k_vector_z_get()};

    // Build determinant matrix D
    for (std::size_t particle = 0; particle < N; ++particle) {
        const double p_x{pos_x[particle]};
        const double p_y{pos_y[particle]};
        const double p_z{pos_z[particle]};

        for (std::size_t orbital = 0; orbital < N; ++orbital) {
            const double k_dot_r = k_x_comp[orbital] * p_x + k_y_comp[orbital] * p_y + k_z_comp[orbital] * p_z;

            det_matrix[index(particle, orbital, N)] = std::cos(k_dot_r);
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

    const double* RESTRICT k_x{k_vector_x_get()};
    const double* RESTRICT k_y{k_vector_y_get()};
    const double* RESTRICT k_z{k_vector_z_get()};

    const double* RESTRICT inv_det{inv_determinant_get()};

    for (std::size_t particle = 0; particle < N; ++particle) {
        const double p_x{pos_x[particle]};
        const double p_y{pos_y[particle]};
        const double p_z{pos_z[particle]};

        double d_Log_det_dx{}, d_Log_det_dy{}, d_Log_det_dz{};

        // Σ_j (D^{-1})_{j,particle} * (∇^2 D_{particle,j})
        double laplace_det_term{};

        for (std::size_t orbital = 0; orbital < N; ++orbital) {
            const double k_x_orbital{k_x[orbital]};
            const double k_y_orbital{k_y[orbital]};
            const double k_z_orbital{k_z[orbital]};

            const double k_dot_R{k_x_orbital * p_x + k_y_orbital * p_y + k_z_orbital * p_z};
            const double cos_term{std::cos(k_dot_R)};
            const double sin_term{std::sin(k_dot_R)};
            // D(particle,orbital) = cos(k·r)
            // ∇_particle D = -sin(k·r) * k
            const double dD_dx{-sin_term * k_x_orbital};
            const double dD_dy{-sin_term * k_y_orbital};
            const double dD_dz{-sin_term * k_z_orbital};

            // ∇^2 D = -cos(k·r) * |k|^2
            const double k_sq{k_x_orbital * k_x_orbital + k_y_orbital * k_y_orbital + k_z_orbital * k_z_orbital};
            const double lap_D{-cos_term * k_sq};
            // weight = (D^{-1})_{orbital,particle}
            // invD_ is row-major, so entry (row=orbital, col=particle):
            const double weight{inv_det[index(orbital, particle, N)]};

            d_Log_det_dx += weight * dD_dx;
            d_Log_det_dy += weight * dD_dy;
            d_Log_det_dz += weight * dD_dz;

            laplace_det_term += weight * lap_D;
        }
        // ∇^2 log det D = Σ_j (D^{-1})_{j,i} ∇^2 D_{i,j}  -  ||∇ log det D||^2
        const double grad_sq{d_Log_det_dx * d_Log_det_dx + d_Log_det_dy * d_Log_det_dy + d_Log_det_dz * d_Log_det_dz};
        const double lap_log_det{laplace_det_term - grad_sq};

        grad_x[particle] += d_Log_det_dx;
        grad_y[particle] += d_Log_det_dy;
        grad_z[particle] += d_Log_det_dz;

        laplacian[particle] += lap_log_det;
    }
}

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
int lowerUpperDecompose(double* lowerUpper, int* pivot, std::size_t N) {
    // Track row swaps
    int swapCount{};

    for (std::size_t row = 0; row < N; ++row)
        pivot[row] = static_cast<int>(row);

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
void lowerUpperSolve(const double* LU, const int* pivot, const double* b, double* x, std::size_t N) {
    // Apply permutation: x = Pb
    // store y in x temporarily
    for (std::size_t row = 0; row < N; ++row) {
        const std::size_t permRow{static_cast<std::size_t>(pivot[row])};
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