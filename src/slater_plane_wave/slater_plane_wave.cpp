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
int lower_upper_decomp(double* lowerUpper, int* pivot, std::size_t N) {
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
void solve_lower_upper(const double* LU, const int* pivot, const double* b, double* x, std::size_t N) {
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

// Canonical representative rule for +-n deduplication:
// The canonical form is: the first nonzero component is positive
// The zero vector (0,0,0) is its own canonical representative
bool is_canonical(int n_x, int n_y, int n_z) {
    if (n_x > 0)
        return true;
    if (n_x < 0)
        return false;
    if (n_y > 0)
        return true;
    if (n_y < 0)
        return false;
    if (n_z > 0)
        return true;
    if (n_z < 0)
        return false;

    // (0,0,0)
    return true;
}

} // namespace

// Real plane-wave basis from complex exponentials
// The complex basis e^{ik·r} = cos(k·r) + i·sin(k·r) gives two
// linearly independent real orbitals per k-vector:
//   φ_cos(r) = cos(k·r)
//   φ_sin(r) = sin(k·r)
//
// Since cos(k·r) = cos(-k·r) and sin(k·r) = -sin(-k·r),
// the ±k pair maps to the same {cos, sin} pair. We deduplicate
// by keeping only the canonical +n representative (first nonzero
// component positive).
//
// k = 0 is special: sin(0) = 0 for all r, so it contributes
// only one orbital (cos(0) = 1).
//
// Orbital count: N = 1 + 2 × (number of nonzero canonical k-vectors)
// Closed shells: N = 1, 7, 19, 27, 33, 57, ...
SlaterPlaneWave::SlaterPlaneWave(std::size_t num_particles, double box_lengthL)
    : num_orbitals_{num_particles}, matrix_size_{num_particles * num_particles}, box_length_{box_lengthL},
      orbital_k_index_(num_particles), orbital_type_(num_particles, 0), int_vectors_{num_particles, NUM_INT_VECTORS_},
      double_vectors_{num_particles, NUM_DOUBLE_VECTORS_}, matrices_{num_particles * num_particles, NUM_MATRIX_} {

    const std::size_t N{num_orbitals_get()};

    // Generate deduplicated canonical n-vectors
    // Only keep canonical representatives: first nonzero component positive
    // k=0 produces 1 orbital (cos = 1)
    // Each nonzero k produces 2 orbitals (cos, sin)

    struct nVectorCandidate {
        int n_cand_x;
        int n_cand_y;
        int n_cand_z;
        int n_mag_sq;
    };

    // Increase vector size to be safe
    const int N_MAX{static_cast<int>(N + SIMD_BYTES)};
    const std::size_t side{static_cast<std::size_t>((2 * N_MAX + 1))};

    // Vector for all possible n-vector candidates
    std::vector<nVectorCandidate> n_candidates{};
    n_candidates.reserve(side * side * side);

    // Fill candidate vector with all possible states from [-N, N]
    for (int new_x = -N_MAX; new_x <= N_MAX; ++new_x) {
        for (int new_y = -N_MAX; new_y <= N_MAX; ++new_y) {
            for (int new_z = -N_MAX; new_z <= N_MAX; ++new_z) {
                if (!is_canonical(new_x, new_y, new_z))
                    continue;

                const int new_mag_sq{new_x * new_x + new_y * new_y + new_z * new_z};
                n_candidates.emplace_back(new_x, new_y, new_z, new_mag_sq);
            }
        }
    }

    // Sort n-vector states to go from smallest magnitude to largest
    std::sort(n_candidates.begin(), n_candidates.end(), [](const nVectorCandidate& a, const nVectorCandidate& b) {
        return std::tie(a.n_mag_sq, a.n_cand_x, a.n_cand_y, a.n_cand_z) < 
           std::tie(b.n_mag_sq, b.n_cand_x, b.n_cand_y, b.n_cand_z);
    });

    // Assign orbitals
    // Orbital 0: k=0 -> cos(0 dot r) = cos(0) = 1
    // For each nonzero canonical k: orbital 2m-1 -> cos(k dot r), orbital 2m -> sin(k dot r)

    int* RESTRICT n_x{n_vector_x_get()};
    int* RESTRICT n_y{n_vector_y_get()};
    int* RESTRICT n_z{n_vector_z_get()};

    auto& orb_k_idx{orbital_k_index_get()};
    auto& orb_type{orbital_type_get()};

    std::size_t orb_idx{};
    std::size_t k_idx{};

    for (std::size_t c = 0; c < n_candidates.size() && orb_idx < N; ++c) {
        const auto& cand{n_candidates[c]};
        const bool mag_not_zero{cand.n_mag_sq != 0};

        n_x[k_idx] = cand.n_cand_x;
        n_y[k_idx] = cand.n_cand_y;
        n_z[k_idx] = cand.n_cand_z;

        // Cos orbital (every k-vector gets one)
        orb_k_idx[orb_idx] = k_idx;
        orb_type[orb_idx] = 0;
        ++orb_idx;

        // Sin orbital (only for non-zero k - sin(0) = 0 is singular)
        if (mag_not_zero && orb_idx < N) {
            orb_k_idx[orb_idx] = k_idx;
            orb_type[orb_idx] = 1;
            ++orb_idx;
        }

        ++k_idx;
    }

    num_unique_k_set() = k_idx;

    double* RESTRICT k_x{k_vector_x_get()};
    double* RESTRICT k_y{k_vector_y_get()};
    double* RESTRICT k_z{k_vector_z_get()};

    const double inv_L{1.0 / box_length_get()};

    // Follows the calculation: K = (2pi/L) * n;
    for (std::size_t i = 0; i < k_idx; ++i) {
        k_x[i] = 2 * std::numbers::pi * inv_L * static_cast<double>(n_x[i]);
        k_y[i] = 2 * std::numbers::pi * inv_L * static_cast<double>(n_y[i]);
        k_z[i] = 2 * std::numbers::pi * inv_L * static_cast<double>(n_z[i]);
    }
};

// Slater matrix D_{i,j} = φ_j(r_i)
// Each orbital j has an associated k-vector (via orbital_k_index)
// and a type (via orbital_type): 0 = cos, 1 = sin.
//
//   type 0: D_{i,j} = cos(k_j · r_i)
//   type 1: D_{i,j} = sin(k_j · r_i)
//
// log|Ψ_Slater| = log|det(D)|
double SlaterPlaneWave::log_abs_det(const Particles& particles) {
    const std::size_t N{num_orbitals_get()};
    const std::size_t padded_N{particles.padding_stride_get()};

    const double* RESTRICT pos_x{particles.pos_x_get()};
    const double* RESTRICT pos_y{particles.pos_y_get()};
    const double* RESTRICT pos_z{particles.pos_z_get()};

    const double* RESTRICT k_x_comp{k_vector_x_get()};
    const double* RESTRICT k_y_comp{k_vector_y_get()};
    const double* RESTRICT k_z_comp{k_vector_z_get()};

    const auto& k_index{orbital_k_index_get()};
    const auto& orb_type{orbital_type_get()};

    double* RESTRICT det_matrix{determinant_get()};
    double* RESTRICT lower_upper_matrix{lower_upper_get()};
    double* RESTRICT inv_det_matrix{inv_determinant_get()};

    int* RESTRICT pivot_vector{pivot_get()};

    double* RESTRICT rhs{rhs_get()};
    double* RESTRICT solution{solution_get()};

    // Build determinant matrix D
    for (std::size_t particle = 0; particle < N; ++particle) {
        const double p_x{pos_x[particle]};
        const double p_y{pos_y[particle]};
        const double p_z{pos_z[particle]};

        for (std::size_t orbital = 0; orbital < N; ++orbital) {
            const std::size_t k_idx{k_index[orbital]};
            const double k_dot_r = k_x_comp[k_idx] * p_x + k_y_comp[k_idx] * p_y + k_z_comp[k_idx] * p_z;

            const double type{static_cast<double>(orb_type[orbital])};

            double cos_term{std::cos(k_dot_r)};
            double sin_term{std::sin(k_dot_r)};

            det_matrix[index(particle, orbital, N)] = cos_term + type * (sin_term - cos_term);
        }
    }

    // Copy determinant matrix into LU storage
    std::copy_n(det_matrix, N * N, lower_upper_matrix);

    // Perform LU decomposition
    (void)lower_upper_decomp(lower_upper_matrix, pivot_vector, N);

    // Compute log|det(D)| = Σ log|U_ii|
    double log_abs_det = 0.0;

    for (std::size_t diag = 0; diag < N; ++diag) {
        const double U_ii = lower_upper_matrix[index(diag, diag, N)];

        const double abs_U_ii = std::abs(U_ii);

        if (abs_U_ii == 0.0 || !std::isfinite(abs_U_ii)) {
            return -std::numeric_limits<double>::infinity();
        }

        log_abs_det += std::log(abs_U_ii);
    }

    for (std::size_t column = 0; column < N; ++column) {
        // Using pointer arithmetic to start at rhs[start]
        // and jump to rhs[end] - used padded_N since
        // SIMD alignment might mess up pointer math
        std::fill(rhs, rhs + padded_N, 0.0);
        rhs[column] = 1.0;

        solve_lower_upper(lower_upper_matrix, pivot_vector, rhs, solution, N);

        for (std::size_t row = 0; row < N; ++row) {
            inv_det_matrix[index(row, column, N)] = solution[row];
        }
    }

    return log_abs_det;
}

double* SlaterPlaneWave::build_row(std::size_t particle, const Particles& particles) noexcept {
    const std::size_t N{num_orbitals_get()};

    const double p_x{particles.pos_x_get()[particle]};
    const double p_y{particles.pos_y_get()[particle]};
    const double p_z{particles.pos_z_get()[particle]};

    const double* RESTRICT k_x{k_vector_x_get()};
    const double* RESTRICT k_y{k_vector_y_get()};
    const double* RESTRICT k_z{k_vector_z_get()};

    const auto& k_index{orbital_k_index_get()};
    const auto& orb_type{orbital_type_get()};

    double* RESTRICT row{new_row_get()};

    for (std::size_t orbital = 0; orbital < N; ++orbital) {
        const std::size_t k_idx{k_index[orbital]};
        const double k_dot_r{k_x[k_idx] * p_x + k_y[k_idx] * p_y + k_z[k_idx] * p_z};

        const double type{static_cast<double>(orb_type[orbital])};

        double cos_term{std::cos(k_dot_r)};
        double sin_term{std::sin(k_dot_r)};

        row[orbital] = cos_term + type * (sin_term - cos_term);
    }

    return row;
}

double SlaterPlaneWave::determinant_ratio(std::size_t particle, const double* new_row) const noexcept {
    const std::size_t N{num_orbitals_get()};
    const double* RESTRICT inv_det{inv_determinant_get()};

    double ratio{};
    for (std::size_t j = 0; j < N; ++j) {
        ratio += new_row[j] * inv_det[index(j, particle, N)];
    }

    return ratio;
}

void SlaterPlaneWave::accept_move(std::size_t particle, const double* new_row, double ratio) noexcept {
    const std::size_t N{num_orbitals_get()};

    double* RESTRICT inv_det{inv_determinant_get()};
    double* RESTRICT det_matrix{determinant_get()};
    double* RESTRICT inv_d_col{inv_d_col_get()};

    const double inv_ratio{1.0 / ratio};

    // Cache particle row column j for inv_D before changing
    for (std::size_t j = 0; j < N; ++j) {
        inv_d_col[j] = inv_det[index(j, particle, N)];
    }

    // Follows Sherman-Morrison update:
    for (std::size_t k = 0; k < N; ++k) {
        if (k == particle) {
            // Special case: column p just scales by 1/R
            for (std::size_t j = 0; j < N; ++j) {
                inv_det[index(j, k, N)] = inv_d_col[j] * inv_ratio;
            }
        } else {
            // Compute s_k = Σ_m new_row[m] * (D⁻¹_old)[m, k]
            double s_k{};
            for (std::size_t m = 0; m < N; ++m) {
                s_k += new_row[m] * inv_det[index(m, k, N)];
            }

            const double factor{s_k * inv_ratio};
            for (std::size_t j = 0; j < N; ++j) {
                inv_det[index(j, k, N)] -= inv_d_col[j] * factor;
            }
        }
    }

    // Patch row `particle` of D to match the new positions:
    for (std::size_t j = 0; j < N; ++j) {
        det_matrix[index(particle, j, N)] = new_row[j];
    }
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

    const auto& k_index{orbital_k_index_get()};
    const auto& o_type{orbital_type_get()};

    const double* RESTRICT inv_det{inv_determinant_get()};

    std::vector<double> k_sq_cache;
    std::vector<double> sin_cache;
    std::vector<double> cos_cache;
    std::size_t num_k_vectors{num_unique_k_get()};
    k_sq_cache.resize(num_k_vectors);
    sin_cache.resize(num_k_vectors);
    cos_cache.resize(num_k_vectors);

    for (std::size_t k = 0; k < num_k_vectors; ++k) {
        k_sq_cache[k] = k_x[k] * k_x[k] + k_y[k] * k_y[k] + k_z[k] * k_z[k];
    }

    for (std::size_t particle = 0; particle < N; ++particle) {
        const double p_x{pos_x[particle]};
        const double p_y{pos_y[particle]};
        const double p_z{pos_z[particle]};

        double d_log_det_dx{}, d_log_det_dy{}, d_log_det_dz{};

        // Σ_j (D^{-1})_{j,particle} * (∇^2 D_{particle,j})
        double laplace_det_term{};

        for (std::size_t orbital = 0; orbital < N; ++orbital) {
            const std::size_t k_idx{k_index[orbital]};

            const double k_x_orbital{k_x[k_idx]};
            const double k_y_orbital{k_y[k_idx]};
            const double k_z_orbital{k_z[k_idx]};

            const double k_dot_R{k_x_orbital * p_x + k_y_orbital * p_y + k_z_orbital * p_z};
            const double k_sq{k_x_orbital * k_x_orbital + k_y_orbital * k_y_orbital + k_z_orbital * k_z_orbital};

            // ∇_i log|det D| = Σ_j (D⁻¹)_{j,i} · ∇_i D_{i,j}
            //
            // ∇²_i log|det D| = Σ_j (D⁻¹)_{j,i} · ∇²_i D_{i,j}
            //                    - ||∇_i log|det D|||²
            //
            // The cross term arises from d/dx (f'/f) = f''/f - (f'/f)²
            // applied to each component of the gradient.

            double dD_dx{}, dD_dy{}, dD_dz{};
            double lap_D{};

            // Replaces the old if else branch with branchless math
            const double type{static_cast<double>(o_type[orbital])};

            double cos_term{std::cos(k_dot_R)};
            double sin_term{std::sin(k_dot_R)};

            const double grad_factor{-sin_term + type * (sin_term + cos_term)};
            const double lap_factor{-cos_term + type * (cos_term - sin_term)};

            // D =       (o_type[orbital] == 0) ? cos(k dot r)          : sin(k dot r)
            // grad(D) = (o_type[orbital] == 0) ? -sin(k dot r) * k     : cos(k dot r) * k
            // lap (D) = (o_type[orbital] == 0) ? -cos(k dot r) * |k|^2 : -sin(k dor r) * |k|^2
            dD_dx = k_x_orbital * grad_factor;
            dD_dy = k_y_orbital * grad_factor;
            dD_dz = k_z_orbital * grad_factor;
            lap_D = k_sq * lap_factor;

            const double weight{inv_det[index(orbital, particle, N)]};

            d_log_det_dx += weight * dD_dx;
            d_log_det_dy += weight * dD_dy;
            d_log_det_dz += weight * dD_dz;

            laplace_det_term += weight * lap_D;
        }
        const double grad_sq{d_log_det_dx * d_log_det_dx + d_log_det_dy * d_log_det_dy + d_log_det_dz * d_log_det_dz};
        const double lap_log_det{laplace_det_term - grad_sq};

        grad_x[particle] += d_log_det_dx;
        grad_y[particle] += d_log_det_dy;
        grad_z[particle] += d_log_det_dz;

        laplacian[particle] += lap_log_det;
    }
}