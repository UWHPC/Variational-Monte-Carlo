#include "slater_plane_wave.hpp"
#include "../utilities/matrix.hpp"
#include "particles/particles.hpp"
#include "utilities/aligned_soa.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <omp.h>
#include <vector>

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
SlaterPlaneWave::SlaterPlaneWave(const Particles& particles, double box_lengthL)
    : num_orbitals_{particles.num_particles_get()},
      matrix_size_{particles.num_particles_get() * particles.num_particles_get()},
      box_length_{box_lengthL}, orbital_k_index_(particles.num_particles_get()),
      orbital_type_(particles.num_particles_get(), 0),
      int_vectors_{particles.num_particles_get(), NUM_INT_VECTORS_},
      double_vectors_{particles.num_particles_get(), NUM_DOUBLE_VECTORS_}, trig_cache_{},
      matrices_{particles.num_particles_get() * particles.num_particles_get(), NUM_MATRIX_} {

    const std::size_t N{num_orbitals_get()};
    const std::size_t num_particles{particles.num_particles_get()};

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
    const int N_MAX{static_cast<int>(std::ceil(std::cbrt(static_cast<double>(N)))) + 2};
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
    std::sort(n_candidates.begin(), n_candidates.end(),
              [](const nVectorCandidate& a, const nVectorCandidate& b) {
                  return std::tie(a.n_mag_sq, a.n_cand_x, a.n_cand_y, a.n_cand_z) <
                         std::tie(b.n_mag_sq, b.n_cand_x, b.n_cand_y, b.n_cand_z);
              });

    // Assign orbitals
    // Orbital 0: k=0 -> cos(0 dot r) = cos(0) = 1
    // For each nonzero canonical k: orbital 2m-1 -> cos(k dot r), orbital 2m -> sin(k dot r)

    int* RESTRICT n_x{n_vector_x_get()};
    int* RESTRICT n_y{n_vector_y_get()};
    int* RESTRICT n_z{n_vector_z_get()};

    ASSUME_ALIGNED(n_x, SIMD_BYTES);
    ASSUME_ALIGNED(n_y, SIMD_BYTES);
    ASSUME_ALIGNED(n_z, SIMD_BYTES);


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

    ASSUME_ALIGNED(k_x, SIMD_BYTES);
    ASSUME_ALIGNED(k_y, SIMD_BYTES);
    ASSUME_ALIGNED(k_z, SIMD_BYTES);

    const double inv_L{1.0 / box_length_get()};

// Follows the calculation: K = (2pi/L) * n;
#pragma omp simd
    for (std::size_t i = 0; i < k_idx; ++i) {
        k_x[i] = 2 * std::numbers::pi * inv_L * static_cast<double>(n_x[i]);
        k_y[i] = 2 * std::numbers::pi * inv_L * static_cast<double>(n_y[i]);
        k_z[i] = 2 * std::numbers::pi * inv_L * static_cast<double>(n_z[i]);
    }

    trig_cache_ = AlignedSoA<double>(num_particles * num_unique_k_get(), NUM_TRIG_ARRAYS_);
    trig_scratch_ = AlignedSoA<double>(num_unique_k_get(), NUM_SCRATCH_TRIG_);

    std::size_t trig_size{trig_cache_.num_elements()};
    std::fill_n(sin_cache_get(), trig_size, 0.0);
    std::fill_n(cos_cache_get(), trig_size, 0.0);
};

void SlaterPlaneWave::save_trig_row(std::size_t particle) noexcept {
    const std::size_t num_k{num_unique_k_get()};
    const std::size_t offset{particle * num_k};

    std::memcpy(trig_scratch_[SIN_SAVED_], sin_cache_get() + offset, num_k * sizeof(double));
    std::memcpy(trig_scratch_[COS_SAVED_], cos_cache_get() + offset, num_k * sizeof(double));
}

void SlaterPlaneWave::restore_trig_row(std::size_t particle) noexcept {
    const std::size_t num_k{num_unique_k_get()};
    const std::size_t offset{particle * num_k};

    std::memcpy(sin_cache_get() + offset, trig_scratch_[SIN_SAVED_], num_k * sizeof(double));
    std::memcpy(cos_cache_get() + offset, trig_scratch_[COS_SAVED_], num_k * sizeof(double));
}

void SlaterPlaneWave::update_trig_cache(std::size_t particle, const Particles& particles) noexcept {
    const std::size_t num_k{num_unique_k_get()};

    const double px{particles.pos_x_get()[particle]};
    const double py{particles.pos_y_get()[particle]};
    const double pz{particles.pos_z_get()[particle]};

    const double* RESTRICT kx{k_vector_x_get()};
    const double* RESTRICT ky{k_vector_y_get()};
    const double* RESTRICT kz{k_vector_z_get()};

    double* RESTRICT c_row{cos_cache_get() + particle * num_k};
    double* RESTRICT s_row{sin_cache_get() + particle * num_k};

    ASSUME_ALIGNED(kx, SIMD_BYTES);
    ASSUME_ALIGNED(ky, SIMD_BYTES);
    ASSUME_ALIGNED(kz, SIMD_BYTES);

    ASSUME_ALIGNED(c_row, SIMD_BYTES);
    ASSUME_ALIGNED(s_row, SIMD_BYTES);

#pragma omp simd
    for (std::size_t k = 0; k < num_k; ++k) {
        const double dot{kx[k] * px + ky[k] * py + kz[k] * pz};

        PORTABLE_SINCOS(dot, &s_row[k], &c_row[k]);
    }
}

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

    const std::size_t num_k{num_unique_k_get()};
    double* RESTRICT cos_cache{cos_cache_get()};
    double* RESTRICT sin_cache{sin_cache_get()};

    ASSUME_ALIGNED(pos_x, SIMD_BYTES);
    ASSUME_ALIGNED(pos_y, SIMD_BYTES);
    ASSUME_ALIGNED(pos_z, SIMD_BYTES);

    ASSUME_ALIGNED(k_x_comp, SIMD_BYTES);
    ASSUME_ALIGNED(k_y_comp, SIMD_BYTES);
    ASSUME_ALIGNED(k_z_comp, SIMD_BYTES);

    ASSUME_ALIGNED(det_matrix, SIMD_BYTES);
    ASSUME_ALIGNED(lower_upper_matrix, SIMD_BYTES);
    ASSUME_ALIGNED(inv_det_matrix, SIMD_BYTES);

    ASSUME_ALIGNED(pivot_vector, SIMD_BYTES);

    ASSUME_ALIGNED(rhs, SIMD_BYTES);
    ASSUME_ALIGNED(solution, SIMD_BYTES);

    ASSUME_ALIGNED(cos_cache, SIMD_BYTES);
    ASSUME_ALIGNED(sin_cache, SIMD_BYTES);


    // Build determinant matrix D
    for (std::size_t particle = 0; particle < N; ++particle) {
        const double p_x{pos_x[particle]};
        const double p_y{pos_y[particle]};
        const double p_z{pos_z[particle]};

        const std::size_t offset{particle * num_k};
        
#pragma omp simd
        for (std::size_t k = 0; k < num_k; ++k) {
            const double dot{k_x_comp[k] * p_x + k_y_comp[k] * p_y + k_z_comp[k] * p_z};
            const std::size_t i{offset + k};

            PORTABLE_SINCOS(dot, &sin_cache[i], &cos_cache[i]);
        }

// Build D from cache
#pragma omp simd
        for (std::size_t orbital = 0; orbital < N; ++orbital) {
            const std::size_t k_idx{k_index[orbital]};

            const double type{static_cast<double>(orb_type[orbital])};
            const double cos_term{cos_cache[offset + k_idx]};
            const double sin_term{sin_cache[offset + k_idx]};

            det_matrix[particle * N + orbital] = cos_term + type * (sin_term - cos_term);
        }
    }

    // Copy determinant matrix into LU storage
    std::copy_n(det_matrix, N * N, lower_upper_matrix);

    // Perform LU decomposition
    (void)lower_upper_decomp(lower_upper_matrix, pivot_vector, N);

    // Compute log|det(D)| = Σ log|U_ii|
    double log_abs_det{};

#pragma omp simd reduction(+ : log_abs_det)
    for (std::size_t diag = 0; diag < N; ++diag) {
        const double U_ii{lower_upper_matrix[diag * N + diag]};
        const double abs_U_ii{std::abs(U_ii)};

        log_abs_det += std::log(abs_U_ii);
    }

    if (!std::isfinite(log_abs_det)) {
        return -std::numeric_limits<double>::infinity();
    }

    for (std::size_t column = 0; column < N; ++column) {
        // Using pointer arithmetic to start at rhs[start]
        // and jump to rhs[end] - used padded_N since
        // SIMD alignment might mess up pointer math
        std::fill(rhs, rhs + padded_N, 0.0);
        rhs[column] = 1.0;

        solve_lower_upper(lower_upper_matrix, pivot_vector, rhs, solution, N);

        for (std::size_t row = 0; row < N; ++row) {
            inv_det_matrix[column * N + row] = solution[row];
        }
    }

    return log_abs_det;
}

double* SlaterPlaneWave::build_row(std::size_t particle) noexcept {
    const std::size_t N{num_orbitals_get()};
    const std::size_t num_K{num_unique_k_get()};

    const auto& k_index{orbital_k_index_get()};
    const auto& orb_type{orbital_type_get()};

    double* RESTRICT row{new_row_get()};
    double* RESTRICT sin_cache{sin_cache_get()};
    double* RESTRICT cos_cache{cos_cache_get()};

    ASSUME_ALIGNED(row, SIMD_BYTES);
    ASSUME_ALIGNED(sin_cache, SIMD_BYTES);
    ASSUME_ALIGNED(cos_cache, SIMD_BYTES);


#pragma omp simd
    for (std::size_t orbital = 0; orbital < N; ++orbital) {
        const std::size_t k_idx{k_index[orbital]};

        const double type{static_cast<double>(orb_type[orbital])};

        const double cos_term{cos_cache[particle * num_K + k_idx]};
        const double sin_term{sin_cache[particle * num_K + k_idx]};

        row[orbital] = cos_term + type * (sin_term - cos_term);
    }

    return row;
}

double SlaterPlaneWave::determinant_ratio(std::size_t particle,
                                          const double* new_row) const noexcept {
    const std::size_t N{num_orbitals_get()};
    const double* RESTRICT inv_det{inv_determinant_get()};
    ASSUME_ALIGNED(inv_det, SIMD_BYTES);

    double ratio{};
#pragma omp simd reduction(+ : ratio)
    for (std::size_t j = 0; j < N; ++j) {
        ratio += new_row[j] * inv_det[particle * N + j];
    }

    return ratio;
}

void SlaterPlaneWave::accept_move(std::size_t particle, const double* new_row,
                                  double ratio) noexcept {
    const std::size_t N{num_orbitals_get()};

    double* RESTRICT inv_det{inv_determinant_get()};
    double* RESTRICT det_matrix{determinant_get()};
    double* RESTRICT inv_d_col{inv_d_col_get()};

    ASSUME_ALIGNED(inv_det, SIMD_BYTES);
    ASSUME_ALIGNED(det_matrix, SIMD_BYTES);
    ASSUME_ALIGNED(inv_d_col, SIMD_BYTES);

    const double inv_ratio{1.0 / ratio};
    const std::size_t p_offset{particle * N}; // Pre-calculate particle row offset

// Cache particle row column j for inv_D before changing
#pragma omp simd
    for (std::size_t j = 0; j < N; ++j) {
        inv_d_col[j] = inv_det[p_offset + j];
    }

    // Follows Sherman-Morrison update (branchless)
    for (std::size_t k = 0; k < N; ++k) {
        // Skip the special row for now
        if (k == particle)
            continue;

        const std::size_t k_offset{k * N}; // Pre-calculate row offset
        double s_k{};

// Compiler can now easily auto-vectorize this
#pragma omp simd reduction(+ : s_k)
        for (std::size_t m = 0; m < N; ++m) {
            s_k += new_row[m] * inv_det[k_offset + m];
        }

        const double factor{s_k * inv_ratio};

#pragma omp simd
        for (std::size_t j = 0; j < N; ++j) {
            inv_det[k_offset + j] -= inv_d_col[j] * factor;
        }
    }

// Special case: handle the particle row outside the loop
#pragma omp simd
    for (std::size_t j = 0; j < N; ++j) {
        inv_det[p_offset + j] = inv_d_col[j] * inv_ratio;
    }

    // Patch row `particle` of D to match the new positions:
#pragma omp simd
    for (std::size_t j = 0; j < N; ++j) {
        det_matrix[p_offset + j] = new_row[j];
    }
}

void SlaterPlaneWave::add_derivatives(double* RESTRICT grad_x, double* RESTRICT grad_y,
                                      double* RESTRICT grad_z,
                                      double* RESTRICT laplacian) const noexcept {
    const std::size_t N{num_orbitals_get()};

    const double* RESTRICT k_x{k_vector_x_get()};
    const double* RESTRICT k_y{k_vector_y_get()};
    const double* RESTRICT k_z{k_vector_z_get()};

    const auto& k_index{orbital_k_index_get()};
    const auto& o_type{orbital_type_get()};

    const double* RESTRICT inv_det{inv_determinant_get()};

    const std::size_t num_k{num_unique_k_get()};

    const double* RESTRICT cos_cache{cos_cache_get()};
    const double* RESTRICT sin_cache{sin_cache_get()};

    ASSUME_ALIGNED(k_x, SIMD_BYTES);
    ASSUME_ALIGNED(k_y, SIMD_BYTES);
    ASSUME_ALIGNED(k_z, SIMD_BYTES);

    ASSUME_ALIGNED(inv_det, SIMD_BYTES);
    ASSUME_ALIGNED(cos_cache, SIMD_BYTES);
    ASSUME_ALIGNED(sin_cache, SIMD_BYTES);

    ASSUME_ALIGNED(grad_x, SIMD_BYTES);
    ASSUME_ALIGNED(grad_y, SIMD_BYTES);
    ASSUME_ALIGNED(grad_z, SIMD_BYTES);
    ASSUME_ALIGNED(laplacian, SIMD_BYTES);

    for (std::size_t particle = 0; particle < N; ++particle) {
        double d_log_det_dx{}, d_log_det_dy{}, d_log_det_dz{};
        double laplace_det_term{};

        const std::size_t offset{particle * num_k};
        const std::size_t p_offset{particle * N};

        // Added reduction clauses so the compiler can safely vectorize the accumulators
#pragma omp simd reduction(+ : d_log_det_dx, d_log_det_dy, d_log_det_dz, laplace_det_term)
        for (std::size_t orbital = 0; orbital < N; ++orbital) {
            const std::size_t k_idx{k_index[orbital]};

            const double k_x_orbital{k_x[k_idx]};
            const double k_y_orbital{k_y[k_idx]};
            const double k_z_orbital{k_z[k_idx]};

            const double k_sq{k_x_orbital * k_x_orbital + k_y_orbital * k_y_orbital +
                              k_z_orbital * k_z_orbital};

            const double type{static_cast<double>(o_type[orbital])};

            const double cos_term{cos_cache[offset + k_idx]};
            const double sin_term{sin_cache[offset + k_idx]};

            const double grad_factor{-sin_term + type * (sin_term + cos_term)};
            const double lap_factor{-cos_term + type * (cos_term - sin_term)};

            const double weight{inv_det[p_offset + orbital]};

            d_log_det_dx += weight * k_x_orbital * grad_factor;
            d_log_det_dy += weight * k_y_orbital * grad_factor;
            d_log_det_dz += weight * k_z_orbital * grad_factor;
            laplace_det_term += weight * k_sq * lap_factor;
        }

        const double grad_sq{d_log_det_dx * d_log_det_dx + d_log_det_dy * d_log_det_dy +
                             d_log_det_dz * d_log_det_dz};

        grad_x[particle] += d_log_det_dx;
        grad_y[particle] += d_log_det_dy;
        grad_z[particle] += d_log_det_dz;

        laplacian[particle] += (laplace_det_term - grad_sq);
    }
}