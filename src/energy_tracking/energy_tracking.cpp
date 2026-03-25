#include "energy_tracking.hpp"

#include <numbers>
#include <omp.h>

EnergyTracker::EnergyTracker(double box_length, double num_particles)
    : box_length_{box_length}, ewald_alpha_{6.0 / box_length}, // 6.0 / L
      ewald_correction_{-6.0 * num_particles /
                        (std::sqrt(std::numbers::pi) * box_length)}, // -6.0 * N / (sqrt(pi) * L)
      ewald_background_{-std::numbers::pi * num_particles * num_particles /
                        (72.0 * box_length)}, // -pi * N^2 / (72L)
      V_recip_{}, V_real_{}, num_g_vectors_{}, data_{} {

    // Constants:
    const double two_pi_over_L{2.0 * std::numbers::pi / box_length};
    const double four_alpha_sq{4.0 * ewald_alpha_ * ewald_alpha_};
    const double cutoff_factor{-std::log(EWALD_RECIPROCAL_TOLERANCE)};

    // Maximums:
    const double g_max_mag_sq{four_alpha_sq * cutoff_factor};
    const int m_max{static_cast<int>(std::ceil(std::sqrt(g_max_mag_sq) / two_pi_over_L)) + 1};

    // Temp. Vectors:
    std::vector<double> tmp_x, tmp_y, tmp_z, tmp_w;

    auto& g_x{tmp_x};
    auto& g_y{tmp_y};
    auto& g_z{tmp_z};
    auto& g_weights{tmp_w};

    // G = 2*pi / L
    for (int m_x = -m_max; m_x <= m_max; ++m_x) {
        for (int m_y = -m_max; m_y <= m_max; ++m_y) {
            for (int m_z = -m_max; m_z <= m_max; ++m_z) {
                // Only keeps "positive half-sphere"
                // Make up for this by 2x the weights
                if (m_x < 0)
                    continue;
                if (m_x == 0 && m_y < 0)
                    continue;
                if (m_x == 0 && m_y == 0 && m_z <= 0)
                    continue;

                const double g_cand_x{two_pi_over_L * static_cast<double>(m_x)};
                const double g_cand_y{two_pi_over_L * static_cast<double>(m_y)};
                const double g_cand_z{two_pi_over_L * static_cast<double>(m_z)};
                const double g_cand_mag_sq{g_cand_x * g_cand_x + g_cand_y * g_cand_y +
                                           g_cand_z * g_cand_z};

                if (g_cand_mag_sq > g_max_mag_sq)
                    continue;

                g_x.emplace_back(g_cand_x);
                g_y.emplace_back(g_cand_y);
                g_z.emplace_back(g_cand_z);
                g_weights.emplace_back(8.0 * std::numbers::pi * std::numbers::pi / g_cand_mag_sq *
                                       std::exp(-g_cand_mag_sq / four_alpha_sq));
            }
        }
    }

    num_g_vectors_set() = g_x.size();

    data_ = AlignedSoA<double>(num_g_vectors_, NUM_ARRAYS_);

    std::copy_n(tmp_x.data(), num_g_vectors_, G_vector_x_get());
    std::copy_n(tmp_y.data(), num_g_vectors_, G_vector_y_get());
    std::copy_n(tmp_z.data(), num_g_vectors_, G_vector_z_get());
    std::copy_n(tmp_w.data(), num_g_vectors_, G_vector_weights_get());
}

void EnergyTracker::initialize_reciprocal_energy() noexcept {
    const double L{box_length_};
    const double prefactor{1.0 / (2.0 * std::numbers::pi * L * L * L)};

    const std::size_t num_G{num_g_vectors_get()};
    const double* RESTRICT g_weights{G_vector_weights_get()};
    const double* RESTRICT sum_real{sum_real_get()};
    const double* RESTRICT sum_imag{sum_imag_get()};

    ASSUME_ALIGNED(g_weights, SIMD_BYTES);
    ASSUME_ALIGNED(sum_real, SIMD_BYTES);
    ASSUME_ALIGNED(sum_imag, SIMD_BYTES);

    double sum{};
    #pragma omp simd reduction(+ : sum)
    for (std::size_t g = 0; g < num_G; ++g) {
        sum += g_weights[g] * (sum_real[g] * sum_real[g] + sum_imag[g] * sum_imag[g]);
    }
    V_recip_set() = prefactor * sum;
}

void EnergyTracker::initialize_real_energy(const Particles& particles) noexcept {
    const std::size_t N{particles.num_particles_get()};
    const double L{box_length_};
    const double neg_L{-1.0 * L};
    const double half_L{0.5 * L};
    const double neg_half_L{-1.0 * half_L};
    const double alpha{ewald_alpha_};

    // Fast erfc constants
    const double p{0.3275911};
    const double a1{0.254829592};
    const double a2{-0.284496736};
    const double a3{1.421413741};
    const double a4{-1.453152027};
    const double a5{1.061405429};

    const double* RESTRICT p_x{particles.pos_x_get()};
    const double* RESTRICT p_y{particles.pos_y_get()};
    const double* RESTRICT p_z{particles.pos_z_get()};

    ASSUME_ALIGNED(p_x, SIMD_BYTES);
    ASSUME_ALIGNED(p_y, SIMD_BYTES);
    ASSUME_ALIGNED(p_z, SIMD_BYTES);

    double sum{};
    #pragma omp parallel for reduction(+ : sum) schedule(dynamic)
    for (std::size_t i = 0; i < N; ++i) {
        double local_sum{};

        #pragma omp simd reduction(+ : local_sum)
        for (std::size_t j = i + 1; j < N; ++j) {
            double dx{p_x[i] - p_x[j]};
            double dy{p_y[i] - p_y[j]};
            double dz{p_z[i] - p_z[j]};

            dx += L * (dx <= neg_half_L) + neg_L * (dx > half_L);
            dy += L * (dy <= neg_half_L) + neg_L * (dy > half_L);
            dz += L * (dz <= neg_half_L) + neg_L * (dz > half_L);

            const double r{std::sqrt(dx * dx + dy * dy + dz * dz)};
            const double inv_r{(r < 1e-12) ? 1.0 : 1.0 / r};

            // Abramowitz & Stegun formula for fast std::erfc approx.
            const double erfc_arg{alpha * r};
            const double t{1.0 / (1.0 + p * erfc_arg)};
            const double tau{t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))))};

            local_sum += tau * std::exp(-erfc_arg * erfc_arg) * inv_r;
        }
        sum += local_sum;
    }
    V_real_set() = sum;
}

void EnergyTracker::initialize_structure_factors(const Particles& particles) noexcept {
    const std::size_t N{particles.num_particles_get()};
    const std::size_t num_G{num_g_vectors_get()};

    const double* RESTRICT p_x{particles.pos_x_get()};
    const double* RESTRICT p_y{particles.pos_y_get()};
    const double* RESTRICT p_z{particles.pos_z_get()};

    const double* RESTRICT g_x{G_vector_x_get()};
    const double* RESTRICT g_y{G_vector_y_get()};
    const double* RESTRICT g_z{G_vector_z_get()};

    double* RESTRICT sum_real{sum_real_get()};
    double* RESTRICT sum_imag{sum_imag_get()};

    ASSUME_ALIGNED(p_x, SIMD_BYTES);
    ASSUME_ALIGNED(p_y, SIMD_BYTES);
    ASSUME_ALIGNED(p_z, SIMD_BYTES);

    ASSUME_ALIGNED(g_x, SIMD_BYTES);
    ASSUME_ALIGNED(g_y, SIMD_BYTES);
    ASSUME_ALIGNED(g_z, SIMD_BYTES);

    ASSUME_ALIGNED(sum_real, SIMD_BYTES);
    ASSUME_ALIGNED(sum_imag, SIMD_BYTES);

    #pragma omp parallel for
    for (std::size_t g = 0; g < num_G; ++g) {
        double cos_sum{};
        double sin_sum{};

        #pragma omp simd reduction(+ : cos_sum, sin_sum)
        for (std::size_t j = 0; j < N; ++j) {
            const double G_dot_r{g_x[g] * p_x[j] + g_y[g] * p_y[j] + g_z[g] * p_z[j]};
            double cos_temp{};
            double sin_temp{};

            PORTABLE_SINCOS(G_dot_r, &sin_temp, &cos_temp);

            cos_sum += cos_temp;
            sin_sum += sin_temp;
        }

        sum_real[g] = cos_sum;
        sum_imag[g] = sin_sum;
    }
}

void EnergyTracker::update_structure_factors(double old_x, double old_y, double old_z, double new_x,
                                             double new_y, double new_z) noexcept {
    const double L{box_length_};

    const std::size_t num_G{num_g_vectors_get()};
    const double prefactor{1.0 / (2.0 * std::numbers::pi * L * L * L)};

    const double* RESTRICT g_weights{G_vector_weights_get()};
    const double* RESTRICT g_x{G_vector_x_get()};
    const double* RESTRICT g_y{G_vector_y_get()};
    const double* RESTRICT g_z{G_vector_z_get()};

    double* RESTRICT sum_real{sum_real_get()};
    double* RESTRICT sum_imag{sum_imag_get()};
    double* RESTRICT d_imag_temp{d_imag_temp_get()};
    double* RESTRICT d_real_temp{d_real_temp_get()};

    ASSUME_ALIGNED(g_weights, SIMD_BYTES);
    ASSUME_ALIGNED(g_x, SIMD_BYTES);
    ASSUME_ALIGNED(g_y, SIMD_BYTES);
    ASSUME_ALIGNED(g_z, SIMD_BYTES);

    ASSUME_ALIGNED(sum_real, SIMD_BYTES);
    ASSUME_ALIGNED(sum_imag, SIMD_BYTES);
    ASSUME_ALIGNED(d_imag_temp, SIMD_BYTES);
    ASSUME_ALIGNED(d_real_temp, SIMD_BYTES);

    double delta{};
    
    #pragma omp simd
    for (std::size_t g = 0; g < num_G; ++g) {
        const double old_dot{g_x[g] * old_x + g_y[g] * old_y + g_z[g] * old_z};
        const double new_dot{g_x[g] * new_x + g_y[g] * new_y + g_z[g] * new_z};

        double new_sin{}, new_cos{};
        double old_sin{}, old_cos{};

        PORTABLE_SINCOS(new_dot, &new_sin, &new_cos);
        PORTABLE_SINCOS(old_dot, &old_sin, &old_cos);

        d_real_temp[g] = new_cos - old_cos;
        d_imag_temp[g] = new_sin - old_sin;
    }

   // Accumulate delta and update sum_real / sum_imag
    #pragma omp simd reduction(+ : delta)
    for (std::size_t g = 0; g < num_G; ++g) {
        const double dr{d_real_temp[g]};
        const double di{d_imag_temp[g]};

        delta += g_weights[g] * (2.0 * (sum_real[g] * dr + sum_imag[g] * di) +
                                 dr * dr + di * di);

        sum_real[g] += dr;
        sum_imag[g] += di;
    }

    V_recip_set() += prefactor * delta;
}

double EnergyTracker::kinetic_energy(const Particles& particles) const noexcept {
    const double* RESTRICT grad_x{particles.grad_log_psi_x_get()};
    const double* RESTRICT grad_y{particles.grad_log_psi_y_get()};
    const double* RESTRICT grad_z{particles.grad_log_psi_z_get()};
    const double* RESTRICT lap{particles.lap_log_psi_get()};

    ASSUME_ALIGNED(grad_x, SIMD_BYTES);
    ASSUME_ALIGNED(grad_y, SIMD_BYTES);
    ASSUME_ALIGNED(grad_z, SIMD_BYTES);
    ASSUME_ALIGNED(lap, SIMD_BYTES);

    // Kinetic
    double T_sum{};
    const std::size_t N{particles.num_particles_get()};

    #pragma omp simd reduction(+ : T_sum)
    for (std::size_t i = 0; i < N; ++i) {
        // Computes ||Grad(logPsi)||^2
        const double grad_sq{grad_x[i] * grad_x[i] + grad_y[i] * grad_y[i] + grad_z[i] * grad_z[i]};

        // Accumulate Lapl(LogPsi) + ||Grad(LogPsi)||^2
        T_sum += (lap[i] + grad_sq);
    }

    return -0.5 * T_sum;
}

void EnergyTracker::update_real_energy(std::size_t moved_idx, double old_x, double old_y,
                                       double old_z, const Particles& particles) noexcept {
    const std::size_t N{particles.num_particles_get()};
    const double L{box_length_};
    const double neg_L{-1.0 * L};
    const double half_L{0.5 * L};
    const double neg_half_L{-1.0 * half_L};

    const double alpha{ewald_alpha_};

    // Fast erfc constants
    const double p{0.3275911};
    const double a1{0.254829592};
    const double a2{-0.284496736};
    const double a3{1.421413741};
    const double a4{-1.453152027};
    const double a5{1.061405429};

    const double* RESTRICT p_x{particles.pos_x_get()};
    const double* RESTRICT p_y{particles.pos_y_get()};
    const double* RESTRICT p_z{particles.pos_z_get()};

    ASSUME_ALIGNED(p_x, SIMD_BYTES);
    ASSUME_ALIGNED(p_y, SIMD_BYTES);
    ASSUME_ALIGNED(p_z, SIMD_BYTES);

    const double new_x{p_x[moved_idx]};
    const double new_y{p_y[moved_idx]};
    const double new_z{p_z[moved_idx]};

    double delta{};

    #pragma omp simd reduction(+ : delta)
    for (std::size_t j = 0; j < N; ++j) {
        // Branchless mask to safely skip the moved particle
        const double valid_mask{(j == moved_idx) ? 0.0 : 1.0};

        // Old pair
        double dx_old{old_x - p_x[j]};
        double dy_old{old_y - p_y[j]};
        double dz_old{old_z - p_z[j]};

        dx_old += L * (dx_old <= neg_half_L) + neg_L * (dx_old > half_L);
        dy_old += L * (dy_old <= neg_half_L) + neg_L * (dy_old > half_L);
        dz_old += L * (dz_old <= neg_half_L) + neg_L * (dz_old > half_L);

        const double r_old{std::sqrt(dx_old * dx_old + dy_old * dy_old + dz_old * dz_old)};

        // Protect against 1.0 / 0.0 generating NaN
        const double inv_r_old{(r_old < 1e-12) ? 1.0 : 1.0 / r_old};

        // Abramowitz & Stegun formula for fast std::erfc approx.
        const double erfc_arg_old{alpha * r_old};
        const double t_old{1.0 / (1.0 + p * erfc_arg_old)};
        const double tau_old{t_old * (a1 + t_old * (a2 + t_old * (a3 + t_old * (a4 + t_old * a5))))};

        double const erfc_old{tau_old * std::exp(-erfc_arg_old * erfc_arg_old) * inv_r_old};

        // New pair
        double dx_new{new_x - p_x[j]};
        double dy_new{new_y - p_y[j]};
        double dz_new{new_z - p_z[j]};

        dx_new += L * (dx_new <= neg_half_L) + neg_L * (dx_new > half_L);
        dy_new += L * (dy_new <= neg_half_L) + neg_L * (dy_new > half_L);
        dz_new += L * (dz_new <= neg_half_L) + neg_L * (dz_new > half_L);

        const double r_new{std::sqrt(dx_new * dx_new + dy_new * dy_new + dz_new * dz_new)};

        // Protect against 1.0 / 0.0 generating NaN
        const double inv_r_new{(r_new < 1e-12) ? 1.0 : 1.0 / r_new};

        // Combine operations and apply the mask

        // Abramowitz & Stegun formula for fast std::erfc approx.
        const double erfc_arg_new{alpha * r_new};
        const double t_new{1.0 / (1.0 + p * erfc_arg_new)};
        const double tau_new{t_new * (a1 + t_new * (a2 + t_new * (a3 + t_new * (a4 + t_new * a5))))};

        double const erfc_new{tau_new * std::exp(-erfc_arg_new * erfc_arg_new) * inv_r_new};

        delta += valid_mask * (erfc_new - erfc_old);
    }

    V_real_ += delta;
}

double EnergyTracker::potential_energy() const noexcept {
    // Ewald constants:
    const double ewald_self_correction_term{ewald_correction_get()};
    const double ewald_background{ewald_background_get()};

    // Potentials calcualted in cache:
    const double V_recip{V_recip_get()};
    const double V_real{V_real_get()};

    // Self + background:
    return V_real + V_recip + ewald_self_correction_term + ewald_background;
}

double EnergyTracker::eval_total_energy(const Particles& particles) const noexcept {
    return kinetic_energy(particles) + potential_energy();
}
