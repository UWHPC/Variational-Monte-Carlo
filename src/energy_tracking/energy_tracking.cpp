#include "energy_tracking.hpp"

#include <cmath>
#include <numbers>

EnergyTracker::EnergyTracker(double box_length, double num_particles)
    : box_length_{box_length}, ewald_alpha_{6.0 / box_length},                              // 6.0 / L
      ewald_correction_{-6.0 * num_particles / (std::sqrt(std::numbers::pi) * box_length)}, // -6.0 * N / (sqrt(pi) * L)
      ewald_background_{-std::numbers::pi * num_particles * num_particles / (72.0 * box_length)} { // -pi * N^2 / (72L)

    // Constants:
    const double two_pi_over_L{2.0 * std::numbers::pi / box_length};
    const double four_alpha_sq{4.0 * ewald_alpha_ * ewald_alpha_};
    const double cutoff_factor{15.0};

    // Maximums:
    const double g_max_mag_sq{four_alpha_sq * cutoff_factor};
    const int m_max{static_cast<int>(std::ceil(std::sqrt(g_max_mag_sq) / two_pi_over_L)) + 1};

    auto& g_x{G_vector_x_get()};
    auto& g_y{G_vector_y_get()};
    auto& g_z{G_vector_z_get()};
    auto& g_weights{G_vector_weights_get()};

    // G = 2*pi / L
    for (int m_x = -m_max; m_x <= m_max; ++m_x) {
        for (int m_y = -m_max; m_y <= m_max; ++m_y) {
            for (int m_z = -m_max; m_z <= m_max; ++m_z) {
                if (m_x == 0 && m_y == 0 && m_z == 0)
                    continue;

                const double g_cand_x{two_pi_over_L * static_cast<double>(m_x)};
                const double g_cand_y{two_pi_over_L * static_cast<double>(m_y)};
                const double g_cand_z{two_pi_over_L * static_cast<double>(m_z)};
                const double g_cand_mag_sq{g_cand_x * g_cand_x + g_cand_y * g_cand_y + g_cand_z * g_cand_z};

                if (g_cand_mag_sq > g_max_mag_sq)
                    continue;

                g_x.emplace_back(g_cand_x);
                g_y.emplace_back(g_cand_y);
                g_z.emplace_back(g_cand_z);
                g_weights.emplace_back(4.0 * std::numbers::pi * std::numbers::pi / g_cand_mag_sq *
                                       std::exp(-g_cand_mag_sq / four_alpha_sq));
            }
        }
    }

    // Size structure factor cache to match G-vector count:
    sum_real_get().resize(num_g_vectors_get(), 0.0);
    sum_imag_get().resize(num_g_vectors_get(), 0.0);
}

void EnergyTracker::initialize_structure_factors(const Particles& particles) noexcept {
    const std::size_t N{particles.num_particles_get()};
    const std::size_t num_G{num_g_vectors_get()};

    const double* RESTRICT p_x{particles.pos_x_get()};
    const double* RESTRICT p_y{particles.pos_y_get()};
    const double* RESTRICT p_z{particles.pos_z_get()};

    const auto& g_x{G_vector_x_get()};
    const auto& g_y{G_vector_y_get()};
    const auto& g_z{G_vector_z_get()};

    auto& sum_real{sum_real_get()};
    auto& sum_imag{sum_imag_get()};

    for (std::size_t g = 0; g < num_G; ++g) {
        double cos_sum{};
        double sin_sum{};

        for (std::size_t j = 0; j < N; ++j) {
            const double G_dot_r{g_x[g] * p_x[j] + g_y[g] * p_y[j] + g_z[g] * p_z[j]};

            cos_sum += std::cos(G_dot_r);
            sin_sum += std::sin(G_dot_r);
        }

        sum_real[g] = cos_sum;
        sum_imag[g] = sin_sum;
    }
}

void EnergyTracker::update_structure_factors(double old_x, double old_y, double old_z, double new_x, double new_y,
                                             double new_z) noexcept {
    const std::size_t num_G{num_g_vectors_get()};

    const auto& g_x{G_vector_x_get()};
    const auto& g_y{G_vector_y_get()};
    const auto& g_z{G_vector_z_get()};

    auto& sum_real{sum_real_get()};
    auto& sum_imag{sum_imag_get()};

    for (std::size_t g = 0; g < num_G; ++g) {
        const double old_dot{g_x[g] * old_x + g_y[g] * old_y + g_z[g] * old_z};
        const double new_dot{g_x[g] * new_x + g_y[g] * new_y + g_z[g] * new_z};

        // Subtract old contribution, add new:
        sum_real[g] += std::cos(new_dot) - std::cos(old_dot);
        sum_imag[g] += std::sin(new_dot) - std::sin(old_dot);
    }
}

double EnergyTracker::kinetic_energy(const Particles& particles) const noexcept {
    const double* RESTRICT grad_x{particles.grad_log_psi_x_get()};
    const double* RESTRICT grad_y{particles.grad_log_psi_y_get()};
    const double* RESTRICT grad_z{particles.grad_log_psi_z_get()};
    const double* RESTRICT lap{particles.lap_log_psi_get()};

    // Kinetic
    double T_sum{};
    const std::size_t N{particles.num_particles_get()};

    for (std::size_t i = 0; i < N; ++i) {
        // Computes ||Grad(logPsi)||^2
        const double grad_sq{grad_x[i] * grad_x[i] + grad_y[i] * grad_y[i] + grad_z[i] * grad_z[i]};

        // Accumulate Lapl(LogPsi) + ||Grad(LogPsi)||^2
        T_sum += (lap[i] + grad_sq);
    }

    return -0.5 * T_sum;
}

double EnergyTracker::potential_energy(const Particles& particles,
                                       const PeriodicBoundaryCondition& pbc) const noexcept {
    const std::size_t N{particles.num_particles_get()};
    const double L{box_length_};
    const double half_L{0.5 * L};
    const double inv_L{1.0 / L};

    const double* RESTRICT p_x{particles.pos_x_get()};
    const double* RESTRICT p_y{particles.pos_y_get()};
    const double* RESTRICT p_z{particles.pos_z_get()};

    const double ewald_alpha{ewald_alpha_get()};
    const double ewald_self_correction_term{ewald_correction_get()};
    const double ewald_background{ewald_background_get()};

    // Real-space sum: erfc(a * r_ij) / r_ij
    double V_real{};
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i + 1; j < N; ++j) {
            // To get around if else statement for branching,
            // ensures does not run if r_ij < 1.0e-12
            double displ_x{p_x[i] - p_x[j]};
            double displ_y{p_y[i] - p_y[j]};
            double displ_z{p_z[i] - p_z[j]};

            displ_x -= L * std::round(displ_x * inv_L);
            displ_y -= L * std::round(displ_y * inv_L);
            displ_z -= L * std::round(displ_z * inv_L);

            // Boolean masks - reduce to 0 if false, and 1 if true.
            displ_x += L * (displ_x <= -half_L) - L * (displ_x > half_L);
            displ_y += L * (displ_y <= -half_L) - L * (displ_y > half_L);
            displ_z += L * (displ_z <= -half_L) - L * (displ_z > half_L);

            const double dist_sq{displ_x * displ_x + displ_y * displ_y + displ_z * displ_z};
            const double dist{std::sqrt(dist_sq)};

            const bool degenerate{dist < 1.0e-12};
            const double mask{degenerate ? 0.0 : 1.0};

            const double r_ij{dist};
            const double inv_r_ij{degenerate ? 1.0 : 1 / r_ij};

            // if degenerate == true,
            // V_sum = 0.0 / 1.0 = 0
            // so the contribution is 0
            V_real += mask * std::erfc(ewald_alpha * r_ij) * inv_r_ij;
        }
    }

    // Reciprocal-space sum:
    const double prefactor{1.0 / (2.0 * std::numbers::pi * L * L * L)};
    const std::size_t num_G{num_g_vectors_get()};

    const auto& G_vec_weights{G_vector_weights_get()};
    auto& sum_real{sum_real_get()};
    auto& sum_imag{sum_imag_get()};

    double V_recip{};
    for (std::size_t i = 0; i < num_G; ++i) {
        V_recip += G_vec_weights[i] * (sum_real[i] * sum_real[i] + sum_imag[i] * sum_imag[i]);
    }
    V_recip *= prefactor;

    // Self + background:
    return V_real + V_recip + ewald_self_correction_term + ewald_background;
}

double EnergyTracker::eval_total_energy(const Particles& particles,
                                        const PeriodicBoundaryCondition& pbc) const noexcept {
    return kinetic_energy(particles) + potential_energy(particles, pbc);
}