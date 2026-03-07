#include "energy_tracking.hpp"

#include <cmath>
#include <numbers>

EnergyTracker::EnergyTracker(double box_length, double num_particles)
    : ewald_alpha_{6.0 / box_length},                                                       // 6.0 / L
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
    const double L{pbc.L_get()};

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
            const double dist{pbc.distance(p_x[i], p_y[i], p_z[i], p_x[j], p_y[j], p_z[j])};
            const bool degenerate{dist < 1.0e-12};
            const double mask{degenerate ? 0.0 : 1.0};

            const double r_ij{degenerate ? 0.0 : dist};
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

    const auto& G_vec_x{G_vector_x_get()};
    const auto& G_vec_y{G_vector_y_get()};
    const auto& G_vec_z{G_vector_z_get()};
    const auto& G_vec_weights{G_vector_weights_get()};

    double V_recip{};
    for (std::size_t i = 0; i < num_G; ++i) {
        const double g_x{G_vec_x[i]};
        const double g_y{G_vec_y[i]};
        const double g_z{G_vec_z[i]};
        const double g_w{G_vec_weights[i]};

        double cos_term{};
        double sin_term{};
        for (std::size_t j = 0; j < N; ++j) {
            const double G_dot_r{g_x * p_x[j] + g_y * p_y[j] + g_z * p_z[j]};
            cos_term += std::cos(G_dot_r);
            sin_term += std::sin(G_dot_r);
        }

        V_recip += g_w * (cos_term * cos_term + sin_term * sin_term);
    }
    V_recip *= prefactor;

    // Self + background:
    return V_real + V_recip + ewald_self_correction_term + ewald_background;
}

double EnergyTracker::eval_total_energy(const Particles& particles,
                                        const PeriodicBoundaryCondition& pbc) const noexcept {
    return kinetic_energy(particles) + potential_energy(particles, pbc);
}