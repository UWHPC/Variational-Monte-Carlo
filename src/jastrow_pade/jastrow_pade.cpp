#include "jastrow_pade.hpp"

#include <cmath>
#include <cstddef>
#include <omp.h>

double JastrowPade::value(const Particles& particles) const noexcept {
    const std::size_t num_particles{particles.num_particles_get()};
    const double L{box_length_};
    const double neg_L{-1.0 * L};
    const double half_L{0.5 * L};
    const double neg_half_L{-1.0 * half_L};

    const double* RESTRICT pos_x{particles.pos_x_get()};
    const double* RESTRICT pos_y{particles.pos_y_get()};
    const double* RESTRICT pos_z{particles.pos_z_get()};

    const double a_local{a_get()};
    const double b_local{b_get()};

    double jastrow_pade{};

    for (std::size_t i = 0; i < num_particles; ++i) {
        double local_jastrow{};

#pragma omp simd reduction(+ : local_jastrow)
        for (std::size_t j = 0; j < num_particles; ++j) {
            const double mask{i == j ? 0.0 : 1.0};

            double displ_x{pos_x[i] - pos_x[j]};
            double displ_y{pos_y[i] - pos_y[j]};
            double displ_z{pos_z[i] - pos_z[j]};

            // Boolean masks - reduce to 0 if false, and 1 if true.
            displ_x += L * (displ_x <= neg_half_L) + neg_L * (displ_x > half_L);
            displ_y += L * (displ_y <= neg_half_L) + neg_L * (displ_y > half_L);
            displ_z += L * (displ_z <= neg_half_L) + neg_L * (displ_z > half_L);

            const double dist_sq{displ_x * displ_x + displ_y * displ_y + displ_z * displ_z};
            const double dist{std::sqrt(dist_sq)};

            // u(r) = a*r / (1 + b*r)
            const double denom{1.0 + b_local * dist};
            const double inv_denom{1.0 / denom};

            local_jastrow += mask * a_local * dist * inv_denom;
        }
        jastrow_pade += local_jastrow;
    }
    // Computed full N x N matrix - jastrow pade was calculated twice
    return 0.5 * jastrow_pade;
}

void JastrowPade::add_derivatives(const Particles& particles, double* RESTRICT grad_x,
                                  double* RESTRICT grad_y, double* RESTRICT grad_z,
                                  double* RESTRICT laplacian) const noexcept {
    // NOTE: assumes gradX/gradY/gradZ/lap are zero-initialized by caller
    const std::size_t num_particles{particles.num_particles_get()};
    const double L{box_length_};
    const double neg_L{-1.0 * L};
    const double half_L{0.5 * L};
    const double neg_half_L{-1.0 * half_L};

    const double* RESTRICT pos_x{particles.pos_x_get()};
    const double* RESTRICT pos_y{particles.pos_y_get()};
    const double* RESTRICT pos_z{particles.pos_z_get()};

    const double a_local{a_get()};
    const double b_local{b_get()};
    const double neg_two_a_b{-2.0 * a_local * b_local};

    for (std::size_t i = 0; i < num_particles; ++i) {
        // i accumulators to prevent SIMD lane collisions
        double d_grad_x{}, d_grad_y{}, d_grad_z{}, d_lap{};

#pragma omp simd reduction(+ : d_grad_x, d_grad_y, d_grad_z, d_lap)
        for (std::size_t j = 0; j < num_particles; ++j) {
            const double valid_idx{i == j ? 0.0 : 1.0};

            double displ_x{pos_x[i] - pos_x[j]};
            double displ_y{pos_y[i] - pos_y[j]};
            double displ_z{pos_z[i] - pos_z[j]};

            // Boolean masks - reduce to 0 if false, and 1 if true.
            displ_x += L * (displ_x <= neg_half_L) + neg_L * (displ_x > half_L);
            displ_y += L * (displ_y <= neg_half_L) + neg_L * (displ_y > half_L);
            displ_z += L * (displ_z <= neg_half_L) + neg_L * (displ_z > half_L);

            const double dist_sq{displ_x * displ_x + displ_y * displ_y + displ_z * displ_z};
            const double dist{std::sqrt(dist_sq)};

            // mask to get around if statement:
            const bool degenerate{dist < 1e-12};
            const double inv_dist{degenerate ? 1.0 : 1.0 / dist};
            const double mask{degenerate ? 0.0 : valid_idx};

            // u(r) = a*r / (1 + b*r)
            // u'(r) = a / (1 + b*r)^2
            // u''(r) = -2ab / (1 + b*r)^3
            const double denom{1.0 / (1.0 + b_local * dist)};
            const double denom_sq{denom * denom};
            const double denom_cb{denom_sq * denom};

            const double first_deriv{a_local * denom_sq};
            const double second_deriv{neg_two_a_b * denom_cb};

            // ∇_i u(r_ij) = u'(r) * (r_vec / r)
            const double grad_factor{mask * first_deriv * inv_dist};

            // ∇^2 u(r) = u''(r) + (2/r) u'(r)
            const double laplacian_pair{mask * (second_deriv + 2.0 * first_deriv * inv_dist)};

            d_grad_x += grad_factor * displ_x;
            d_grad_y += grad_factor * displ_y;
            d_grad_z += grad_factor * displ_z;

            d_lap += laplacian_pair;
        }

        // Apply the i values once per row
        grad_x[i] += d_grad_x;
        grad_y[i] += d_grad_y;
        grad_z[i] += d_grad_z;
        laplacian[i] += d_lap;
    }
}

double JastrowPade::delta_value(const Particles& particles, std::size_t moved, double old_x,
                                double old_y, double old_z) const noexcept {
    const std::size_t num_particles{particles.num_particles_get()};
    const double L{box_length_};
    const double neg_L{-1.0 * L};
    const double half_L{0.5 * L};
    const double neg_half_L{-1.0 * half_L};

    const double* RESTRICT pos_x{particles.pos_x_get()};
    const double* RESTRICT pos_y{particles.pos_y_get()};
    const double* RESTRICT pos_z{particles.pos_z_get()};

    const double a_local{a_get()};
    const double b_local{b_get()};

    const double new_x{pos_x[moved]};
    const double new_y{pos_y[moved]};
    const double new_z{pos_z[moved]};

    double delta{};

#pragma omp simd reduction(+ : delta)
    for (std::size_t j = 0; j < num_particles; ++j) {
        // mask to skip the moved particle safely
        const double valid_mask{(j == moved) ? 0.0 : 1.0};

        // Old pair:
        double displ_old_x{old_x - pos_x[j]};
        double displ_old_y{old_y - pos_y[j]};
        double displ_old_z{old_z - pos_z[j]};

        displ_old_x += L * (displ_old_x <= neg_half_L) + neg_L * (displ_old_x > half_L);
        displ_old_y += L * (displ_old_y <= neg_half_L) + neg_L * (displ_old_y > half_L);
        displ_old_z += L * (displ_old_z <= neg_half_L) + neg_L * (displ_old_z > half_L);

        const double dist_old{std::sqrt(displ_old_x * displ_old_x + displ_old_y * displ_old_y +
                                        displ_old_z * displ_old_z)};
        const double denom_old{1.0 / (1.0 + b_local * dist_old)};

        // New pair:
        double displ_new_x{new_x - pos_x[j]};
        double displ_new_y{new_y - pos_y[j]};
        double displ_new_z{new_z - pos_z[j]};

        displ_new_x += L * (displ_new_x <= neg_half_L) + neg_L * (displ_new_x > half_L);
        displ_new_y += L * (displ_new_y <= neg_half_L) + neg_L * (displ_new_y > half_L);
        displ_new_z += L * (displ_new_z <= neg_half_L) + neg_L * (displ_new_z > half_L);

        const double dist_new{std::sqrt(displ_new_x * displ_new_x + displ_new_y * displ_new_y +
                                        displ_new_z * displ_new_z)};
        const double denom_new{1.0 / (1.0 + b_local * dist_new)};

        delta += valid_mask * a_local * (dist_new * denom_new - dist_old * denom_old);
    }

    return delta;
}

// Used to incrementally update, faster than recomputing
void JastrowPade::update_derivatives_for_move(const Particles& particles, std::size_t moved,
                                              double old_x, double old_y, double old_z,
                                              double* RESTRICT grad_x, double* RESTRICT grad_y,
                                              double* RESTRICT grad_z,
                                              double* RESTRICT laplacian) const noexcept {
    grad_x[moved] = 0.0;
    grad_y[moved] = 0.0;
    grad_z[moved] = 0.0;
    laplacian[moved] = 0.0;

    const std::size_t num_particles{particles.num_particles_get()};
    const double L{box_length_};
    const double neg_L{-1.0 * L};

    const double half_L{0.5 * L};
    const double neg_half_L{-1.0 * half_L};

    const double* RESTRICT pos_x{particles.pos_x_get()};
    const double* RESTRICT pos_y{particles.pos_y_get()};
    const double* RESTRICT pos_z{particles.pos_z_get()};

    const double a_local{a_get()};
    const double b_local{b_get()};
    const double m2ab{-2.0 * a_local * b_local};

    const double new_x{pos_x[moved]};
    const double new_y{pos_y[moved]};
    const double new_z{pos_z[moved]};

    // Moved variables out of the loop
    double m_grad_x{}, m_grad_y{}, m_grad_z{}, m_lap{};

#pragma omp simd reduction(+ : m_grad_x, m_grad_y, m_grad_z, m_lap)
    for (std::size_t j = 0; j < num_particles; ++j) {
        // Combined branchless masks
        const bool is_moved{j == moved};

        // Old pair:
        double displ_old_x{old_x - pos_x[j]};
        double displ_old_y{old_y - pos_y[j]};
        double displ_old_z{old_z - pos_z[j]};

        displ_old_x += L * (displ_old_x <= neg_half_L) + neg_L * (displ_old_x > half_L);
        displ_old_y += L * (displ_old_y <= neg_half_L) + neg_L * (displ_old_y > half_L);
        displ_old_z += L * (displ_old_z <= neg_half_L) + neg_L * (displ_old_z > half_L);

        const double dist_old{std::sqrt(displ_old_x * displ_old_x + displ_old_y * displ_old_y +
                                        displ_old_z * displ_old_z)};

        const double inv_dist_old{(dist_old < 1e-12) ? 1.0 : 1.0 / dist_old};
        const double mask_old{(is_moved || dist_old < 1e-12) ? 0.0 : 1.0};

        const double inv_denom_old{1.0 / (1.0 + b_local * dist_old)};
        const double inv_denom_sq_old{inv_denom_old * inv_denom_old};

        const double first_deriv_old{a_local * inv_denom_sq_old};
        const double second_deriv_old{m2ab * inv_denom_sq_old * inv_denom_old};

        const double grad_factor_old{mask_old * first_deriv_old * inv_dist_old};
        const double lap_pair_old{mask_old *
                                  (second_deriv_old + 2.0 * first_deriv_old * inv_dist_old)};

        // New pair:
        double displ_new_x{new_x - pos_x[j]};
        double displ_new_y{new_y - pos_y[j]};
        double displ_new_z{new_z - pos_z[j]};

        displ_new_x += L * (displ_new_x <= neg_half_L) + neg_L * (displ_new_x > half_L);
        displ_new_y += L * (displ_new_y <= neg_half_L) + neg_L * (displ_new_y > half_L);
        displ_new_z += L * (displ_new_z <= neg_half_L) + neg_L * (displ_new_z > half_L);

        const double dist_new{std::sqrt(displ_new_x * displ_new_x + displ_new_y * displ_new_y +
                                        displ_new_z * displ_new_z)};

        const double inv_dist_new{(dist_new < 1e-12) ? 1.0 : 1.0 / dist_new};
        const double mask_new{(is_moved || dist_new < 1e-12) ? 0.0 : 1.0};

        const double inv_denom_new{1.0 / (1.0 + b_local * dist_new)};
        const double inv_denom_sq_new{inv_denom_new * inv_denom_new};

        const double first_deriv_new{a_local * inv_denom_sq_new};
        const double second_deriv_new{m2ab * inv_denom_sq_new * inv_denom_new};

        const double grad_factor_new{mask_new * first_deriv_new * inv_dist_new};
        const double lap_pair_new{mask_new *
                                  (second_deriv_new + 2.0 * first_deriv_new * inv_dist_new)};

        // Apply accumulations safely
        m_grad_x += grad_factor_new * displ_new_x;
        m_grad_y += grad_factor_new * displ_new_y;
        m_grad_z += grad_factor_new * displ_new_z;
        m_lap += lap_pair_new;

        // Combine the j-array updates to save memory writes
        grad_x[j] += grad_factor_old * displ_old_x - grad_factor_new * displ_new_x;
        grad_y[j] += grad_factor_old * displ_old_y - grad_factor_new * displ_new_y;
        grad_z[j] += grad_factor_old * displ_old_z - grad_factor_new * displ_new_z;

        laplacian[j] += lap_pair_new - lap_pair_old;
    }

    grad_x[moved] = m_grad_x;
    grad_y[moved] = m_grad_y;
    grad_z[moved] = m_grad_z;
    laplacian[moved] = m_lap;
}