#include "jastrow_pade.hpp"

#include <cmath>
#include <cstddef>

double JastrowPade::value(const Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept {
    const std::size_t num_particles{particles.num_particles_get()};
    const double L{box_length_};
    const double half_L{0.5 * L};
    const double inv_L{1.0 / L};

    const double* RESTRICT pos_x{particles.pos_x_get()};
    const double* RESTRICT pos_y{particles.pos_y_get()};
    const double* RESTRICT pos_z{particles.pos_z_get()};

    const double a_local{a_get()};
    const double b_local{b_get()};

    double jastrow_pade{};

    for (std::size_t i = 0; i < num_particles; ++i) {
        for (std::size_t j = i + 1; j < num_particles; ++j) {
            double displ_x{pos_x[i] - pos_x[j]};
            double displ_y{pos_y[i] - pos_y[j]};
            double displ_z{pos_z[i] - pos_z[j]};

            displ_x -= L * std::round(displ_x * inv_L);
            displ_y -= L * std::round(displ_y * inv_L);
            displ_z -= L * std::round(displ_z * inv_L);

            // Boolean masks - reduce to 0 if false, and 1 if true.
            displ_x += L * (displ_x <= -half_L) - L * (displ_x > half_L);
            displ_y += L * (displ_y <= -half_L) - L * (displ_y > half_L);
            displ_z += L * (displ_z <= -half_L) - L * (displ_z > half_L);

            const double dist_sq{displ_x * displ_x + displ_y * displ_y + displ_z * displ_z};
            const double dist{std::sqrt(dist_sq)};

            // u(r) = a*r / (1 + b*r)
            const double denom{1.0 + b_local * dist};
            const double inv_denom{1.0 / denom};

            jastrow_pade += a_local * dist * inv_denom;
        }
    }
    return jastrow_pade;
}

void JastrowPade::add_derivatives(const Particles& particles, const PeriodicBoundaryCondition& pbc,
                                  double* RESTRICT grad_x, double* RESTRICT grad_y, double* RESTRICT grad_z,
                                  double* RESTRICT laplacian) const noexcept {
    // NOTE: assumes gradX/gradY/gradZ/lap are zero-initialized by caller
    const std::size_t num_particles{particles.num_particles_get()};
    const double L{box_length_};
    const double half_L{0.5 * L};
    const double inv_L{1.0 / L};

    const double* RESTRICT pos_x{particles.pos_x_get()};
    const double* RESTRICT pos_y{particles.pos_y_get()};
    const double* RESTRICT pos_z{particles.pos_z_get()};

    const double a_local{a_get()};
    const double b_local{b_get()};

    for (std::size_t i = 0; i < num_particles; ++i) {
        for (std::size_t j = i + 1; j < num_particles; ++j) {
            double displ_x{pos_x[i] - pos_x[j]};
            double displ_y{pos_y[i] - pos_y[j]};
            double displ_z{pos_z[i] - pos_z[j]};

            displ_x -= L * std::round(displ_x * inv_L);
            displ_y -= L * std::round(displ_y * inv_L);
            displ_z -= L * std::round(displ_z * inv_L);

            // Boolean masks - reduce to 0 if false, and 1 if true.
            displ_x += L * (displ_x <= -half_L) - L * (displ_x > half_L);
            displ_y += L * (displ_y <= -half_L) - L * (displ_y > half_L);
            displ_z += L * (displ_z <= -half_L) - L * (displ_z > half_L);

            const double dist_sq{displ_x * displ_x + displ_y * displ_y + displ_z * displ_z};
            const double dist{std::sqrt(dist_sq)};

            // mask to get around if statement:
            const bool degenerate{dist < 1e-12};
            const double inv_dist{degenerate ? 1.0 : 1.0 / dist};
            const double mask{degenerate ? 0.0 : 1.0};

            // u(r) = a*r / (1 + b*r)
            // u'(r) = a / (1 + b*r)^2
            // u''(r) = -2ab / (1 + b*r)^3
            const double denom{1.0 + b_local * dist};
            const double denom_sq{denom * denom};
            const double denom_cb{denom_sq * denom};

            const double first_deriv{a_local / denom_sq};
            const double second_deriv{-2.0 * a_local * b_local / denom_cb};

            // ∇_i u(r_ij) = u'(r) * (r_vec / r)
            const double grad_factor{first_deriv * inv_dist};
            const double grad_contribution_x{grad_factor * displ_x};
            const double grad_contribution_y{grad_factor * displ_y};
            const double grad_contribution_z{grad_factor * displ_z};

            grad_x[i] += mask * grad_contribution_x;
            grad_y[i] += mask * grad_contribution_y;
            grad_z[i] += mask * grad_contribution_z;

            grad_x[j] -= mask * grad_contribution_x;
            grad_y[j] -= mask * grad_contribution_y;
            grad_z[j] -= mask * grad_contribution_z;

            // ∇^2 u(r) = u''(r) + (2/r) u'(r)
            const double laplacian_pair{second_deriv + 2.0 * first_deriv * inv_dist};

            laplacian[i] += mask * laplacian_pair;
            laplacian[j] += mask * laplacian_pair;
        }
    }
}
