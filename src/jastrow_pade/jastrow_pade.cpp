#include "jastrow_pade.hpp"

#include <cmath>
#include <cstddef>

double JastrowPade::value(const Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept {
    const std::size_t num_particles{particles.num_particles_ptr()};

    const double* RESTRICT POS_X{particles.pos_x_ptr()};
    const double* RESTRICT POS_Y{particles.pos_y_ptr()};
    const double* RESTRICT POS_Z{particles.pos_z_ptr()};

    const double A_LOCAL{a_ptr()};
    const double B_LOCAL{b_ptr()};

    double jastrow_pade{0.0};

    for (std::size_t i = 0; i < num_particles; ++i) {
        for (std::size_t j = i + 1; j < num_particles; ++j) {
            double displ_x{}, displ_y{}, displ_z{};

            pbc.displacement(POS_X[i], POS_Y[i], POS_Z[i], POS_X[j], POS_Y[j], POS_Z[j], displ_x, displ_y, displ_z);
            const double DIST_SQ{displ_x * displ_x + displ_y * displ_y + displ_z * displ_z};
            const double DIST{std::sqrt(DIST_SQ)};

            // MASK to get around if statement
            const bool DEGENERATE{DIST < 1e-12};
            const double MASK{DEGENERATE ? 0.0 : 1.0};

            // u(r) = a*r / (1 + b*r)
            const double DENOM{1.0 + B_LOCAL * DIST};

            jastrow_pade += MASK * (A_LOCAL * DIST) / DENOM;
        }
    }
    return jastrow_pade;
}

void JastrowPade::add_derivatives(const Particles& particles, const PeriodicBoundaryCondition& pbc,
                                  double* RESTRICT grad_x, double* RESTRICT grad_y, double* RESTRICT grad_z,
                                  double* RESTRICT laplacian) const noexcept {
    // NOTE: assumes gradX/gradY/gradZ/lap are zero-initialized by caller
    const std::size_t num_particles{particles.num_particles_ptr()};

    const double* RESTRICT POS_X{particles.pos_x_ptr()};
    const double* RESTRICT POS_Y{particles.pos_y_ptr()};
    const double* RESTRICT POS_Z{particles.pos_z_ptr()};

    const double A_LOCAL{a_ptr()};
    const double B_LOCAL{b_ptr()};

    for (std::size_t i = 0; i < num_particles; ++i) {
        for (std::size_t j = i + 1; j < num_particles; ++j) {
            double displ_x{}, displ_y{}, displ_z{};

            pbc.displacement(POS_X[i], POS_Y[i], POS_Z[i], POS_X[j], POS_Y[j], POS_Z[j], displ_x, displ_y, displ_z);
            const double DIST_SQ{displ_x * displ_x + displ_y * displ_y + displ_z * displ_z};
            const double DIST{std::sqrt(DIST_SQ)};

            // MASK to get around if statement:
            const bool DEGENERATE{DIST < 1e-12};
            const double INV_DIST{DEGENERATE ? 1.0 : 1.0 / DIST};
            const double MASK{DEGENERATE ? 0.0 : 1.0};

            // u(r) = a*r / (1 + b*r)
            // u'(r) = a / (1 + b*r)^2
            // u''(r) = -2ab / (1 + b*r)^3
            const double DENOM{1.0 + B_LOCAL * DIST};
            const double DENOM_SQ{DENOM * DENOM};
            const double DENOM_CB{DENOM_SQ * DENOM};

            const double FIRST_DERIV{A_LOCAL / DENOM_SQ};
            const double SECOND_DERIV{-2.0 * A_LOCAL * B_LOCAL / DENOM_CB};

            // ∇_i u(r_ij) = u'(r) * (r_vec / r)
            const double GRAD_FACTOR{FIRST_DERIV * INV_DIST};
            const double GRAD_CONTRIBUTION_X{GRAD_FACTOR * displ_x};
            const double GRAD_CONTRIBUTION_Y{GRAD_FACTOR * displ_y};
            const double GRAD_CONTRIBUTION_Z{GRAD_FACTOR * displ_z};

            grad_x[i] += MASK * GRAD_CONTRIBUTION_X;
            grad_y[i] += MASK * GRAD_CONTRIBUTION_Y;
            grad_z[i] += MASK * GRAD_CONTRIBUTION_Z;

            grad_x[j] -= MASK * GRAD_CONTRIBUTION_X;
            grad_y[j] -= MASK * GRAD_CONTRIBUTION_Y;
            grad_z[j] -= MASK * GRAD_CONTRIBUTION_Z;

            // ∇^2 u(r) = u''(r) + (2/r) u'(r)
            const double LAPLACIAN_PAIR{SECOND_DERIV + 2.0 * FIRST_DERIV * INV_DIST};

            laplacian[i] += MASK * LAPLACIAN_PAIR;
            laplacian[j] += MASK * LAPLACIAN_PAIR;
        }
    }
}
