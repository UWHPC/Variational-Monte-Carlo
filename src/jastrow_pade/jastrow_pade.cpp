#include "jastrow_pade.hpp"

#include <cmath>
#include <cstddef>

double JastrowPade::value(
    const Particles& particles,
    const PeriodicBoundaryCondition& pbc
) const noexcept {
    const std::size_t N{particles.numParticles()};

    const double* RESTRICT p_x{particles.posX()};
    const double* RESTRICT p_y{particles.posY()};
    const double* RESTRICT p_z{particles.posZ()};

    const double a_local{a()};
    const double b_local{b()};

    double jastrowPade{};

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i + 1; j < N; ++j) {
            double dx{}, dy{}, dz{};

            pbc.displacement(
                p_x[i], p_y[i], p_z[i],
                p_x[j], p_y[j], p_z[j],
                dx, dy, dz
            );
            const double r_sq{dx * dx + dy * dy + dz * dz};
            const double r{std::sqrt(r_sq)};

            // Mask to get around if statement
            const bool degenerate{r < 1e-12};
            const double mask{degenerate ? 0.0 : 1.0};

            // u(r) = a*r / (1 + b*r)
            const double denom{1.0 + b_local * r};
            jastrowPade += mask * (a_local * r) / denom;
        }
    }
    return jastrowPade;
}

void JastrowPade::addDerivatives(
    const Particles& particles,
    const PeriodicBoundaryCondition& pbc,
    double* RESTRICT gradX,
    double* RESTRICT gradY,
    double* RESTRICT gradZ,
    double* RESTRICT lap
) const noexcept {
    // NOTE: assumes gradX/gradY/gradZ/lap are zero-initialized by caller
    const std::size_t N{particles.numParticles()};

    const double* RESTRICT p_x{particles.posX()};
    const double* RESTRICT p_y{particles.posY()};
    const double* RESTRICT p_z{particles.posZ()};

    const double a_local{a()};
    const double b_local{b()};

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i + 1; j < N; ++j) {
            double dx{}, dy{}, dz{};

            pbc.displacement(
                p_x[i], p_y[i], p_z[i],
                p_x[j], p_y[j], p_z[j],
                dx, dy, dz
            );
            const double r_sq{dx * dx + dy * dy + dz * dz};
            const double r{std::sqrt(r_sq)};

            // Mask to get around if statement:
            const bool degenerate{r < 1e-12};
            const double inv_r{degenerate ? 1.0 : 1 / r};
            const double mask{degenerate ? 0.0 : 1.0};

            // u(r) = a*r / (1 + b*r)
            // u'(r) = a / (1 + b*r)^2
            // u''(r) = -2ab / (1 + b*r)^3
            const double denom{1.0 + b_local * r};
            const double denom_sq{denom * denom};
            const double denom_cb{denom_sq * denom};

            const double uprime{a_local / denom_sq};
            const double usecond{-2.0 * a_local * b_local / denom_cb};

            // ∇_i u(r_ij) = u'(r) * (r_vec / r)
            const double fx{uprime * dx * inv_r};
            const double fy{uprime * dy * inv_r};
            const double fz{uprime * dz * inv_r};

            gradX[i] += mask * fx;
            gradY[i] += mask * fy;
            gradZ[i] += mask * fz;

            gradX[j] -= mask * fx;
            gradY[j] -= mask * fy;
            gradZ[j] -= mask * fz;

            // ∇^2 u(r) = u''(r) + (2/r) u'(r)
            const double lap_pair{usecond + 2.0 * uprime * inv_r};

            lap[i] += mask * lap_pair;
            lap[j] += mask * lap_pair;
        }
    }
}
