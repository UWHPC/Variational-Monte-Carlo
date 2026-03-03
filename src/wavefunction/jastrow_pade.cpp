#include "jastrow_pade.hpp"

#include <cmath>
#include <cstddef>

double JastrowPade::value(const Particles& particles, const periodicBoundaryCondition& pbc) const noexcept {
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
            const double r_sq{dx*dx + dy*dy + dz*dz};
            if (r_sq < 1e-24) continue; // defensive: avoid pathological coincident particles

            const double r{std::sqrt(r_sq)};

            // u(r) = a*r / (1 + b*r)
            const double denom{1.0 + b_local*r};
            jastrowPade += (a_local * r) / denom;
        }
    }

    return jastrowPade;
}

void JastrowPade::addDerivatives(
    const Particles& particles,
    const periodicBoundaryCondition& pbc,
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
        for (std::size_t j = i + 1; j < N; ++j ) {
            double dx{}, dy{}, dz{};
            
            pbc.displacement(
                p_x[i], p_y[i], p_z[i], 
                p_x[j], p_y[j], p_z[j], 
                dx, dy, dz
            );
            const double r_sq{dx*dx + dy*dy + dz*dz};
            if (r_sq < 1e-24) continue; // 1/r will blow up if r2 < 1e-24.

            const double r{std::sqrt(r_sq)};
            const double inv_r{1.0/r};

            // u(r) = a*r / (1 + b*r)
            // u'(r) = a / (1 + b*r)^2
            // u''(r) = -2ab / (1 + b*r)^3
            const double denom{1.0 + b_local*r};
            const double denom_sq{denom * denom};
            const double denom_cb{denom_sq * denom};

            const double uprime{a_local / denom_sq};
            const double usecond{-2.0*a_local*b_local / denom_cb};

            // ∇_i u(r_ij) = u'(r) * (r_vec / r)
            const double fx{uprime * dx * inv_r};
            const double fy{uprime * dy * inv_r};
            const double fz{uprime * dz * inv_r};

            gradX[i] += fx;
            gradY[i] += fy;
            gradZ[i] += fz;

            gradX[j] -= fx;
            gradY[j] -= fy;
            gradZ[j] -= fz;

            // ∇^2 u(r) = u''(r) + (2/r) u'(r)
            const double lap_pair{usecond + 2.0*uprime*inv_r};
            
            lap[i] += lap_pair;
            lap[j] += lap_pair;
        }
    }
}