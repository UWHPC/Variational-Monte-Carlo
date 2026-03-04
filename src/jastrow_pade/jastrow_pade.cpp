#include "jastrow_pade.hpp"

#include <cmath>
#include <cstddef>

double JastrowPade::value(const Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept {
    const std::size_t numParticles{particles.numParticles()};

    const double* RESTRICT positionX{particles.posX()};
    const double* RESTRICT positionY{particles.posY()};
    const double* RESTRICT positionZ{particles.posZ()};

    const double aParameter{a()};
    const double bParameter{b()};

    double jastrowPade{0.0};

    for (std::size_t particleI = 0; particleI < numParticles; ++particleI) {
        for (std::size_t particleJ = particleI + 1; particleJ < numParticles; ++particleJ) {
            double displacementX{};
            double displacementY{};
            double displacementZ{};

            pbc.displacement(positionX[particleI], positionY[particleI], positionZ[particleI], 
                             positionX[particleJ], positionY[particleJ], positionZ[particleJ], 
                            displacementX, displacementY, displacementZ);
            const double distanceSquared{
                displacementX * displacementX +
                displacementY * displacementY +
                displacementZ * displacementZ
            };
            const double distance{std::sqrt(distanceSquared)};

            // Mask to get around if statement
            const bool particlesCoincident{distance < 1e-12};
            const double mask{particlesCoincident ? 0.0 : 1.0};

            // u(r) = a*r / (1 + b*r)
            const double denominator{1.0 + bParameter * distance};
            jastrowPade += mask * (aParameter * distance) / denominator;
        }
    }
    return jastrowPade;
}

void JastrowPade::addDerivatives(const Particles& particles, const PeriodicBoundaryCondition& pbc,
                                 double* RESTRICT gradX, double* RESTRICT gradY, double* RESTRICT gradZ,
                                 double* RESTRICT lap) const noexcept {
    // NOTE: assumes gradX/gradY/gradZ/lap are zero-initialized by caller
    const std::size_t numParticles{particles.numParticles()};

    const double* RESTRICT positionX{particles.posX()};
    const double* RESTRICT positionY{particles.posY()};
    const double* RESTRICT positionZ{particles.posZ()};

    const double a_local{a()};
    const double b_local{b()};

    for (std::size_t i = 0; i < numParticles; ++i) {
        for (std::size_t j = i + 1; j < numParticles; ++j) {
            double displacementX{};
            double displacementY{};
            double displacementZ{};

            pbc.displacement(positionX[i], positionY[i], positionZ[i], 
                             positionX[j], positionY[j], positionZ[j], 
                            displacementX, displacementY, displacementZ);
            const double distanceSquared{displacementX * displacementX + 
                                         displacementY * displacementY + 
                                         displacementZ * displacementZ};
            const double distance{std::sqrt(distanceSquared)};

            // Mask to get around if statement:
            const bool particlesCoincident{distance < 1e-12};
            const double inv_r{particlesCoincident ? 1.0 : 1 / distance};
            const double mask{particlesCoincident? 0.0 : 1.0};

            // u(r) = a*r / (1 + b*r)
            // u'(r) = a / (1 + b*r)^2
            // u''(r) = -2ab / (1 + b*r)^3
            const double denom{1.0 + b_local * distance};
            const double denom_sq{denom * denom};
            const double denom_cb{denom_sq * denom};

            const double uprime{a_local / denom_sq};
            const double usecond{-2.0 * a_local * b_local / denom_cb};

            // ∇_i u(r_ij) = u'(r) * (r_vec / r)
            const double fx{uprime * displacementX * inv_r};
            const double fy{uprime * displacementY * inv_r};
            const double fz{uprime * displacementZ * inv_r};

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
