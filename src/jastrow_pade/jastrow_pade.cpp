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
                                 double* RESTRICT gradientX, double* RESTRICT gradientY, double* RESTRICT gradientZ,
                                 double* RESTRICT laplacian) const noexcept {
    // NOTE: assumes gradX/gradY/gradZ/lap are zero-initialized by caller
    const std::size_t numParticles{particles.numParticles()};

    const double* RESTRICT positionX{particles.posX()};
    const double* RESTRICT positionY{particles.posY()};
    const double* RESTRICT positionZ{particles.posZ()};

    const double aParameter{a()};
    const double bParameter{b()};

    for (std::size_t particleI = 0; particleI < numParticles; ++particleI) {
        for (std::size_t particleJ = particleI + 1; particleJ < numParticles; ++particleJ) {
            double displacementX{};
            double displacementY{};
            double displacementZ{};

            pbc.displacement(positionX[particleI], positionY[particleI], positionZ[particleI], 
                             positionX[particleJ], positionY[particleJ], positionZ[particleJ], 
                            displacementX, displacementY, displacementZ);
            const double distanceSquared{displacementX * displacementX + 
                                         displacementY * displacementY + 
                                         displacementZ * displacementZ};
            const double distance{std::sqrt(distanceSquared)};

            // Mask to get around if statement:
            const bool particlesCoincident{distance < 1e-12};
            const double inverseDistance{particlesCoincident ? 1.0 : 1.0 / distance};
            const double mask{particlesCoincident ? 0.0 : 1.0};

            // u(r) = a*r / (1 + b*r)
            // u'(r) = a / (1 + b*r)^2
            // u''(r) = -2ab / (1 + b*r)^3
            const double denominator{1.0 + bParameter * distance};
            const double denominatorSquared{denominator * denominator};
            const double denominatorCubed{denominatorSquared * denominator};

            const double firstDerivative{aParameter / denominatorSquared};
            const double secondDerivative{-2.0 * aParameter * bParameter / denominatorCubed};

            // ∇_i u(r_ij) = u'(r) * (r_vec / r)
            const double gradientFactor{firstDerivative * inverseDistance};
            const double gradientContributionX{gradientFactor * displacementX};
            const double gradientContributionY{gradientFactor * displacementY};
            const double gradientContributionZ{gradientFactor * displacementZ};

            gradientX[particleI] += mask * gradientContributionX;
            gradientY[particleI] += mask * gradientContributionY;
            gradientZ[particleI] += mask * gradientContributionZ;

            gradientX[particleJ] -= mask * gradientContributionX;
            gradientY[particleJ] -= mask * gradientContributionY;
            gradientZ[particleJ] -= mask * gradientContributionZ;

            // ∇^2 u(r) = u''(r) + (2/r) u'(r)
            const double laplacianPair{secondDerivative + 2.0 * firstDerivative * inverseDistance};

            laplacian[particleI] += mask * laplacianPair;
            laplacian[particleJ] += mask * laplacianPair;
        }
    }
}
