#include "wavefunction.hpp"
#include <algorithm> // needed for fill_n

void WaveFunction::evaluateLogPsi(Particles<>& particles, const periodicBoundaryCondition& pbc) const noexcept
{
    particles.logPsi()[0] = j_.value(particles, pbc);
}

void WaveFunction::evaluateDerivatives(Particles<>& particles, const periodicBoundaryCondition& pbc) const noexcept
{
    const std::size_t stride = particles.paddingStride();

    std::fill_n(particles.gradLogPsiX(), stride, 0.0);
    std::fill_n(particles.gradLogPsiY(), stride, 0.0);
    std::fill_n(particles.gradLogPsiZ(), stride, 0.0);
    std::fill_n(particles.laplogPsi(),   stride, 0.0);

    j_.addDerivatives(particles, pbc, particles.gradLogPsiX(), particles.gradLogPsiY(), particles.gradLogPsiZ(), particles.laplogPsi());
}