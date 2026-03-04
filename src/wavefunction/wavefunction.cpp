#include "wavefunction.hpp"

#include <algorithm>

void WaveFunction::evaluateLogPsi(
    Particles& particles,
    const PeriodicBoundaryCondition& pbc
    ) {
    const double logDet{slaterPlaneWave_.logAbsDet(particles, pbc)};
    const double jastrowPade{jastrowPade_.value(particles, pbc)};
    particles.logPsi()[0] = logDet + jastrowPade;
}

void WaveFunction::evaluateDerivatives(Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept {
    const std::size_t paddedStride{particles.paddingStride()};

    std::fill_n(particles.gradLogPsiX(), paddedStride, 0.0);
    std::fill_n(particles.gradLogPsiY(), paddedStride, 0.0);
    std::fill_n(particles.gradLogPsiZ(), paddedStride, 0.0);
    std::fill_n(particles.laplogPsi(), paddedStride, 0.0);

    slaterPlaneWave_.addDerivatives(particles, pbc, particles.gradLogPsiX(), particles.gradLogPsiY(), 
                                    particles.gradLogPsiZ(), particles.laplogPsi());

    jastrowPade_.addDerivatives(particles, pbc, particles.gradLogPsiX(), particles.gradLogPsiY(),
                                particles.gradLogPsiZ(), particles.laplogPsi());
}