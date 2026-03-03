#include "wavefunction.hpp"

#include <algorithm>

// to do!!!
// void WaveFunction::evaluateLogPsi(
//     Particles& particles, 
//     const PeriodicBoundaryCondition& pbc
//     ) const noexcept {
//     particles.logPsi()[0] = jastrowPade_.value(particles, pbc);
// }

void WaveFunction::evaluateDerivatives(Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept {
    const std::size_t paddedStride{particles.paddingStride()};

    std::fill_n(particles.gradLogPsiX(), paddedStride, 0.0);
    std::fill_n(particles.gradLogPsiY(), paddedStride, 0.0);
    std::fill_n(particles.gradLogPsiZ(), paddedStride, 0.0);
    std::fill_n(particles.laplogPsi(),   paddedStride, 0.0);

    jastrowPade_.addDerivatives(
        particles, 
        pbc, 
        particles.gradLogPsiX(),
        particles.gradLogPsiY(), 
        particles.gradLogPsiZ(), 
        particles.laplogPsi()
    );
}