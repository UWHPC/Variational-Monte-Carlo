#pragma once

#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"

#include <cstddef>
#include <vector>

class SlaterPlaneWave {
private:
    std::size_t N_{};
    double L_{};
    // to do: cached std::vector matrices
    // to do: k-vectors
public:
    [[nodiscard]] double logAbsDet(const Particles& particles, const PeriodicBoundaryCondition& pbc);
    
    void addDerivatives(
        const Particles& particles,
        const PeriodicBoundaryCondition& pbc,
        double* gradX,
        double* gradY,
        double* gradZ,
        double* la
    ) const noexcept;
};