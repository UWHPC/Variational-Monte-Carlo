#pragma once
#include <cstddef>
#include <vector>
#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"

class SlaterPlaneWave {
private:
    std::size_t N_{};
    double L_{};
    // to do: cached std::vector matrices
    // to do: k-vectors
public:
    [[nodiscard]] double logAbsDet(const Particles<>& particles, const periodicBoundaryCondition& pbc);
    
    void addDerivatives(const Particles<>& p, const periodicBoundaryCondition& pbc, double* gradX, double* gradY, double* gradZ, double* lap) const noexcept;
};