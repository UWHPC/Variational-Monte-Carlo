#pragma once

#include "../utilities/aligned_soa.hpp"
#include "../utilities/memory.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>

class Particles {
private:
    // All sub-arrays owned by Particles:
    enum ArrayIndex : std::size_t {
        POS_X_,
        POS_Y_,
        POS_Z_,
        GRAD_X_,
        GRAD_Y_,
        GRAD_Z_,
        LOG_PSI_,
        LAP_LOG_PSI_,
        NUM_SUB_ARRAYS_
    };
    // Number of particles:
    std::size_t numParticles_;

    // Aligned memory block:
    AlignedSoA<double> particleData_;

public:
    explicit Particles(std::size_t numParticles)
        : particleData_{numParticles, NUM_SUB_ARRAYS_}, numParticles_{numParticles} {}

    // Physical number of particles
    [[nodiscard]] std::size_t numParticles() const { return numParticles_; }

    // Length of the padded stride
    [[nodiscard]] std::size_t paddingStride() const { return particleData_.stride(); }

    // Raw Pointers:
    // X position of particle
    [[nodiscard]] double* posX() { return particleData_[POS_X_]; }
    [[nodiscard]] double const* posX() const { return particleData_[POS_X_]; }

    // Y position of particle
    [[nodiscard]] double* posY() { return particleData_[POS_Y_]; }
    [[nodiscard]] double const* posY() const { return particleData_[POS_Y_]; }

    // Z position of particle
    [[nodiscard]] double* posZ() { return particleData_[POS_Z_]; }
    [[nodiscard]] double const* posZ() const { return particleData_[POS_Z_]; }

    // X component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiX() { return particleData_[GRAD_X_]; }
    [[nodiscard]] double const* gradLogPsiX() const { return particleData_[GRAD_X_]; }

    // Y component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiY() { return particleData_[GRAD_Y_]; }
    [[nodiscard]] double const* gradLogPsiY() const { return particleData_[GRAD_Y_]; }

    // Z component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiZ() { return particleData_[GRAD_Z_]; }
    [[nodiscard]] double const* gradLogPsiZ() const { return particleData_[GRAD_Z_]; }

    // Log|PSI|
    [[nodiscard]] double* logPsi() { return particleData_[LOG_PSI_]; }
    [[nodiscard]] double const* logPsi() const { return particleData_[LOG_PSI_]; }

    // Laplacian of Log|PSI|
    [[nodiscard]] double* laplogPsi() { return particleData_[LAP_LOG_PSI_]; }
    [[nodiscard]] double const* laplogPsi() const { return particleData_[LAP_LOG_PSI_]; }
};
