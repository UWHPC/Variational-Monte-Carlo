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
    AlignedSoA<double[]> particleData_;

public:
    explicit Particles(std::size_t numParticles)
        : particleData_{numParticles, NUM_SUB_ARRAYS_}, numParticles_{numParticles} {}

    // Physical number of particles
    [[nodiscard]] std::size_t numParticles() const { return numParticles_; }

    // Length of the padded stride
    [[nodiscard]] std::size_t paddingStride() const { return particleData_.stride(); }

    // Raw Pointers - Mutable:
    // MUT - X position of particle
    [[nodiscard]] double* posX() { return particleData_[POS_X_]; }

    // MUT - Y position of particle
    [[nodiscard]] double* posY() { return particleData_[POS_Y_]; }

    // MUT - Z position of particle
    [[nodiscard]] double* posZ() { return particleData_[POS_Z_]; }

    // MUT - X component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiX() { return particleData_[GRAD_X_]; }

    // MUT - Y component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiY() { return particleData_[GRAD_Y_]; }

    // MUT - Z component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiZ() { return particleData_[GRAD_Z_]; }

    // MUT - Log|PSI|
    [[nodiscard]] double* logPsi() { return particleData_[LOG_PSI_]; }

    // MUT - Laplacian of Log|PSI|
    [[nodiscard]] double* laplogPsi() { return particleData_[LAP_LOG_PSI_]; }

    // Raw Pointers - Immutable:
    // IMMUT - X position of particle
    [[nodiscard]] double const* posX() const { return particleData_[POS_X_]; }

    // IMMUT - Y position of particle
    [[nodiscard]] double const* posY() const { return particleData_[POS_Y_]; }

    // IMMUT - Z position of particle
    [[nodiscard]] double const* posZ() const { return particleData_[POS_Z_]; }

    // IMMUT - X component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiX() const { return particleData_[GRAD_X_]; }

    // IMMUT - Y component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiY() const { return particleData_[GRAD_Y_]; }

    // IMMUT - Z component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiZ() const { return particleData_[GRAD_Z_]; }

    // IMMUT - Log|PSI|
    [[nodiscard]] double const* logPsi() const { return particleData_[LOG_PSI_]; }

    // IMMUT - Laplacian of Log|PSI|
    [[nodiscard]] double const* laplogPsi() const { return particleData_[LAP_LOG_PSI_]; }
};
