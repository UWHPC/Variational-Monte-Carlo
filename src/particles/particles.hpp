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
    static constexpr std::size_t POS_X{0}, POS_Y{1}, POS_Z{2};
    static constexpr std::size_t GRAD_X{3}, GRAD_Y{4}, GRAD_Z{5};
    static constexpr std::size_t LOG_PSI{6}, LAP_LOG_PSI{7};

    // Number of sub-arrays:
    static constexpr std::size_t NUM_SUB_ARRAYS{8};

    // Number of particles:
    std::size_t numParticles_;

    // Aligned memory block:
    AlignedSoA particleData_;

public:
    explicit Particles(std::size_t numParticles)
        : particleData_{numParticles, NUM_SUB_ARRAYS}, numParticles_{numParticles} {}

    // Physical number of particles
    [[nodiscard]] std::size_t numParticles() const { return numParticles_; }

    // Length of the padded stride
    [[nodiscard]] std::size_t paddingStride() const { return particleData_.stride(); }

    // Raw Pointers - Mutable:
    // MUT - X position of particle
    [[nodiscard]] double* posX() { return particleData_[POS_X]; }

    // MUT - Y position of particle
    [[nodiscard]] double* posY() { return particleData_[POS_Y]; }

    // MUT - Z position of particle
    [[nodiscard]] double* posZ() { return particleData_[POS_Z]; }

    // MUT - X component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiX() { return particleData_[GRAD_X]; }

    // MUT - Y component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiY() { return particleData_[GRAD_Y]; }

    // MUT - Z component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiZ() { return particleData_[GRAD_Z]; }

    // MUT - Log|PSI|
    [[nodiscard]] double* logPsi() { return particleData_[LOG_PSI]; }

    // MUT - Laplacian of Log|PSI|
    [[nodiscard]] double* laplogPsi() { return particleData_[LAP_LOG_PSI]; }

    // Raw Pointers - Immutable:
    // IMMUT - X position of particle
    [[nodiscard]] double const* posX() const { return particleData_[POS_X]; }

    // IMMUT - Y position of particle
    [[nodiscard]] double const* posY() const { return particleData_[POS_Y]; }

    // IMMUT - Z position of particle
    [[nodiscard]] double const* posZ() const { return particleData_[POS_Z]; }

    // IMMUT - X component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiX() const { return particleData_[GRAD_X]; }

    // IMMUT - Y component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiY() const { return particleData_[GRAD_Y]; }

    // IMMUT - Z component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiZ() const { return particleData_[GRAD_Z]; }

    // IMMUT - Log|PSI|
    [[nodiscard]] double const* logPsi() const { return particleData_[LOG_PSI]; }

    // IMMUT - Laplacian of Log|PSI|
    [[nodiscard]] double const* laplogPsi() const { return particleData_[LAP_LOG_PSI]; }
};
