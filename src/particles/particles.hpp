#pragma once

#include "../memory/memory.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>

class Particles {
private:
    // Number of components
    static constexpr std::size_t numVectorComponents_{8};

    // SIMD byte alignment
    static constexpr std::size_t alignmentBytes_{SIMD_BYTES};

    // Ensures sub-arrays are byte aligned
    static constexpr std::size_t doublesPerAlignment_{SIMD_BYTES / sizeof(double)};

    // Number of particles
    std::size_t numParticles_;

    // Ensures all sub-arrays are aligned
    std::size_t alignmentPadding_;

    // Memory block size
    std::unique_ptr<double[], AlignedDeleter> memoryBlock_;

public:
    explicit Particles(std::size_t numParticles)
        : numParticles_{numParticles},
          alignmentPadding_{(numParticles_ + doublesPerAlignment_ - 1) & ~(doublesPerAlignment_ - 1)} {
        // Round to nearest multiple of SIMD width - needed to ensure sub-arrays are byte aligned:

        // Determine size of the memory block in bytes:
        std::size_t const memoryBlockSize{numVectorComponents_ * alignmentPadding_};
        std::size_t const blockSizeBytes{memoryBlockSize * sizeof(double)};

        // Proper byte alignment and allocation:
        double* ptr{static_cast<double*>(alignedAlloc(alignmentBytes_, blockSizeBytes))};
        if (!ptr)
            throw std::bad_alloc();

        // Zero initialization:
        std::fill_n(ptr, memoryBlockSize, 0.0);

        // Pointer owner transfership:
        memoryBlock_.reset(ptr);
    }

    // Logical number of particles
    [[nodiscard]] std::size_t numParticles() const { return numParticles_; }

    // Memory stride for flat arrays
    [[nodiscard]] std::size_t paddingStride() const { return alignmentPadding_; }

    // Raw Pointers - Mutable:
    // MUT - X position of particle
    [[nodiscard]] double* posX() { return memoryBlock_.get(); }

    // MUT - Y position of particle
    [[nodiscard]] double* posY() { return memoryBlock_.get() + paddingStride(); }

    // MUT - Z position of particle
    [[nodiscard]] double* posZ() { return memoryBlock_.get() + 2 * paddingStride(); }

    // MUT - X component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiX() { return memoryBlock_.get() + 3 * paddingStride(); }

    // MUT - Y component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiY() { return memoryBlock_.get() + 4 * paddingStride(); }

    // MUT - Z component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiZ() { return memoryBlock_.get() + 5 * paddingStride(); }

    // MUT - Log|PSI|
    [[nodiscard]] double* logPsi() { return memoryBlock_.get() + 6 * paddingStride(); }

    // MUT - Laplacian of Log|PSI|
    [[nodiscard]] double* laplogPsi() { return memoryBlock_.get() + 7 * paddingStride(); }

    // Raw Pointers - Immutable:
    // IMMUT - X position of particle
    [[nodiscard]] double const* posX() const { return memoryBlock_.get(); }

    // IMMUT - Y position of particle
    [[nodiscard]] double const* posY() const { return memoryBlock_.get() + paddingStride(); }

    // IMMUT - Z position of particle
    [[nodiscard]] double const* posZ() const { return memoryBlock_.get() + 2 * paddingStride(); }

    // IMMUT - X component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiX() const { return memoryBlock_.get() + 3 * paddingStride(); }

    // IMMUT - Y component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiY() const { return memoryBlock_.get() + 4 * paddingStride(); }

    // IMMUT - Z component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiZ() const { return memoryBlock_.get() + 5 * paddingStride(); }

    // IMMUT - Log|PSI|
    [[nodiscard]] double const* logPsi() const { return memoryBlock_.get() + 6 * paddingStride(); }

    // IMMUT - Laplacian of Log|PSI|
    [[nodiscard]] double const* laplogPsi() const { return memoryBlock_.get() + 7 * paddingStride(); }
};
