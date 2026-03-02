#pragma once

#include "../memory/memory.hpp"

#include <cstdlib>
#include <memory>
#include <cstring>
#include <stdexcept>

class Particles {
private:
    static constexpr std::size_t numVectorComponents_{8};   // Number of components
    static constexpr std::size_t alignmentBytes_{64};       // 64 byte alignment
    static constexpr std::size_t doublesPerAlignment_{8};   // Ensures sub-arrays are 64 byte aligned

    std::size_t numParticles_;                              // Number of particles 
    std::size_t alignmentPadding_;
    std::unique_ptr<double[], AlignedDeleter> memoryBlock_; // Memory block size

public:
    explicit Particles(std::size_t numParticles)
    : numParticles_{numParticles}
    , alignmentPadding_{(numParticles_ + doublesPerAlignment_ - 1) & ~(doublesPerAlignment_ - 1)} {
        // Round to nearest multiple of 8 - needed to ensure sub-arrays are 64 byte aligned:

        // Determine size of the memory block in bytes:
        std::size_t const memoryBlockSize{numVectorComponents_*alignmentPadding_};
        std::size_t const blockSizeBytes{memoryBlockSize*sizeof(double)};
        
        // 64 byte alignment and allocation:
        double* ptr = static_cast<double*>(alignedAlloc(alignmentBytes_, blockSizeBytes));
        if (!ptr) { throw std::bad_alloc(); }

        // Zero initialization:
        std::fill_n(ptr, memoryBlockSize, 0.0);

        // Pointer owner transfership:
        memoryBlock_.reset(ptr);
    }

    // Number of Particles:
    [[nodiscard]] std::size_t numParticles() const { return numParticles_; }                           // Logical number of particles
    [[nodiscard]] std::size_t paddingStride() const { return alignmentPadding_; }                      // Memory stride for flat ararys

    // Raw Pointers - Mutable:
    [[nodiscard]] double* posX() { return memoryBlock_.get(); }                                        // MUT - X position of particle
    [[nodiscard]] double* posY() { return memoryBlock_.get() + paddingStride(); }                      // MUT - Y position of particle
    [[nodiscard]] double* posZ() { return memoryBlock_.get() + 2*paddingStride(); }                    // MUT - Z position of particle

    [[nodiscard]] double* gradLogPsiX() { return memoryBlock_.get() + 3*paddingStride(); }             // MUT - X component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiY() { return memoryBlock_.get() + 4*paddingStride(); }             // MUT - Y component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiZ() { return memoryBlock_.get() + 5*paddingStride(); }             // MUT - Z component of gradient( log|PSI| )

    [[nodiscard]] double* logPsi() { return memoryBlock_.get() + 6*paddingStride(); }                  // MUT - Log|PSI|
    [[nodiscard]] double* laplogPsi() { return memoryBlock_.get() + 7*paddingStride(); }               // MUT - Laplacian of Log|PSI|

    // Raw Pointers - Immutable:
    [[nodiscard]] double const* posX() const { return memoryBlock_.get(); }                            // IMMUT - X position of particle
    [[nodiscard]] double const* posY() const { return memoryBlock_.get() + paddingStride(); }          // IMMUT - Y position of particle
    [[nodiscard]] double const* posZ() const { return memoryBlock_.get() + 2*paddingStride(); }        // IMMUT - Z position of particle

    [[nodiscard]] double const* gradLogPsiX() const { return memoryBlock_.get() + 3*paddingStride(); } // IMMUT - X component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiY() const { return memoryBlock_.get() + 4*paddingStride(); } // IMMUT - Y component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiZ() const { return memoryBlock_.get() + 5*paddingStride(); } // IMMUT - Z component of gradient( log|PSI| )

    [[nodiscard]] double const* logPsi() const { return memoryBlock_.get() + 6*paddingStride(); }      // IMMUT - Log|PSI|
    [[nodiscard]] double const* laplogPsi() const { return memoryBlock_.get() + 7*paddingStride(); }   // IMMUT - Laplacian of Log|PSI|
};