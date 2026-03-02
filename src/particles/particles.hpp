#pragma once

#include "../memory/memory.hpp"

#include <cstdlib>
#include <memory>
#include <cstring>

class Particles {
private:
    static constexpr std::size_t numVectorComponents{8};    // Number of components

    std::size_t numParticles_;                              // Number of particles 
    std::unique_ptr<double[], AlignedDeleter> memoryBlock_; // Memory block size

public:
    explicit Particles(std::size_t numParticles)
    : numParticles_{numParticles} {
        // Determine size of the memory block in bytes:
        std::size_t const memoryBlockSize{numVectorComponents*numParticles};
        std::size_t const blockSizeBytes{memoryBlockSize*sizeof(double)};
        
        // 64 byte alignment and allocation:
        std::size_t const alignment{64};
        double* ptr = static_cast<double*>(alignedAlloc(alignment, blockSizeBytes));

        // Zero initialization:
        std::memset(ptr, 0, blockSizeBytes);

        // Pointer owner transfership:
        memoryBlock_.reset(ptr);
    }

    // Number of Particles:
    [[nodiscard]] std::size_t N() const { return numParticles_; }

    // Raw Pointers - Mutable:
    [[nodiscard]] double* posX() { return memoryBlock_.get(); }                            // MUT - X position of particle
    [[nodiscard]] double* posY() { return memoryBlock_.get() + N(); }                      // MUT - Y position of particle
    [[nodiscard]] double* posZ() { return memoryBlock_.get() + 2*N(); }                    // MUT - Z position of particle

    [[nodiscard]] double* gradLogPsiX() { return memoryBlock_.get() + 3*N(); }             // MUT - X component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiY() { return memoryBlock_.get() + 4*N(); }             // MUT - Y component of gradient( log|PSI| )
    [[nodiscard]] double* gradLogPsiZ() { return memoryBlock_.get() + 5*N(); }             // MUT - Z component of gradient( log|PSI| )

    [[nodiscard]] double* logPsi() { return memoryBlock_.get() + 6*N(); }                  // MUT - Log|PSI|
    [[nodiscard]] double* laplogPsi() { return memoryBlock_.get() + 7*N(); }               // MUT - Laplacian of Log|PSI|

    // Raw Pointers - Immutable:
    [[nodiscard]] double const* posX() const { return memoryBlock_.get(); }                // IMMUT - X position of particle
    [[nodiscard]] double const* posY() const { return memoryBlock_.get() + N(); }          // IMMUT - Y position of particle
    [[nodiscard]] double const* posZ() const { return memoryBlock_.get() + 2*N(); }        // IMMUT - Z position of particle

    [[nodiscard]] double const* gradLogPsiX() const { return memoryBlock_.get() + 3*N(); } // IMMUT - X component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiY() const { return memoryBlock_.get() + 4*N(); } // IMMUT - Y component of gradient( log|PSI| )
    [[nodiscard]] double const* gradLogPsiZ() const { return memoryBlock_.get() + 5*N(); } // IMMUT - Z component of gradient( log|PSI| )

    [[nodiscard]] double const* logPsi() const { return memoryBlock_.get() + 6*N(); }      // IMMUT - Log|PSI|
    [[nodiscard]] double const* laplogPsi() const { return memoryBlock_.get() + 7*N(); }   // IMMUT - Laplacian of Log|PSI|
};