#pragma once

#include <cstdlib>
#include <memory>

class Particles {
private:
    std::size_t static constexpr numVectorComponents{3}; // Number of 

    std::size_t numParticles_;              // Number of particles 
    std::unique_ptr<double[]> memoryBlock_; // Memory block size

public:
    explicit Particles(std::size_t numParticles)
    : numParticles_{numParticles} {
        memoryBlock_ = std::make_unique<double[]>(numVectorComponents*numParticles);
    }

    // Number of Particles:
    [[nodiscard]] std::size_t N() const { return numParticles_; }

    // Raw Pointers - Mutable:
    double* posX() { return memoryBlock_.get(); }
    double* posY() { return memoryBlock_.get() + N(); }
    double* posZ() { return memoryBlock_.get() + 2*N(); }

    // Raw Pointers - Immutable:
    double const* posX() const { return memoryBlock_.get(); }
    double const* posY() const { return memoryBlock_.get() + N(); }
    double const* posZ() const { return memoryBlock_.get() + 2*N(); }
};