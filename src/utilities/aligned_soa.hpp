#pragma once

#include "memory.hpp"
#include <cstddef>
#include <memory>

template <typename T> class AlignedSoA {
private:
    // SIMD byte alignment
    static constexpr std::size_t alignmentBytes{SIMD_BYTES};

    // Ensures sub-arrays are byte aligned
    static constexpr std::size_t doublesPerAlignment{SIMD_BYTES / sizeof(double)};

    std::size_t numElements_;
    std::size_t strideLength_;
    std::size_t numArrays_;
    std::unique_ptr<T[], AlignedDeleter> memoryBlock_;

    // Round up to nearest factor of SIMD bytes:
    std::size_t roundUp(std::size_t unpadded) const {
        return (unpadded + doublesPerAlignment - 1) & ~(doublesPerAlignment - 1);
    }

public:
    AlignedSoA(std::size_t numElements, std::size_t numArrays)
        : numElements_{numElements}, strideLength_{roundUp(numElements)}, numArrays_{numArrays} {
        // Determine the total doubles and bytes needed for padded memory block:
        std::size_t totalDoubles{numArrays_ * strideLength_};
        std::size_t totalBytes{totalDoubles * sizeof(double)};

        // Allocate aligned bytes and check:
        double* ptr{static_cast<double*>(alignedAlloc(alignmentBytes, totalBytes))};
        if (!ptr)
            throw std::bad_alloc();

        // Initialize the pointer to all 0.0 and transfer ownership to memoryBlock_:
        std::fill_n(ptr, totalDoubles, 0.0);
        memoryBlock_.reset(ptr);
    }

    // Getters:
    // Padded stride length:
    [[nodiscard]] std::size_t stride() const { return strideLength_; }

    // Number of elements:
    [[nodiscard]] std::size_t numElements() const { return numElements_; }

    // Raw pointer accessors:
    // Mutable:
    double* operator[](std::size_t arrayIndex) { return memoryBlock_.get() + arrayIndex * stride(); }

    // Immutable:
    const double* operator[](std::size_t arrayIndex) const { return memoryBlock_.get() + arrayIndex * stride(); }
};