#pragma once

#include "memory.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>

template <typename T> class AlignedSoA {
private:
    // SIMD byte alignment
    static constexpr std::size_t ALIGNMENT_BYTES{SIMD_BYTES};

    // Ensures sub-arrays are byte aligned
    static constexpr std::size_t ELEMENTS_PER_ALIGNMENTS{SIMD_BYTES / sizeof(T)};

    std::size_t num_elements_;
    std::size_t stride_length_;
    std::size_t num_arrays_;
    std::unique_ptr<T[], AlignedDeleter> memory_block_;

    // Round up to nearest factor of SIMD bytes:
    std::size_t round_up(std::size_t unpadded) const {
        return (unpadded + ELEMENTS_PER_ALIGNMENTS - 1) & ~(ELEMENTS_PER_ALIGNMENTS - 1);
    }

public:
    AlignedSoA(std::size_t num_elements, std::size_t num_arrays)
        : num_elements_{num_elements}, stride_length_{round_up(num_elements)}, num_arrays_{num_arrays} {
        // Determine the total elements and bytes needed for padded memory block:
        std::size_t total_elements{num_arrays_ * stride_length_};
        std::size_t total_bytes{total_elements * sizeof(T)};

        // Allocate aligned bytes and check:
        T* ptr{static_cast<T*>(alignedAlloc(ALIGNMENT_BYTES, total_bytes))};
        if (!ptr)
            throw std::bad_alloc();

        // Initialize the pointer to all default and transfer ownership to memory_block_:
        std::fill_n(ptr, total_elements, T{});
        memory_block_.reset(ptr);
    }

    // Getters:
    // Padded stride length:
    [[nodiscard]] std::size_t stride() const { return stride_length_; }

    // Number of elements:
    [[nodiscard]] std::size_t num_elements() const { return num_elements_; }

    // Raw pointer accessors:
    // Mutable:
    T* operator[](std::size_t array_index) { return memory_block_.get() + array_index * stride(); }

    // Immutable:
    const T* operator[](std::size_t array_index) const { return memory_block_.get() + array_index * stride(); }
};
