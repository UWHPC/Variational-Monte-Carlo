#pragma once

#include <memory>
#include <cstdlib>

// Aligned allocation:
#if defined(_WIN32)
    #include <malloc.h>
    inline void* alignedAlloc(std::size_t alignment, std::size_t size) {
        return _aligned_malloc(size, alignment);
    }
    inline void alignedFree(void* ptr) {
        _aligned_free(ptr);
    }
#elif defined(__APPLE__)
    inline void* alignedAlloc(std::size_t alignment, std::size_t size) {
        void* ptr = nullptr;
        posix_memalign(&ptr, alignment, size);
        return ptr;
    }
    inline void alignedFree(void* ptr) {
        std::free(ptr);
    }
#else
    inline void* alignedAlloc(std::size_t alignment, std::size_t size) {
        return std::aligned_alloc(alignment, size);
    }
    inline void alignedFree(void* ptr) {
        std::free(ptr);
    }
#endif

// Deletion of memory for aligned_alloc:
struct AlignedDeleter {
    void operator()(double* ptr) const { alignedFree(ptr); }
};

// Detect SIMD width:
#if defined(__AVX512F__)
    constexpr std::size_t SIMD_BYTES{64};
#elif defined(__AVX2__)
    constexpr std::size_t SIMD_BYTES{32};
#elif defined(__AVX__)
    constexpr std::size_t SIMD_BYTES{32};
#elif defined(__SSE2__)
    constexpr std::size_t SIMD_BYTES{16};
#else
    constexpr std::size_t SIMD_BYTES{sizeof(double)};
#endif