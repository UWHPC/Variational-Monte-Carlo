#pragma once

#include "macros.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <new>

#if defined(_WIN32)
#include <malloc.h>
#endif

inline void* aligned_alloc_backend(std::size_t alignment, std::size_t size) {
#if defined(_WIN32)
  return _aligned_malloc(size, alignment);
#elif defined(__APPLE__) || defined(__GNUC__) || defined(__clang__)
  void* ptr{};
  alignment = std::max(alignment, sizeof(void*));
  if (posix_memalign(&ptr, alignment, size) != 0) {
    return nullptr;
  }
  return ptr;
#else
  return std::aligned_alloc(alignment, size);
#endif
}

inline void aligned_free_backend(void* ptr) {
#if defined(_WIN32)
  _aligned_free(ptr);
#else
  std::free(ptr);
#endif
}

struct AlignedDeleter {
  template <typename T>
  void operator()(T* ptr) const {
    aligned_free_backend(ptr);
  }
};

template <arithmetic T> class AlignedSoA {
private:
  static constexpr std::size_t elements_per_align_{SIMD_BYTES / sizeof(T)};
  std::size_t num_elements_{};
  std::size_t stride_length_{};
  std::unique_ptr<T[], AlignedDeleter> memory_block_;

public:
  AlignedSoA() = default;
  AlignedSoA(AlignedSoA&&) noexcept = default;
  AlignedSoA& operator=(AlignedSoA&&) noexcept = default;

  AlignedSoA(std::size_t num_elements, std::size_t num_arrays)
  : num_elements_{num_elements}
  , stride_length_{round_up(num_elements)}
  {
    const std::size_t total_elements{num_arrays * stride_length_};
    const std::size_t total_bytes{total_elements * sizeof(T)};

    T* ptr{static_cast<T*>(aligned_alloc_backend(SIMD_BYTES, total_bytes))};
    if (!ptr) {
      throw std::bad_alloc();
    }

    std::fill_n(ptr, total_elements, T{});
    memory_block_.reset(ptr);
  }

  [[nodiscard]] std::size_t stride() const { return stride_length_; }
  [[nodiscard]] std::size_t num_elements() const { return num_elements_; }

  [[nodiscard]] T* operator[](std::size_t array_index) {
    return memory_block_.get() + array_index * stride();
  }
  [[nodiscard]] const T* operator[](std::size_t array_index) const {
    return memory_block_.get() + array_index * stride();
  }

  [[nodiscard]] static constexpr std::size_t round_up(std::size_t unpadded) {
    return (unpadded + elements_per_align_ - 1) & ~(elements_per_align_ - 1);
  }
};
