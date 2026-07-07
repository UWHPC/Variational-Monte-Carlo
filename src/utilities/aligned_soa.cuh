#pragma once

#include "macros.cuh"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <new>

inline void* backend_alloc(std::size_t alignment, std::size_t size) {
#ifdef VMC_CUDA_BACKEND
  static_cast<void>(alignment);

  void* ptr{};
  CUDA_CHECK(cudaMallocManaged(&ptr, size));

  return ptr;
#elif defined(_MSC_VER)
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

inline void backend_free(void* ptr) {
#ifdef VMC_CUDA_BACKEND
  cudaFree(ptr);
#elif defined(_MSC_VER)
  _aligned_free(ptr);
#else
  std::free(ptr);
#endif
}

struct SoADeleter {
  template <typename T>
  void operator()(T* ptr) const {
    backend_free(ptr);
  }
};

template <arithmetic T>
class AlignedSoA {
private:
  static constexpr std::size_t elements_per_align_{SIMD_BYTES / sizeof(T)};
  std::size_t num_elements_;
  std::size_t stride_length_;
  std::unique_ptr<T[], SoADeleter> data_;

public:
  AlignedSoA() : num_elements_{}, stride_length_{}, data_{} {}
  AlignedSoA(AlignedSoA&&) noexcept = default;
  AlignedSoA& operator=(AlignedSoA&&) noexcept = default;

  AlignedSoA(std::size_t num_elements, std::size_t num_arrays)
  : num_elements_{num_elements}
  , stride_length_{round_up(num_elements)} {
    #ifdef VMC_CUDA_BACKEND
    stride_length_ = num_elements;
    #endif
    
    const std::size_t total_elements{num_arrays * stride_length_};
    const std::size_t total_bytes{total_elements * sizeof(T)};

    T* ptr{static_cast<T*>(backend_alloc(SIMD_BYTES, total_bytes))};
    if (!ptr) { throw std::bad_alloc(); }

    // For once cudaMalloc is used:
    // #ifdef VMC_CUDA_BACKEND
    //   CUDA_CHECK(cudaMemset(ptr, 0, total_bytes));
    // #else
    //   std::fill_n(ptr, total_elements, T{});
    // #endif

    std::fill_n(ptr, total_elements, T{});
    data_.reset(ptr);
  }

  [[nodiscard]] std::size_t num_elements() const {
    return num_elements_;
  }

  [[nodiscard]]
  std::size_t stride() const {
    return stride_length_;
  }

  [[nodiscard]]
  T* operator[](std::size_t array_index) {
    return data_.get() + array_index * stride();
  }

  [[nodiscard]]
  const T* operator[](std::size_t array_index) const {
    return data_.get() + array_index * stride();
  }

  [[nodiscard]]
  static constexpr std::size_t round_up(std::size_t unpadded) {
    return (unpadded + elements_per_align_ - 1) & 
          ~(elements_per_align_ - 1);
  }
};
