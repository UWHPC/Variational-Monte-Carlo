#pragma once

#include "macros.cuh"

template <typename T> struct [[nodiscard]] Ptr3D {
  T* RESTRICT x_;
  T* RESTRICT y_;
  T* RESTRICT z_;

  [[nodiscard]]
  Ptr3D align() const noexcept {
    Ptr3D p{x_, y_, z_};

    ASSUME_ALIGNED(p.x_, SIMD_BYTES);
    ASSUME_ALIGNED(p.y_, SIMD_BYTES);
    ASSUME_ALIGNED(p.z_, SIMD_BYTES);

    return p;
  }
};