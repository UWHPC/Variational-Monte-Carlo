#pragma once

#include "macros.cuh"

#if defined(VMC_FAST_MATH) && defined(FP_64)
  #error "VMC_FAST_MATH requires single precision: build with FP_64=OFF."
#endif

namespace vmc {

CUDA_CALLABLE
inline void sincos(real_t theta, real_t* s, real_t* c) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  #ifdef FP_64
    ::sincos(theta, s, c);
  #elif defined(VMC_FAST_MATH)
    __sincosf(theta, s, c);
  #else
    sincosf(theta, s, c);
  #endif
#elif defined(__APPLE__)
  #ifdef FP_64
    __sincos(theta, s, c);
  #else
    __sincosf(theta, s, c);
  #endif
#elif (defined(__GNUC__) || defined(__clang__)) && defined(_GNU_SOURCE)
  #ifdef FP_64
    ::sincos(theta, s, c);
  #else
    sincosf(theta, s, c);
  #endif
#else
  *s = std::sin(theta);
  *c = std::cos(theta);
#endif
}

CUDA_CALLABLE [[nodiscard]]
inline real_t sqrt(real_t arg) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  #ifdef FP_64
    return ::sqrt(arg);
  #elif defined(VMC_FAST_MATH)
    return __fsqrt_rn(arg);
  #else
    return sqrtf(arg);
  #endif
#else
  return std::sqrt(arg);
#endif
}

CUDA_CALLABLE [[nodiscard]]
inline real_t exp(real_t arg) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  #ifdef FP_64
    return ::exp(arg);
  #elif defined(VMC_FAST_MATH)
    return __expf(arg);
  #else
    return expf(arg);
  #endif
#else
  return std::exp(arg);
#endif
}

CUDA_CALLABLE [[nodiscard]]
inline real_t pow(real_t arg, real_t power) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  #ifdef FP_64
    return ::pow(arg, power);
  #elif defined(VMC_FAST_MATH)
    return __powf(arg, power);;
  #else
    return powf(arg, power);
  #endif
#else
  return std::pow(arg, power);
#endif
}

CUDA_CALLABLE [[nodiscard]]
inline real_t cbrt(real_t arg) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  #ifdef FP_64
    return ::cbrt(arg);
  #elif defined(VMC_FAST_MATH)
    return copysignf(__powf(fabsf(arg), 1.0f / 3.0f), arg);
  #else
    return cbrtf(arg);
  #endif
#else
  return std::cbrt(arg);
#endif
}

CUDA_CALLABLE [[nodiscard]]
inline real_t log(real_t arg) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  #ifdef FP_64
    return ::log(arg);
  #elif defined(VMC_FAST_MATH)
    return __logf(arg);
  #else
    return logf(arg);
  #endif
#else
  return std::log(arg);
#endif
}

CUDA_CALLABLE [[nodiscard]]
inline real_t abs(real_t arg) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  #ifdef FP_64
    return fabs(arg);
  #else
    return fabsf(arg);
  #endif
#else
  return std::abs(arg);
#endif
}

CUDA_CALLABLE [[nodiscard]]
inline real_t floor(real_t arg) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  #ifdef FP_64
    return ::floor(arg);
  #else
    return floorf(arg);
  #endif
#else
  return std::floor(arg);
#endif
}

CUDA_CALLABLE [[nodiscard]]
inline real_t ceil(real_t arg) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  #ifdef FP_64
    return ::ceil(arg);
  #else
    return ceilf(arg);
  #endif
#else
  return std::ceil(arg);
#endif
}

} // namespace