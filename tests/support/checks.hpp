#pragma once

#include <catch2/catch_test_macros.hpp>

#include "utilities/macros.hpp"

#include <cmath>

#if defined(FP_64)

inline constexpr auto test_tolerance{1e-12_fp};

#else

inline constexpr auto test_tolerance{1e-6_fp};

#endif

inline void require_near(
  fp_t actual,
  fp_t expected,
  fp_t tolerance = test_tolerance
) {
  REQUIRE(std::abs(actual - expected) <= tolerance);
}
