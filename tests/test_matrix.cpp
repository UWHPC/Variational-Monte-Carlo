#include <catch2/catch_test_macros.hpp>

#include "utilities/matrix.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#ifdef FP_64
constexpr real_t MATRIX_TOLERANCE{1e-9_r};
#else
constexpr real_t MATRIX_TOLERANCE{1e-4_r};
#endif

TEST_CASE("LU decomposition solves a non-singular system", "[matrix]") {
  constexpr std::size_t N{3};
  std::array<real_t, N * N> A{
      2.0_r, 1.0_r, 1.0_r, 4.0_r, 3.0_r, 3.0_r, 8.0_r, 7.0_r, 9.0_r,
  };
  std::array<real_t, N> b{7.0_r, 19.0_r, 49.0_r};
  std::array<int, N> pivot{};
  std::array<real_t, N> x{};

  (void)lower_upper_decomp(A.data(), pivot.data(), N, N);
  solve_lower_upper(A.data(), pivot.data(), b.data(), x.data(), N, N);

  REQUIRE(std::abs(x[0] - 1.0_r) < MATRIX_TOLERANCE);
  REQUIRE(std::abs(x[1] - 2.0_r) < MATRIX_TOLERANCE);
  REQUIRE(std::abs(x[2] - 3.0_r) < MATRIX_TOLERANCE);
}

TEST_CASE("LU back-substitution stays finite on a singular matrix", "[matrix]") {
  constexpr std::size_t N{3};
  std::array<real_t, N * N> A{
      1.0_r, 2.0_r, 3.0_r, 2.0_r, 4.0_r, 6.0_r, 1.0_r, 1.0_r, 1.0_r,
  };
  std::array<real_t, N> b{1.0_r, 1.0_r, 1.0_r};
  std::array<int, N> pivot{};
  std::array<real_t, N> x{};

  (void)lower_upper_decomp(A.data(), pivot.data(), N, N);
  solve_lower_upper(A.data(), pivot.data(), b.data(), x.data(), N, N);

  for (std::size_t i = 0; i < N; ++i) {
    REQUIRE(std::isfinite(x[i]));
  }
}

TEST_CASE("is_canonical keeps one representative of each +/-n pair", "[matrix]") {
  REQUIRE(is_canonical(0, 0, 0));

  REQUIRE(is_canonical(1, 0, 0));
  REQUIRE_FALSE(is_canonical(-1, 0, 0));

  REQUIRE(is_canonical(0, 2, -3));
  REQUIRE_FALSE(is_canonical(0, -2, 3));

  REQUIRE(is_canonical(0, 0, 5));
  REQUIRE_FALSE(is_canonical(0, 0, -5));
}
