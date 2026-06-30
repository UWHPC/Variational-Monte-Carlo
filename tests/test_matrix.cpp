#include <catch2/catch_test_macros.hpp>

#include "utilities/matrix.hpp"

#include <array>
#include <cmath>
#include <cstddef>

TEST_CASE("LU decomposition solves a non-singular system", "[matrix]") {
  constexpr std::size_t N{3};
  std::array<double, N * N> A{
      2.0, 1.0, 1.0, 4.0, 3.0, 3.0, 8.0, 7.0, 9.0,
  };
  std::array<double, N> b{7.0, 19.0, 49.0};
  std::array<int, N> pivot{};
  std::array<double, N> x{};

  (void)lower_upper_decomp(A.data(), pivot.data(), N, N);
  solve_lower_upper(A.data(), pivot.data(), b.data(), x.data(), N, N);

  REQUIRE(std::abs(x[0] - 1.0) < 1e-9);
  REQUIRE(std::abs(x[1] - 2.0) < 1e-9);
  REQUIRE(std::abs(x[2] - 3.0) < 1e-9);
}

TEST_CASE("LU back-substitution stays finite on a singular matrix", "[matrix]") {
  constexpr std::size_t N{3};
  std::array<double, N * N> A{
      1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 1.0, 1.0, 1.0,
  };
  std::array<double, N> b{1.0, 1.0, 1.0};
  std::array<int, N> pivot{};
  std::array<double, N> x{};

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
