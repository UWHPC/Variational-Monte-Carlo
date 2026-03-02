#include <catch2/catch_test_macros.hpp>

TEST_CASE("Sanity test: build and test runner works", "[sanity]") {
    REQUIRE(1 + 1 == 2);
}

TEST_CASE("Floating point sanity", "[sanity]") {
    const double x = 0.1 + 0.2;
    REQUIRE(x > 0.0);
}
