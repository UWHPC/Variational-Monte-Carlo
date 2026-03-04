#include <catch2/catch_test_macros.hpp>

#include "blocking_analysis/blocking_analysis.hpp"

#include <cstddef>
#include <type_traits>
#include <utility>

TEST_CASE("BlockingAnalysis API is present", "[blocking]") {
    REQUIRE(std::is_constructible_v<BlockingAnalysis, std::size_t>);
    REQUIRE(std::is_same_v<decltype(std::declval<BlockingAnalysis&>().add(0.0)), void>);
    REQUIRE(std::is_same_v<decltype(std::declval<const BlockingAnalysis&>().ready()), bool>);
}
