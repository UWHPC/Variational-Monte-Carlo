#include <catch2/catch_test_macros.hpp>

#include "simulation/simulation.hpp"

#include <type_traits>
#include <utility>

TEST_CASE("Simulation API is present", "[simulation]") {
    REQUIRE(std::is_constructible_v<Simulation, Config>);
    REQUIRE(std::is_same_v<decltype(std::declval<Simulation&>().run()), void>);
}
