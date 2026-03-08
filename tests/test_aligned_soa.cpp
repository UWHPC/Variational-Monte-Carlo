#include <catch2/catch_test_macros.hpp>

#include "utilities/aligned_soa.hpp"

#include <cstddef>
#include <cstdint>

TEST_CASE("AlignedSoA rounds stride to SIMD alignment and zero-initializes", "[aligned_soa]") {
    AlignedSoA<double> soa{3U, 2U};

    const std::size_t stride{soa.stride()};
    const std::size_t doubles_per_alignment{SIMD_BYTES / sizeof(double)};

    REQUIRE(soa.num_elements() == 3U);
    REQUIRE(stride >= 3U);
    REQUIRE(stride % doubles_per_alignment == 0U);

    const auto address{reinterpret_cast<std::uintptr_t>(soa[0])};
    REQUIRE(address % SIMD_BYTES == 0U);

    for (std::size_t i = 0; i < stride; ++i) {
        REQUIRE(soa[0][i] == 0.0);
        REQUIRE(soa[1][i] == 0.0);
    }
}

TEST_CASE("AlignedSoA exposes non-overlapping fixed-stride slices", "[aligned_soa]") {
    AlignedSoA<double> soa{5U, 3U};

    const std::ptrdiff_t stride{static_cast<std::ptrdiff_t>(soa.stride())};
    REQUIRE(soa[1] - soa[0] == stride);
    REQUIRE(soa[2] - soa[1] == stride);

    soa[0][0] = 1.0;
    soa[1][0] = 2.0;
    soa[2][0] = 3.0;

    REQUIRE(soa[0][0] == 1.0);
    REQUIRE(soa[1][0] == 2.0);
    REQUIRE(soa[2][0] == 3.0);
}

TEST_CASE("AlignedSoA const operator[] reads previously written values", "[aligned_soa]") {
    AlignedSoA<double> soa{4U, 1U};
    soa[0][0] = 9.5;
    soa[0][1] = -1.25;

    const AlignedSoA<double>& const_soa{soa};
    REQUIRE(const_soa[0][0] == 9.5);
    REQUIRE(const_soa[0][1] == -1.25);
}
