#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"
#include "pbc/pbc.hpp"
#include "slater_plane_wave/slater_plane_wave.hpp"

#include <cmath>
#include <numbers>
#include <tuple>
#include <vector>

namespace {

void requireNearSlater(double actual, double expected, double tolerance = 1e-10) {
    REQUIRE(std::abs(actual - expected) <= tolerance);
}

std::size_t matrixIndex(std::size_t row, std::size_t col, std::size_t n) { return row * n + col; }

} // namespace

TEST_CASE("SlaterPlaneWave constructor initializes correctly", "[slater]") {
    constexpr std::size_t N{3U};
    SlaterPlaneWave slater{N, 5.0};

    REQUIRE(slater.num_orbitals_get() == N);
    REQUIRE(slater.box_length_get() == 5.0);
    REQUIRE(slater.matrix_size_get() == N * N);
}

TEST_CASE("log_abs_det handles the N=1 constant orbital case", "[slater]") {
    SlaterPlaneWave slater{1U, 10.0};
    Particles particles{1U};

    slater.k_vector_x_get()[0] = 0.0;
    slater.k_vector_y_get()[0] = 0.0;
    slater.k_vector_z_get()[0] = 0.0;

    particles.pos_x_get()[0] = 3.25;
    particles.pos_y_get()[0] = 1.50;
    particles.pos_z_get()[0] = 7.75;

    const double logDet{slater.log_abs_det(particles)};

    requireNearSlater(logDet, 0.0);
    requireNearSlater(slater.determinant_get()[0], 1.0);
    requireNearSlater(slater.inv_determinant_get()[0], 1.0);
}

TEST_CASE("log_abs_det returns negative infinity for singular matrices", "[slater]") {
    SlaterPlaneWave slater{2U, 10.0};
    Particles particles{2U};

    slater.k_vector_x_get()[0] = 0.0;
    slater.k_vector_y_get()[0] = 0.0;
    slater.k_vector_z_get()[0] = 0.0;
    slater.k_vector_x_get()[1] = 0.0;
    slater.k_vector_y_get()[1] = 0.0;
    slater.k_vector_z_get()[1] = 0.0;

    particles.pos_x_get()[0] = 0.0;
    particles.pos_y_get()[0] = 0.0;
    particles.pos_z_get()[0] = 0.0;
    particles.pos_x_get()[1] = 1.0;
    particles.pos_y_get()[1] = 0.0;
    particles.pos_z_get()[1] = 0.0;

    const double logDet{slater.log_abs_det(particles)};

    REQUIRE(std::isinf(logDet));
    REQUIRE(logDet < 0.0);
}

TEST_CASE("log_abs_det computes determinant and inverse for a known 2x2 system", "[slater]") {
    SlaterPlaneWave slater{2U, 10.0};
    Particles particles{2U};

    slater.k_vector_x_get()[0] = 0.0;
    slater.k_vector_y_get()[0] = 0.0;
    slater.k_vector_z_get()[0] = 0.0;

    slater.k_vector_x_get()[1] = 1.0;
    slater.k_vector_y_get()[1] = 0.0;
    slater.k_vector_z_get()[1] = 0.0;

    particles.pos_x_get()[0] = 0.0;
    particles.pos_y_get()[0] = 0.0;
    particles.pos_z_get()[0] = 0.0;

    particles.pos_x_get()[1] = std::numbers::pi_v<double> * 0.5;
    particles.pos_y_get()[1] = 0.0;
    particles.pos_z_get()[1] = 0.0;

    const double logDet{slater.log_abs_det(particles)};

    requireNearSlater(logDet, 0.0, 1e-9);

    requireNearSlater(slater.determinant_get()[0], 1.0);
    requireNearSlater(slater.determinant_get()[1], 1.0);
    requireNearSlater(slater.determinant_get()[2], 1.0);
    requireNearSlater(slater.determinant_get()[3], 0.0, 1e-12);

    requireNearSlater(slater.inv_determinant_get()[0], 0.0, 1e-9);
    requireNearSlater(slater.inv_determinant_get()[1], 1.0, 1e-9);
    requireNearSlater(slater.inv_determinant_get()[2], 1.0, 1e-9);
    requireNearSlater(slater.inv_determinant_get()[3], -1.0, 1e-9);
}

TEST_CASE("log_abs_det computes an inverse satisfying D*invD = I", "[slater]") {
    constexpr std::size_t N{3U};
    SlaterPlaneWave slater{N, 11.0};
    Particles particles{N};

    slater.k_vector_x_get()[0] = 0.0;
    slater.k_vector_y_get()[0] = 0.0;
    slater.k_vector_z_get()[0] = 0.0;

    slater.k_vector_x_get()[1] = 1.0;
    slater.k_vector_y_get()[1] = 0.5;
    slater.k_vector_z_get()[1] = 0.0;

    slater.k_vector_x_get()[2] = -0.25;
    slater.k_vector_y_get()[2] = 1.1;
    slater.k_vector_z_get()[2] = 0.7;

    particles.pos_x_get()[0] = 0.3;
    particles.pos_y_get()[0] = 0.4;
    particles.pos_z_get()[0] = 0.5;

    particles.pos_x_get()[1] = 1.7;
    particles.pos_y_get()[1] = 0.2;
    particles.pos_z_get()[1] = 2.1;

    particles.pos_x_get()[2] = 2.2;
    particles.pos_y_get()[2] = 1.8;
    particles.pos_z_get()[2] = 0.9;

    const double logDet{slater.log_abs_det(particles)};
    REQUIRE(std::isfinite(logDet));

    for (std::size_t row = 0; row < N; ++row) {
        for (std::size_t col = 0; col < N; ++col) {
            double value{};
            for (std::size_t k = 0; k < N; ++k) {
                value += slater.determinant_get()[matrixIndex(row, k, N)] *
                         slater.inv_determinant_get()[matrixIndex(k, col, N)];
            }
            const double expected{row == col ? 1.0 : 0.0};
            requireNearSlater(value, expected, 1e-9);
        }
    }
}

TEST_CASE("Slater derivatives match analytic N=1 formulas", "[slater]") {
    SlaterPlaneWave slater{1U, 10.0};
    Particles particles{1U};

    slater.k_vector_x_get()[0] = 1.2;
    slater.k_vector_y_get()[0] = -0.4;
    slater.k_vector_z_get()[0] = 0.7;

    particles.pos_x_get()[0] = 0.3;
    particles.pos_y_get()[0] = 0.5;
    particles.pos_z_get()[0] = 0.9;

    const double logDet{slater.log_abs_det(particles)};
    REQUIRE(std::isfinite(logDet));

    const std::size_t stride{particles.padding_stride_get()};
    std::vector<double> gradX(stride, 2.0);
    std::vector<double> gradY(stride, 2.0);
    std::vector<double> gradZ(stride, 2.0);
    std::vector<double> lap(stride, -3.0);

    slater.add_derivatives(particles, gradX.data(), gradY.data(), gradZ.data(), lap.data());

    const double kdotr{slater.k_vector_x_get()[0] * particles.pos_x_get()[0] +
                       slater.k_vector_y_get()[0] * particles.pos_y_get()[0] +
                       slater.k_vector_z_get()[0] * particles.pos_z_get()[0]};

    const double tangent{std::tan(kdotr)};
    const double kSquared{slater.k_vector_x_get()[0] * slater.k_vector_x_get()[0] +
                          slater.k_vector_y_get()[0] * slater.k_vector_y_get()[0] +
                          slater.k_vector_z_get()[0] * slater.k_vector_z_get()[0]};

    const double cosine{std::cos(kdotr)};

    const double expectedGradX{-tangent * slater.k_vector_x_get()[0]};
    const double expectedGradY{-tangent * slater.k_vector_y_get()[0]};
    const double expectedGradZ{-tangent * slater.k_vector_z_get()[0]};
    const double expectedLap{-kSquared / (cosine * cosine)};

    requireNearSlater(gradX[0], 2.0 + expectedGradX, 1e-10);
    requireNearSlater(gradY[0], 2.0 + expectedGradY, 1e-10);
    requireNearSlater(gradZ[0], 2.0 + expectedGradZ, 1e-10);
    requireNearSlater(lap[0], -3.0 + expectedLap, 1e-10);

    for (std::size_t i = 1; i < stride; ++i) {
        requireNearSlater(gradX[i], 2.0);
        requireNearSlater(gradY[i], 2.0);
        requireNearSlater(gradZ[i], 2.0);
        requireNearSlater(lap[i], -3.0);
    }
}

// ── Shell-filling tests ──

TEST_CASE("Shell filling produces (0,0,0) as the first n-vector", "[slater]") {
    SlaterPlaneWave slater{1U, 5.0};

    REQUIRE(slater.n_vector_x_get()[0] == 0);
    REQUIRE(slater.n_vector_y_get()[0] == 0);
    REQUIRE(slater.n_vector_z_get()[0] == 0);
}

TEST_CASE("Shell filling for N=7 gives a non-singular Slater determinant", "[slater]") {
    constexpr std::size_t N{7U};
    constexpr double L{10.0};
    SlaterPlaneWave slater{N, L};
    Particles particles{N};

    // Positions spread across the box in all three dimensions
    particles.pos_x_get()[0] = 1.0;
    particles.pos_y_get()[0] = 2.0;
    particles.pos_z_get()[0] = 3.0;
    particles.pos_x_get()[1] = 4.5;
    particles.pos_y_get()[1] = 0.5;
    particles.pos_z_get()[1] = 7.2;
    particles.pos_x_get()[2] = 8.1;
    particles.pos_y_get()[2] = 6.3;
    particles.pos_z_get()[2] = 1.4;
    particles.pos_x_get()[3] = 2.7;
    particles.pos_y_get()[3] = 8.8;
    particles.pos_z_get()[3] = 5.5;
    particles.pos_x_get()[4] = 6.0;
    particles.pos_y_get()[4] = 3.7;
    particles.pos_z_get()[4] = 9.1;
    particles.pos_x_get()[5] = 0.3;
    particles.pos_y_get()[5] = 5.1;
    particles.pos_z_get()[5] = 4.8;
    particles.pos_x_get()[6] = 7.4;
    particles.pos_y_get()[6] = 9.2;
    particles.pos_z_get()[6] = 0.6;

    const double logDet{slater.log_abs_det(particles)};
    REQUIRE(std::isfinite(logDet));
}

TEST_CASE("Shell filling sorts by magnitude then lexicographically", "[slater]") {
    SlaterPlaneWave slater{7U, 10.0};

    const int* n_x{slater.n_vector_x_get()};
    const int* n_y{slater.n_vector_y_get()};
    const int* n_z{slater.n_vector_z_get()};

    for (std::size_t i = 0; i + 1 < 7; ++i) {
        const int mag_sq_a{n_x[i] * n_x[i] + n_y[i] * n_y[i] + n_z[i] * n_z[i]};
        const int mag_sq_b{n_x[i + 1] * n_x[i + 1] + n_y[i + 1] * n_y[i + 1] + n_z[i + 1] * n_z[i + 1]};

        const bool magnitude_ok{mag_sq_a <= mag_sq_b};
        REQUIRE(magnitude_ok);

        if (mag_sq_a == mag_sq_b) {
            const bool lex_ok{std::tie(n_x[i], n_y[i], n_z[i]) <= std::tie(n_x[i + 1], n_y[i + 1], n_z[i + 1])};
            REQUIRE(lex_ok);
        }
    }
}

TEST_CASE("Shell filling k-vectors match 2pi/L times n-vectors", "[slater]") {
    constexpr std::size_t N{7U};
    constexpr double L{8.0};
    SlaterPlaneWave slater{N, L};

    const double TWO_PI_OVER_L{2.0 * std::numbers::pi / L};

    for (std::size_t i = 0; i < N; ++i) {
        requireNearSlater(slater.k_vector_x_get()[i], TWO_PI_OVER_L * static_cast<double>(slater.n_vector_x_get()[i]));
        requireNearSlater(slater.k_vector_y_get()[i], TWO_PI_OVER_L * static_cast<double>(slater.n_vector_y_get()[i]));
        requireNearSlater(slater.k_vector_z_get()[i], TWO_PI_OVER_L * static_cast<double>(slater.n_vector_z_get()[i]));
    }
}

TEST_CASE("Shell filling for N=7 gives a non-singular Slater determinant", "[slater]") {
    constexpr std::size_t N{7U};
    SlaterPlaneWave slater{N, 10.0};
    Particles particles{N};

    for (std::size_t i = 0; i < N; ++i) {
        particles.pos_x_get()[i] = static_cast<double>(i) * 1.3;
        particles.pos_y_get()[i] = static_cast<double>(i) * 0.7;
        particles.pos_z_get()[i] = static_cast<double>(i) * 0.9;
    }

    const double logDet{slater.log_abs_det(particles)};
    REQUIRE(std::isfinite(logDet));
}