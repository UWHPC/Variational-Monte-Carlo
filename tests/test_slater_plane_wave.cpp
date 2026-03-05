#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"
#include "pbc/pbc.hpp"
#include "slater_plane_wave/slater_plane_wave.hpp"

#include <cmath>
#include <numbers>
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

    REQUIRE(slater.num_orbitals_ptr() == N);
    REQUIRE(slater.box_length_ptr() == 5.0);
    REQUIRE(slater.matrix_size_ptr() == N * N);
}

TEST_CASE("log_abs_det handles the N=1 constant orbital case", "[slater]") {
    SlaterPlaneWave slater{1U, 10.0};
    Particles particles{1U};

    slater.k_vector_x_ptr()[0] = 0.0;
    slater.k_vector_y_ptr()[0] = 0.0;
    slater.k_vector_z_ptr()[0] = 0.0;

    particles.pos_x_ptr()[0] = 3.25;
    particles.pos_y_ptr()[0] = 1.50;
    particles.pos_z_ptr()[0] = 7.75;

    const double logDet{slater.log_abs_det(particles)};

    requireNearSlater(logDet, 0.0);
    requireNearSlater(slater.determinant_ptr()[0], 1.0);
    requireNearSlater(slater.inv_determinant_ptr()[0], 1.0);
}

TEST_CASE("log_abs_det returns negative infinity for singular matrices", "[slater]") {
    SlaterPlaneWave slater{2U, 10.0};
    Particles particles{2U};

    slater.k_vector_x_ptr()[0] = 0.0;
    slater.k_vector_y_ptr()[0] = 0.0;
    slater.k_vector_z_ptr()[0] = 0.0;
    slater.k_vector_x_ptr()[1] = 0.0;
    slater.k_vector_y_ptr()[1] = 0.0;
    slater.k_vector_z_ptr()[1] = 0.0;

    particles.pos_x_ptr()[0] = 0.0;
    particles.pos_y_ptr()[0] = 0.0;
    particles.pos_z_ptr()[0] = 0.0;
    particles.pos_x_ptr()[1] = 1.0;
    particles.pos_y_ptr()[1] = 0.0;
    particles.pos_z_ptr()[1] = 0.0;

    const double logDet{slater.log_abs_det(particles)};

    REQUIRE(std::isinf(logDet));
    REQUIRE(logDet < 0.0);
}

TEST_CASE("log_abs_det computes determinant and inverse for a known 2x2 system", "[slater]") {
    SlaterPlaneWave slater{2U, 10.0};
    Particles particles{2U};

    slater.k_vector_x_ptr()[0] = 0.0;
    slater.k_vector_y_ptr()[0] = 0.0;
    slater.k_vector_z_ptr()[0] = 0.0;

    slater.k_vector_x_ptr()[1] = 1.0;
    slater.k_vector_y_ptr()[1] = 0.0;
    slater.k_vector_z_ptr()[1] = 0.0;

    particles.pos_x_ptr()[0] = 0.0;
    particles.pos_y_ptr()[0] = 0.0;
    particles.pos_z_ptr()[0] = 0.0;

    particles.pos_x_ptr()[1] = std::numbers::pi_v<double> * 0.5;
    particles.pos_y_ptr()[1] = 0.0;
    particles.pos_z_ptr()[1] = 0.0;

    const double logDet{slater.log_abs_det(particles)};

    requireNearSlater(logDet, 0.0, 1e-9);

    requireNearSlater(slater.determinant_ptr()[0], 1.0);
    requireNearSlater(slater.determinant_ptr()[1], 1.0);
    requireNearSlater(slater.determinant_ptr()[2], 1.0);
    requireNearSlater(slater.determinant_ptr()[3], 0.0, 1e-12);

    requireNearSlater(slater.inv_determinant_ptr()[0], 0.0, 1e-9);
    requireNearSlater(slater.inv_determinant_ptr()[1], 1.0, 1e-9);
    requireNearSlater(slater.inv_determinant_ptr()[2], 1.0, 1e-9);
    requireNearSlater(slater.inv_determinant_ptr()[3], -1.0, 1e-9);
}

TEST_CASE("log_abs_det computes an inverse satisfying D*invD = I", "[slater]") {
    constexpr std::size_t N{3U};
    SlaterPlaneWave slater{N, 11.0};
    Particles particles{N};

    slater.k_vector_x_ptr()[0] = 0.0;
    slater.k_vector_y_ptr()[0] = 0.0;
    slater.k_vector_z_ptr()[0] = 0.0;

    slater.k_vector_x_ptr()[1] = 1.0;
    slater.k_vector_y_ptr()[1] = 0.5;
    slater.k_vector_z_ptr()[1] = 0.0;

    slater.k_vector_x_ptr()[2] = -0.25;
    slater.k_vector_y_ptr()[2] = 1.1;
    slater.k_vector_z_ptr()[2] = 0.7;

    particles.pos_x_ptr()[0] = 0.3;
    particles.pos_y_ptr()[0] = 0.4;
    particles.pos_z_ptr()[0] = 0.5;

    particles.pos_x_ptr()[1] = 1.7;
    particles.pos_y_ptr()[1] = 0.2;
    particles.pos_z_ptr()[1] = 2.1;

    particles.pos_x_ptr()[2] = 2.2;
    particles.pos_y_ptr()[2] = 1.8;
    particles.pos_z_ptr()[2] = 0.9;

    const double logDet{slater.log_abs_det(particles)};
    REQUIRE(std::isfinite(logDet));

    for (std::size_t row = 0; row < N; ++row) {
        for (std::size_t col = 0; col < N; ++col) {
            double value{};
            for (std::size_t k = 0; k < N; ++k) {
                value += slater.determinant_ptr()[matrixIndex(row, k, N)] *
                         slater.inv_determinant_ptr()[matrixIndex(k, col, N)];
            }
            const double expected{row == col ? 1.0 : 0.0};
            requireNearSlater(value, expected, 1e-9);
        }
    }
}

TEST_CASE("Slater derivatives match analytic N=1 formulas", "[slater]") {
    SlaterPlaneWave slater{1U, 10.0};
    Particles particles{1U};

    slater.k_vector_x_ptr()[0] = 1.2;
    slater.k_vector_y_ptr()[0] = -0.4;
    slater.k_vector_z_ptr()[0] = 0.7;

    particles.pos_x_ptr()[0] = 0.3;
    particles.pos_y_ptr()[0] = 0.5;
    particles.pos_z_ptr()[0] = 0.9;

    const double logDet{slater.log_abs_det(particles)};
    REQUIRE(std::isfinite(logDet));

    const std::size_t stride{particles.padding_stride_ptr()};
    std::vector<double> gradX(stride, 2.0);
    std::vector<double> gradY(stride, 2.0);
    std::vector<double> gradZ(stride, 2.0);
    std::vector<double> lap(stride, -3.0);

    slater.add_derivatives(particles, gradX.data(), gradY.data(), gradZ.data(), lap.data());

    const double kdotr{slater.k_vector_x_ptr()[0] * particles.pos_x_ptr()[0] +
                       slater.k_vector_y_ptr()[0] * particles.pos_y_ptr()[0] +
                       slater.k_vector_z_ptr()[0] * particles.pos_z_ptr()[0]};

    const double tangent{std::tan(kdotr)};
    const double kSquared{slater.k_vector_x_ptr()[0] * slater.k_vector_x_ptr()[0] +
                          slater.k_vector_y_ptr()[0] * slater.k_vector_y_ptr()[0] +
                          slater.k_vector_z_ptr()[0] * slater.k_vector_z_ptr()[0]};

    const double cosine{std::cos(kdotr)};

    const double expectedGradX{-tangent * slater.k_vector_x_ptr()[0]};
    const double expectedGradY{-tangent * slater.k_vector_y_ptr()[0]};
    const double expectedGradZ{-tangent * slater.k_vector_z_ptr()[0]};
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
