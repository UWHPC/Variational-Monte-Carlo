#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"
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
    Particles particles{N};
    SlaterPlaneWave slater{particles, 5.0};

    REQUIRE(slater.num_orbitals_get() == N);
    REQUIRE(slater.box_length_get() == 5.0);
    REQUIRE(slater.matrix_size_get() == N * N);
}

TEST_CASE("log_abs_det handles the N=1 constant orbital case", "[slater]") {
    // N=1: orbital 0 is k=0, cos → D = cos(0) = 1
    Particles particles{1U};
    SlaterPlaneWave slater{particles, 10.0};

    particles.pos_x_get()[0] = 3.25;
    particles.pos_y_get()[0] = 1.50;
    particles.pos_z_get()[0] = 7.75;

    const double logDet{slater.log_abs_det(particles)};

    requireNearSlater(logDet, 0.0);
    requireNearSlater(slater.determinant_get()[0], 1.0);
    requireNearSlater(slater.inv_determinant_get()[0], 1.0);
}

TEST_CASE("log_abs_det computes an inverse satisfying D*invD = I", "[slater]") {
    // N=3: orbital 0 = cos(0)=1, orbital 1 = cos(k1·r), orbital 2 = sin(k1·r)
    constexpr std::size_t N{3U};
    Particles particles{N};
    SlaterPlaneWave slater{particles, 11.0};

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
                         slater.inv_determinant_get()[matrixIndex(col, k, N)];
            }
            const double expected{row == col ? 1.0 : 0.0};
            requireNearSlater(value, expected, 1e-9);
        }
    }
}

TEST_CASE("N=3 determinant matrix uses cos/sin basis correctly", "[slater]") {
    // N=3: orbital 0 = cos(0)=1, orbital 1 = cos(k1·r), orbital 2 = sin(k1·r)
    // k1 is the first nonzero canonical k-vector
    constexpr std::size_t N{3U};
    constexpr double L{10.0};
    Particles particles{N};
    SlaterPlaneWave slater{particles, L};

    particles.pos_x_get()[0] = 1.0;
    particles.pos_y_get()[0] = 2.0;
    particles.pos_z_get()[0] = 3.0;

    particles.pos_x_get()[1] = 4.0;
    particles.pos_y_get()[1] = 5.0;
    particles.pos_z_get()[1] = 6.0;

    particles.pos_x_get()[2] = 7.0;
    particles.pos_y_get()[2] = 8.0;
    particles.pos_z_get()[2] = 9.0;

    slater.log_abs_det(particles);

    const auto& k_index{slater.orbital_k_index_get()};
    const auto& o_type{slater.orbital_type_get()};

    // Orbital 0 should be cos type with k=0
    REQUIRE(o_type[0] == 0);
    requireNearSlater(slater.k_vector_x_get()[k_index[0]], 0.0);
    requireNearSlater(slater.k_vector_y_get()[k_index[0]], 0.0);
    requireNearSlater(slater.k_vector_z_get()[k_index[0]], 0.0);

    // First column should all be 1.0 (cos(0))
    for (std::size_t i = 0; i < N; ++i) {
        requireNearSlater(slater.determinant_get()[matrixIndex(i, 0, N)], 1.0);
    }

    // Orbital 1 should be cos, orbital 2 should be sin, same k-vector
    REQUIRE(o_type[1] == 0);
    REQUIRE(o_type[2] == 1);
    REQUIRE(k_index[1] == k_index[2]);

    // Verify D entries for orbitals 1 and 2
    const std::size_t ki{k_index[1]};
    const double kx{slater.k_vector_x_get()[ki]};
    const double ky{slater.k_vector_y_get()[ki]};
    const double kz{slater.k_vector_z_get()[ki]};

    for (std::size_t i = 0; i < N; ++i) {
        const double k_dot_r{kx * particles.pos_x_get()[i] + ky * particles.pos_y_get()[i] +
                             kz * particles.pos_z_get()[i]};
        requireNearSlater(slater.determinant_get()[matrixIndex(i, 1, N)], std::cos(k_dot_r));
        requireNearSlater(slater.determinant_get()[matrixIndex(i, 2, N)], std::sin(k_dot_r));
    }
}

TEST_CASE("N=7 determinant is nonzero with cos/sin basis", "[slater]") {
    // N=7 is a closed shell: 1 (k=0) + 3 pairs × 2 = 7
    constexpr std::size_t N{7U};
    Particles particles{N};
    SlaterPlaneWave slater{particles, 10.0};

    // Spread particles around the box
    for (std::size_t i = 0; i < N; ++i) {
        particles.pos_x_get()[i] = 1.0 + static_cast<double>(i) * 1.1;
        particles.pos_y_get()[i] = 0.5 + static_cast<double>(i) * 0.7;
        particles.pos_z_get()[i] = 0.3 + static_cast<double>(i) * 1.3;
    }

    const double logDet{slater.log_abs_det(particles)};
    REQUIRE(std::isfinite(logDet));

    // Verify D * D^{-1} = I
    for (std::size_t row = 0; row < N; ++row) {
        for (std::size_t col = 0; col < N; ++col) {
            double value{};
            for (std::size_t k = 0; k < N; ++k) {
                value += slater.determinant_get()[matrixIndex(row, k, N)] *
                         slater.inv_determinant_get()[matrixIndex(col, k, N)];
            }
            const double expected{row == col ? 1.0 : 0.0};
            requireNearSlater(value, expected, 1e-9);
        }
    }
}

TEST_CASE("Slater derivatives match finite-difference for N=3 cos/sin basis", "[slater]") {
    constexpr std::size_t N{3U};
    constexpr double L{10.0};
    Particles particles{N};
    SlaterPlaneWave slater{particles, L};

    particles.pos_x_get()[0] = 1.1;
    particles.pos_y_get()[0] = 2.3;
    particles.pos_z_get()[0] = 0.7;

    particles.pos_x_get()[1] = 4.2;
    particles.pos_y_get()[1] = 1.8;
    particles.pos_z_get()[1] = 3.5;

    particles.pos_x_get()[2] = 7.6;
    particles.pos_y_get()[2] = 5.1;
    particles.pos_z_get()[2] = 8.9;

    // Compute analytic derivatives
    slater.log_abs_det(particles);
    const std::size_t stride{particles.padding_stride_get()};
    std::vector<double> gradX(stride, 0.0);
    std::vector<double> gradY(stride, 0.0);
    std::vector<double> gradZ(stride, 0.0);
    std::vector<double> lap(stride, 0.0);
    slater.add_derivatives(gradX.data(), gradY.data(), gradZ.data(), lap.data());

    // Finite-difference check for particle 0
    const double h{1e-5};
    const double center{slater.log_abs_det(particles)};

    auto shift_and_eval = [&](std::size_t p, double dx, double dy, double dz) {
        Particles shifted{N};
        for (std::size_t i = 0; i < N; ++i) {
            shifted.pos_x_get()[i] = particles.pos_x_get()[i];
            shifted.pos_y_get()[i] = particles.pos_y_get()[i];
            shifted.pos_z_get()[i] = particles.pos_z_get()[i];
        }
        shifted.pos_x_get()[p] += dx;
        shifted.pos_y_get()[p] += dy;
        shifted.pos_z_get()[p] += dz;
        return slater.log_abs_det(shifted);
    };

    for (std::size_t p = 0; p < N; ++p) {
        const double fd_gx{(shift_and_eval(p, h, 0.0, 0.0) - shift_and_eval(p, -h, 0.0, 0.0)) / (2.0 * h)};
        const double fd_gy{(shift_and_eval(p, 0.0, h, 0.0) - shift_and_eval(p, 0.0, -h, 0.0)) / (2.0 * h)};
        const double fd_gz{(shift_and_eval(p, 0.0, 0.0, h) - shift_and_eval(p, 0.0, 0.0, -h)) / (2.0 * h)};

        const double fd_lx{(shift_and_eval(p, h, 0.0, 0.0) - 2.0 * center + shift_and_eval(p, -h, 0.0, 0.0)) / (h * h)};
        const double fd_ly{(shift_and_eval(p, 0.0, h, 0.0) - 2.0 * center + shift_and_eval(p, 0.0, -h, 0.0)) / (h * h)};
        const double fd_lz{(shift_and_eval(p, 0.0, 0.0, h) - 2.0 * center + shift_and_eval(p, 0.0, 0.0, -h)) / (h * h)};

        requireNearSlater(gradX[p], fd_gx, 1e-6);
        requireNearSlater(gradY[p], fd_gy, 1e-6);
        requireNearSlater(gradZ[p], fd_gz, 1e-6);
        requireNearSlater(lap[p], fd_lx + fd_ly + fd_lz, 5e-4);
    }
}

// ── Shell-filling tests ──

TEST_CASE("Shell filling produces (0,0,0) as the first n-vector", "[slater]") {
    Particles p{1U};
    SlaterPlaneWave slater{p, 5.0};

    REQUIRE(slater.n_vector_x_get()[0] == 0);
    REQUIRE(slater.n_vector_y_get()[0] == 0);
    REQUIRE(slater.n_vector_z_get()[0] == 0);
}

TEST_CASE("Shell filling for N=7 uses canonical n-vectors with 4 unique k-vectors", "[slater]") {
    // N=7: 1 (k=0) + 3 nonzero k-vectors × 2 (cos,sin) = 7
    Particles p{7U};
    SlaterPlaneWave slater{p, 10.0};

    REQUIRE(slater.num_unique_k_get() == 4U);

    const int* n_x{slater.n_vector_x_get()};
    const int* n_y{slater.n_vector_y_get()};
    const int* n_z{slater.n_vector_z_get()};

    // k-index 0 should be (0,0,0)
    REQUIRE(n_x[0] == 0);
    REQUIRE(n_y[0] == 0);
    REQUIRE(n_z[0] == 0);

    // k-indices 1,2,3 should have |n|^2 = 1 and be canonical (first nonzero > 0)
    for (std::size_t i = 1; i < 4; ++i) {
        const int mag_sq{n_x[i] * n_x[i] + n_y[i] * n_y[i] + n_z[i] * n_z[i]};
        REQUIRE(mag_sq == 1);

        // Canonical: first nonzero component is positive
        if (n_x[i] != 0) {
            REQUIRE(n_x[i] > 0);
        } else if (n_y[i] != 0) {
            REQUIRE(n_y[i] > 0);
        } else {
            REQUIRE(n_z[i] > 0);
        }
    }
}

TEST_CASE("Shell filling orbital types alternate cos/sin for nonzero k", "[slater]") {
    Particles p{7U};
    SlaterPlaneWave slater{p, 10.0};

    const auto& o_type{slater.orbital_type_get()};
    const auto& k_index{slater.orbital_k_index_get()};

    // Orbital 0: cos (k=0)
    REQUIRE(o_type[0] == 0);

    // Orbitals 1-6: pairs of (cos, sin)
    for (std::size_t i = 1; i < 7; i += 2) {
        REQUIRE(o_type[i] == 0);               // cos
        REQUIRE(o_type[i + 1] == 1);           // sin
        REQUIRE(k_index[i] == k_index[i + 1]); // same k-vector
    }
}

TEST_CASE("Shell filling n-vectors are sorted by magnitude then lexicographically", "[slater]") {
    Particles p{7U};
    SlaterPlaneWave slater{p, 10.0};

    const std::size_t num_k{slater.num_unique_k_get()};
    const int* n_x{slater.n_vector_x_get()};
    const int* n_y{slater.n_vector_y_get()};
    const int* n_z{slater.n_vector_z_get()};

    for (std::size_t i = 0; i + 1 < num_k; ++i) {
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
    Particles p{N};
    SlaterPlaneWave slater{p, L};

    const double TWO_PI_OVER_L{2.0 * std::numbers::pi / L};
    const std::size_t num_k{slater.num_unique_k_get()};

    for (std::size_t i = 0; i < num_k; ++i) {
        requireNearSlater(slater.k_vector_x_get()[i], TWO_PI_OVER_L * static_cast<double>(slater.n_vector_x_get()[i]));
        requireNearSlater(slater.k_vector_y_get()[i], TWO_PI_OVER_L * static_cast<double>(slater.n_vector_y_get()[i]));
        requireNearSlater(slater.k_vector_z_get()[i], TWO_PI_OVER_L * static_cast<double>(slater.n_vector_z_get()[i]));
    }
}