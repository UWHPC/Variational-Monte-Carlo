#include <catch2/catch_test_macros.hpp>

#include "jastrow_pade/jastrow_pade.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace {

void requireNearJastrow(double actual, double expected, double tolerance = 1e-12) {
    REQUIRE(std::abs(actual - expected) <= tolerance);
}

Particles copyParticlePositions(const Particles& source) {
    Particles copy{source.num_particles_get()};
    const std::size_t numParticles{source.num_particles_get()};
    std::copy_n(source.pos_x_get(), numParticles, copy.pos_x_get());
    std::copy_n(source.pos_y_get(), numParticles, copy.pos_y_get());
    std::copy_n(source.pos_z_get(), numParticles, copy.pos_z_get());
    return copy;
}

double valueAtOffset(const JastrowPade& jastrow, const Particles& reference, std::size_t particle,
                     double dx, double dy, double dz) {
    Particles shifted{copyParticlePositions(reference)};
    shifted.pos_x_get()[particle] += dx;
    shifted.pos_y_get()[particle] += dy;
    shifted.pos_z_get()[particle] += dz;
    return jastrow.value(shifted);
}

} // namespace

TEST_CASE("Jastrow value uses minimum-image pair distances", "[jastrow]") {
    const JastrowPade jastrow{10.0, 0.25, 1.0};
    Particles particles{2U};

    particles.pos_x_get()[0] = 0.1;
    particles.pos_y_get()[0] = 0.0;
    particles.pos_z_get()[0] = 0.0;

    particles.pos_x_get()[1] = 9.9;
    particles.pos_y_get()[1] = 0.0;
    particles.pos_z_get()[1] = 0.0;

    const double r{0.2};
    const double expected{(0.25 * r) / (1.0 + r)};
    requireNearJastrow(jastrow.value(particles), expected);
}

TEST_CASE("Jastrow value skips degenerate pairs", "[jastrow]") {
    const JastrowPade jastrow{10.0, 0.25, 1.0};
    Particles particles{2U};

    particles.pos_x_get()[0] = 1.0;
    particles.pos_y_get()[0] = 2.0;
    particles.pos_z_get()[0] = 3.0;

    particles.pos_x_get()[1] = 1.0;
    particles.pos_y_get()[1] = 2.0;
    particles.pos_z_get()[1] = 3.0;

    requireNearJastrow(jastrow.value(particles), 0.0);
}

TEST_CASE("Jastrow derivatives match the analytic two-particle result", "[jastrow]") {
    const JastrowPade jastrow{100.0, 0.25, 1.0};
    Particles particles{2U};

    particles.pos_x_get()[0] = 0.0;
    particles.pos_y_get()[0] = 0.0;
    particles.pos_z_get()[0] = 0.0;

    particles.pos_x_get()[1] = 1.0;
    particles.pos_y_get()[1] = 0.0;
    particles.pos_z_get()[1] = 0.0;

    const std::size_t stride{particles.padding_stride_get()};
    std::vector<double> gradX(stride, 0.0);
    std::vector<double> gradY(stride, 0.0);
    std::vector<double> gradZ(stride, 0.0);
    std::vector<double> lap(stride, 0.0);

    jastrow.add_derivatives(particles, gradX.data(), gradY.data(), gradZ.data(), lap.data());

    requireNearJastrow(gradX[0], -0.125);
    requireNearJastrow(gradX[1], 0.125);
    requireNearJastrow(gradX[0] + gradX[1], 0.0);

    requireNearJastrow(gradY[0], 0.0);
    requireNearJastrow(gradY[1], 0.0);
    requireNearJastrow(gradZ[0], 0.0);
    requireNearJastrow(gradZ[1], 0.0);

    requireNearJastrow(lap[0], 0.125);
    requireNearJastrow(lap[1], 0.125);

    for (std::size_t i = 2; i < stride; ++i) {
        requireNearJastrow(gradX[i], 0.0);
        requireNearJastrow(gradY[i], 0.0);
        requireNearJastrow(gradZ[i], 0.0);
        requireNearJastrow(lap[i], 0.0);
    }
}

TEST_CASE("Jastrow derivatives are unchanged for degenerate pairs", "[jastrow]") {
    const JastrowPade jastrow{10.0, 0.25, 1.0};
    Particles particles{2U};

    particles.pos_x_get()[0] = 4.0;
    particles.pos_y_get()[0] = 5.0;
    particles.pos_z_get()[0] = 6.0;

    particles.pos_x_get()[1] = 4.0;
    particles.pos_y_get()[1] = 5.0;
    particles.pos_z_get()[1] = 6.0;

    const std::size_t stride{particles.padding_stride_get()};
    std::vector<double> gradX(stride, 3.0);
    std::vector<double> gradY(stride, -2.0);
    std::vector<double> gradZ(stride, 1.5);
    std::vector<double> lap(stride, 7.0);

    jastrow.add_derivatives(particles, gradX.data(), gradY.data(), gradZ.data(), lap.data());

    for (std::size_t i = 0; i < stride; ++i) {
        requireNearJastrow(gradX[i], 3.0);
        requireNearJastrow(gradY[i], -2.0);
        requireNearJastrow(gradZ[i], 1.5);
        requireNearJastrow(lap[i], 7.0);
    }
}

TEST_CASE("Jastrow derivatives match finite-difference gradients and Laplacians", "[jastrow]") {
    const JastrowPade jastrow{20.0, 0.6, 0.9};
    Particles particles{3U};

    particles.pos_x_get()[0] = 1.1;
    particles.pos_y_get()[0] = 2.2;
    particles.pos_z_get()[0] = 0.7;

    particles.pos_x_get()[1] = 3.8;
    particles.pos_y_get()[1] = 1.4;
    particles.pos_z_get()[1] = 2.5;

    particles.pos_x_get()[2] = 0.2;
    particles.pos_y_get()[2] = 4.1;
    particles.pos_z_get()[2] = 3.3;

    const std::size_t stride{particles.padding_stride_get()};
    std::vector<double> gradX(stride, 0.0);
    std::vector<double> gradY(stride, 0.0);
    std::vector<double> gradZ(stride, 0.0);
    std::vector<double> lap(stride, 0.0);
    jastrow.add_derivatives(particles, gradX.data(), gradY.data(), gradZ.data(), lap.data());

    const double h{1e-5};
    const double valueCenter{jastrow.value(particles)};

    const double dJdx{(valueAtOffset(jastrow, particles, 0U, h, 0.0, 0.0) -
                       valueAtOffset(jastrow, particles, 0U, -h, 0.0, 0.0)) /
                      (2.0 * h)};
    const double dJdy{(valueAtOffset(jastrow, particles, 0U, 0.0, h, 0.0) -
                       valueAtOffset(jastrow, particles, 0U, 0.0, -h, 0.0)) /
                      (2.0 * h)};
    const double dJdz{(valueAtOffset(jastrow, particles, 0U, 0.0, 0.0, h) -
                       valueAtOffset(jastrow, particles, 0U, 0.0, 0.0, -h)) /
                      (2.0 * h)};

    const double d2Jdx2{(valueAtOffset(jastrow, particles, 0U, h, 0.0, 0.0) - 2.0 * valueCenter +
                         valueAtOffset(jastrow, particles, 0U, -h, 0.0, 0.0)) /
                        (h * h)};
    const double d2Jdy2{(valueAtOffset(jastrow, particles, 0U, 0.0, h, 0.0) - 2.0 * valueCenter +
                         valueAtOffset(jastrow, particles, 0U, 0.0, -h, 0.0)) /
                        (h * h)};
    const double d2Jdz2{(valueAtOffset(jastrow, particles, 0U, 0.0, 0.0, h) - 2.0 * valueCenter +
                         valueAtOffset(jastrow, particles, 0U, 0.0, 0.0, -h)) /
                        (h * h)};

    requireNearJastrow(gradX[0], dJdx, 1e-7);
    requireNearJastrow(gradY[0], dJdy, 1e-7);
    requireNearJastrow(gradZ[0], dJdz, 1e-7);
    requireNearJastrow(lap[0], d2Jdx2 + d2Jdy2 + d2Jdz2, 2e-4);

    requireNearJastrow(gradX[0] + gradX[1] + gradX[2], 0.0, 1e-12);
    requireNearJastrow(gradY[0] + gradY[1] + gradY[2], 0.0, 1e-12);
    requireNearJastrow(gradZ[0] + gradZ[1] + gradZ[2], 0.0, 1e-12);
}