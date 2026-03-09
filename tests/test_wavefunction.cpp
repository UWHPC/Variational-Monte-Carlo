#include "wavefunction/wavefunction.hpp"
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>

namespace {
void requireNearWave(double actual, double expected, double tolerance = 1e-12) {
    REQUIRE(std::abs(actual - expected) <= tolerance);
}
} // namespace

TEST_CASE("WaveFunction evaluateDerivatives clears buffers and delegates to Jastrow", "[wavefunction]") {
    Particles particles{2U};

    particles.pos_x_get()[0] = 0.0;
    particles.pos_y_get()[0] = 0.0;
    particles.pos_z_get()[0] = 0.0;

    particles.pos_x_get()[1] = 1.0;
    particles.pos_y_get()[1] = 0.0;
    particles.pos_z_get()[1] = 0.0;

    const std::size_t stride{particles.padding_stride_get()};
    for (std::size_t i = 0; i < stride; ++i) {
        particles.grad_log_psi_x_get()[i] = 999.0;
        particles.grad_log_psi_y_get()[i] = 999.0;
        particles.grad_log_psi_z_get()[i] = 999.0;
        particles.lap_log_psi_get()[i] = 999.0;
    }

    WaveFunction waveFunction{2U, 10.0, 0.5, 1.0};

    // Compute expected: zero-init then add both Slater and Jastrow
    std::vector<double> expectedX(stride, 0.0);
    std::vector<double> expectedY(stride, 0.0);
    std::vector<double> expectedZ(stride, 0.0);
    std::vector<double> expectedLap(stride, 0.0);

    // Need to call log_abs_det first to populate the inverse
    waveFunction.slater_plane_wave_get().log_abs_det(particles);
    waveFunction.slater_plane_wave_get().add_derivatives(particles, expectedX.data(), expectedY.data(),
                                                         expectedZ.data(), expectedLap.data());
    waveFunction.jastrow_pade_get().add_derivatives(particles, expectedX.data(), expectedY.data(), expectedZ.data(),
                                                    expectedLap.data());

    waveFunction.evaluate_derivatives(particles);

    for (std::size_t i = 0; i < stride; ++i) {
        requireNearWave(particles.grad_log_psi_x_get()[i], expectedX[i]);
        requireNearWave(particles.grad_log_psi_y_get()[i], expectedY[i]);
        requireNearWave(particles.grad_log_psi_z_get()[i], expectedZ[i]);
        requireNearWave(particles.lap_log_psi_get()[i], expectedLap[i]);
    }
}

TEST_CASE("WaveFunction evaluate_log_psi returns finite for N=1", "[wavefunction]") {
    // N=1: orbital 0 is cos(0·r) = 1, so log|det| = 0
    // Jastrow with 1 particle has no pairs, so J = 0
    // Total log_psi = 0
    Particles particles{1U};

    particles.pos_x_get()[0] = 0.4;
    particles.pos_y_get()[0] = 0.7;
    particles.pos_z_get()[0] = 0.9;

    WaveFunction waveFunction{1U, 10.0, 0.5, 1.0};

    const double logPsi{waveFunction.evaluate_log_psi(particles)};
    requireNearWave(logPsi, 0.0);
}