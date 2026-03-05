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
    const PeriodicBoundaryCondition pbc{10.0};

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

    std::vector<double> expectedX(stride, 0.0);
    std::vector<double> expectedY(stride, 0.0);
    std::vector<double> expectedZ(stride, 0.0);
    std::vector<double> expectedLap(stride, 0.0);

    waveFunction.jastrow_pade_ptr().add_derivatives(particles, pbc, expectedX.data(), expectedY.data(),
                                                    expectedZ.data(), expectedLap.data());

    waveFunction.evaluate_derivatives(particles, pbc);

    for (std::size_t i = 0; i < stride; ++i) {
        requireNearWave(particles.grad_log_psi_x_get()[i], expectedX[i]);
        requireNearWave(particles.grad_log_psi_y_get()[i], expectedY[i]);
        requireNearWave(particles.grad_log_psi_z_get()[i], expectedZ[i]);
        requireNearWave(particles.lap_log_psi_get()[i], expectedLap[i]);
    }
}

TEST_CASE("WaveFunction evaluate_log_psi updates particle log_psi", "[wavefunction]") {
    Particles particles{1U};
    const PeriodicBoundaryCondition pbc{10.0};

    particles.pos_x_get()[0] = 0.4;
    particles.pos_y_get()[0] = 0.7;
    particles.pos_z_get()[0] = 0.9;

    WaveFunction waveFunction{1U, 10.0, 0.5, 1.0};

    waveFunction.slater_plane_wave_ptr().k_vector_x_get()[0] = 0.3;
    waveFunction.slater_plane_wave_ptr().k_vector_y_get()[0] = -0.2;
    waveFunction.slater_plane_wave_ptr().k_vector_z_get()[0] = 0.5;

    waveFunction.evaluate_log_psi(particles, pbc);

    const double k_dot_r{waveFunction.slater_plane_wave_ptr().k_vector_x_get()[0] * particles.pos_x_get()[0] +
                         waveFunction.slater_plane_wave_ptr().k_vector_y_get()[0] * particles.pos_y_get()[0] +
                         waveFunction.slater_plane_wave_ptr().k_vector_z_get()[0] * particles.pos_z_get()[0]};

    const double expectedLogPsi{std::log(std::abs(std::cos(k_dot_r)))};

    requireNearWave(particles.log_psi_get()[0], expectedLogPsi);
}