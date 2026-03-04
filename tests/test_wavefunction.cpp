#include <catch2/catch_test_macros.hpp>

#include "wavefunction/wavefunction.hpp"

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

    particles.posX()[0] = 0.0;
    particles.posY()[0] = 0.0;
    particles.posZ()[0] = 0.0;

    particles.posX()[1] = 1.0;
    particles.posY()[1] = 0.0;
    particles.posZ()[1] = 0.0;

    const std::size_t stride{particles.paddingStride()};
    for (std::size_t i = 0; i < stride; ++i) {
        particles.gradLogPsiX()[i] = 999.0;
        particles.gradLogPsiY()[i] = 999.0;
        particles.gradLogPsiZ()[i] = 999.0;
        particles.laplogPsi()[i] = 999.0;
    }

    const JastrowPade jastrow{0.5, 1.0};
    SlaterPlaneWave slater{2U, 10.0};
    WaveFunction waveFunction{jastrow, slater};

    std::vector<double> expectedX(stride, 0.0);
    std::vector<double> expectedY(stride, 0.0);
    std::vector<double> expectedZ(stride, 0.0);
    std::vector<double> expectedLap(stride, 0.0);

    jastrow.addDerivatives(
        particles, pbc,
        expectedX.data(),
        expectedY.data(),
        expectedZ.data(),
        expectedLap.data()
    );

    waveFunction.evaluateDerivatives(particles, pbc);

    for (std::size_t i = 0; i < stride; ++i) {
        requireNearWave(particles.gradLogPsiX()[i], expectedX[i]);
        requireNearWave(particles.gradLogPsiY()[i], expectedY[i]);
        requireNearWave(particles.gradLogPsiZ()[i], expectedZ[i]);
        requireNearWave(particles.laplogPsi()[i], expectedLap[i]);
    }
}

TEST_CASE("WaveFunction evaluateLogPsi updates particle logPsi", "[wavefunction]") {
    Particles particles{1U};
    const PeriodicBoundaryCondition pbc{10.0};

    particles.posX()[0] = 0.4;
    particles.posY()[0] = 0.7;
    particles.posZ()[0] = 0.9;

    const JastrowPade jastrow{0.5, 1.0};
    SlaterPlaneWave slater{1U, 10.0};

    slater.kVectorX()[0] = 0.3;
    slater.kVectorY()[0] = -0.2;
    slater.kVectorZ()[0] = 0.5;

    WaveFunction waveFunction{jastrow, slater};
    waveFunction.evaluateLogPsi(particles, pbc);

    const double kdotr{
        slater.kVectorX()[0] * particles.posX()[0] +
        slater.kVectorY()[0] * particles.posY()[0] +
        slater.kVectorZ()[0] * particles.posZ()[0]
    };

    const double expectedLogPsi{
        std::log(std::abs(std::cos(kdotr)))
    };

    requireNearWave(particles.logPsi()[0], expectedLogPsi);
}