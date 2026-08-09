#include "test_utilities.hpp"

#include <catch2/catch_message.hpp>

#include "../src/wavefunction/wavefunction.hpp"

#include <vector>

TEST_CASE("WaveFunction evaluateDerivatives clears buffers and delegates to Jastrow",
          "[wavefunction]") {
  Particles particles{2U};

  particles.pos().x_[0] = 0.0_r;
  particles.pos().y_[0] = 0.0_r;
  particles.pos().z_[0] = 0.0_r;

  particles.pos().x_[1] = 1.0_r;
  particles.pos().y_[1] = 0.0_r;
  particles.pos().z_[1] = 0.0_r;

  const std::size_t stride{particles.p_stride()};
  for (std::size_t i = 0; i < stride; ++i) {
    particles.grad_log_psi().x_[i] = 999.0_r;
    particles.grad_log_psi().y_[i] = 999.0_r;
    particles.grad_log_psi().z_[i] = 999.0_r;
    particles.lap_log_psi()[i] = 999.0_r;
  }

  WaveFunction waveFunction{particles, 10.0_r, 0.5_r, 1.0_r};

  // Compute expected: zero-init then add both Slater and Jastrow
  std::vector<real_t> expectedX(stride, 0.0_r);
  std::vector<real_t> expectedY(stride, 0.0_r);
  std::vector<real_t> expectedZ(stride, 0.0_r);
  std::vector<real_t> expectedLap(stride, 0.0_r);

  // Need to call log_abs_det first to populate the inverse
  waveFunction.slater_plane_wave().log_abs_det(particles);
  waveFunction.slater_plane_wave().add_derivatives(
    expectedX.data(), expectedY.data(),
    expectedZ.data(), expectedLap.data()
  );
  waveFunction.jastrow_pade().add_derivatives(
    particles,
    expectedX.data(), expectedY.data(),
    expectedZ.data(), expectedLap.data()
  );

  waveFunction.evaluate_derivatives(particles);

  for (std::size_t i = 0; i < stride; ++i) {
    require_near(particles.grad_log_psi().x_[i], expectedX[i]);
    require_near(particles.grad_log_psi().y_[i], expectedY[i]);
    require_near(particles.grad_log_psi().z_[i], expectedZ[i]);
    require_near(particles.lap_log_psi()[i], expectedLap[i]);
  }
}

TEST_CASE("WaveFunction evaluate_log_psi returns finite for N=1", "[wavefunction]") {
  // N=1: orbital 0 is cos(0·r) = 1, so log|det| = 0
  // Jastrow with 1 particle has no pairs, so J = 0
  // Total log_psi = 0
  Particles particles{1U};

  particles.pos().x_[0] = 0.4_r;
  particles.pos().y_[0] = 0.7_r;
  particles.pos().z_[0] = 0.9_r;

  WaveFunction waveFunction{particles, 10.0_r, 0.5_r, 1.0_r};

  const real_t logPsi{waveFunction.evaluate_log_psi(particles)};
  require_near(logPsi, 0.0_r);
}

TEST_CASE("WaveFunction default Jastrow parameter matches fully polarized cusp choice",
          "[wavefunction]") {
  Particles particles{2U};
  WaveFunction waveFunction{particles, 9.0_r};

  const real_t actualA{waveFunction.jastrow_pade().a()};
  const real_t actualB{waveFunction.jastrow_pade().b()};

  INFO("Default WaveFunction Jastrow parameters should match the fully polarized HEG choice.");
  CAPTURE(actualA, actualB);
  require_near(actualA, 0.25_r);
  require_near(actualB, 1.0_r);
}