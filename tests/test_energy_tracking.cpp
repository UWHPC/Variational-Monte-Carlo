#include "test_utilities.hpp"

#include <catch2/catch_message.hpp>

#include "energy_tracking/energy_tracking.hpp"

#include <cstddef>

TEST_CASE("EnergyTracker total energy changes by expected kinetic contribution", "[energy]") {
  constexpr std::size_t n{3U};
  constexpr real_t L{8.0_r};
  const EnergyTracker tracker{L, static_cast<real_t>(n)};

  Particles reference{n};
  reference.pos().x_[0] = 0.3_r;
  reference.pos().y_[0] = 1.5_r;
  reference.pos().z_[0] = 2.4_r;
  reference.pos().x_[1] = 3.7_r;
  reference.pos().y_[1] = 0.2_r;
  reference.pos().z_[1] = 5.1_r;
  reference.pos().x_[2] = 6.9_r;
  reference.pos().y_[2] = 7.2_r;
  reference.pos().z_[2] = 1.1_r;

  Particles with_derivatives{n};
  copy_positions(reference, with_derivatives);

  with_derivatives.grad_log_psi().x_[0] = 0.2_r;
  with_derivatives.grad_log_psi().y_[0] = -0.1_r;
  with_derivatives.grad_log_psi().z_[0] = 0.3_r;
  with_derivatives.lap_log_psi()[0] = 0.5_r;

  with_derivatives.grad_log_psi().x_[1] = -0.4_r;
  with_derivatives.grad_log_psi().y_[1] = 0.0_r;
  with_derivatives.grad_log_psi().z_[1] = 0.1_r;
  with_derivatives.lap_log_psi()[1] = -0.2_r;

  with_derivatives.grad_log_psi().x_[2] = 0.3_r;
  with_derivatives.grad_log_psi().y_[2] = -0.2_r;
  with_derivatives.grad_log_psi().z_[2] = -0.5_r;
  with_derivatives.lap_log_psi()[2] = 0.7_r;

  const real_t energy_without_derivatives{tracker.eval_total_energy(reference)};
  const real_t energy_with_derivatives{tracker.eval_total_energy(with_derivatives)};

  real_t expected_kinetic{};
  for (std::size_t i = 0; i < n; ++i) {
    const real_t gx{with_derivatives.grad_log_psi().x_[i]};
    const real_t gy{with_derivatives.grad_log_psi().y_[i]};
    const real_t gz{with_derivatives.grad_log_psi().z_[i]};
    const real_t lap{with_derivatives.lap_log_psi()[i]};
    expected_kinetic += -0.5_r * (lap + gx * gx + gy * gy + gz * gz);
  }

  require_near(energy_with_derivatives - energy_without_derivatives, expected_kinetic);
}

TEST_CASE("EnergyTracker is invariant under box-periodic particle translations", "[energy]") {
  constexpr std::size_t n{3U};
  constexpr real_t L{7.5_r};
  const EnergyTracker tracker{L, static_cast<real_t>(n)};

  Particles particles{n};
  particles.pos().x_[0] = 1.1_r;
  particles.pos().y_[0] = 2.3_r;
  particles.pos().z_[0] = 3.7_r;
  particles.pos().x_[1] = 6.4_r;
  particles.pos().y_[1] = 5.6_r;
  particles.pos().z_[1] = 0.4_r;
  particles.pos().x_[2] = 2.9_r;
  particles.pos().y_[2] = 1.2_r;
  particles.pos().z_[2] = 7.0_r;

  particles.grad_log_psi().x_[0] = 0.12_r;
  particles.grad_log_psi().y_[0] = -0.08_r;
  particles.grad_log_psi().z_[0] = 0.05_r;
  particles.lap_log_psi()[0] = 0.2_r;

  particles.grad_log_psi().x_[1] = -0.15_r;
  particles.grad_log_psi().y_[1] = 0.11_r;
  particles.grad_log_psi().z_[1] = -0.03_r;
  particles.lap_log_psi()[1] = -0.1_r;

  particles.grad_log_psi().x_[2] = 0.07_r;
  particles.grad_log_psi().y_[2] = 0.09_r;
  particles.grad_log_psi().z_[2] = -0.04_r;
  particles.lap_log_psi()[2] = 0.05_r;

  Particles translated{n};
  copy_positions(particles, translated);
  copy_derivatives(particles, translated);

  translated.pos().x_[0] += L;
  translated.pos().y_[0] -= 2.0_r * L;
  translated.pos().z_[1] += 3.0_r * L;
  translated.pos().x_[2] -= L;

  const real_t baseline{tracker.eval_total_energy(particles)};
  const real_t shifted{tracker.eval_total_energy(translated)};
  require_near(shifted, baseline, 1e-9_r);
}

TEST_CASE("EnergyTracker handles degenerate positions and permutation symmetry", "[energy]") {
  constexpr std::size_t n{2U};
  constexpr real_t L{6.0_r};
  const EnergyTracker tracker{L, static_cast<real_t>(n)};

  Particles particles{n};
  particles.pos().x_[0] = 2.0_r;
  particles.pos().y_[0] = 2.0_r;
  particles.pos().z_[0] = 2.0_r;
  particles.pos().x_[1] = 2.0_r;
  particles.pos().y_[1] = 2.0_r;
  particles.pos().z_[1] = 2.0_r;

  particles.grad_log_psi().x_[0] = 0.1_r;
  particles.grad_log_psi().y_[0] = 0.2_r;
  particles.grad_log_psi().z_[0] = 0.3_r;
  particles.lap_log_psi()[0] = 0.4_r;
  particles.grad_log_psi().x_[1] = -0.3_r;
  particles.grad_log_psi().y_[1] = -0.1_r;
  particles.grad_log_psi().z_[1] = 0.2_r;
  particles.lap_log_psi()[1] = -0.5_r;

  const real_t degenerate_energy{tracker.eval_total_energy(particles)};
  REQUIRE(std::isfinite(degenerate_energy));

  Particles permuted{n};
  permuted.pos().x_[0] = particles.pos().x_[1];
  permuted.pos().y_[0] = particles.pos().y_[1];
  permuted.pos().z_[0] = particles.pos().z_[1];
  permuted.pos().x_[1] = particles.pos().x_[0];
  permuted.pos().y_[1] = particles.pos().y_[0];
  permuted.pos().z_[1] = particles.pos().z_[0];

  permuted.grad_log_psi().x_[0] = particles.grad_log_psi().x_[1];
  permuted.grad_log_psi().y_[0] = particles.grad_log_psi().y_[1];
  permuted.grad_log_psi().z_[0] = particles.grad_log_psi().z_[1];
  permuted.lap_log_psi()[0] = particles.lap_log_psi()[1];
  permuted.grad_log_psi().x_[1] = particles.grad_log_psi().x_[0];
  permuted.grad_log_psi().y_[1] = particles.grad_log_psi().y_[0];
  permuted.grad_log_psi().z_[1] = particles.grad_log_psi().z_[0];
  permuted.lap_log_psi()[1] = particles.lap_log_psi()[0];

  const real_t permuted_energy{tracker.eval_total_energy(permuted)};
  require_near(permuted_energy, degenerate_energy, 1e-10_r);
}

TEST_CASE("EnergyTracker incremental reciprocal and real-energy updates match full recomputation",
          "[energy]") {
  constexpr std::size_t n{3U};
  constexpr real_t L{8.5_r};
  constexpr std::size_t moved{1U};

  Particles initial{n};
  initial.pos().x_[0] = 0.7_r;
  initial.pos().y_[0] = 1.1_r;
  initial.pos().z_[0] = 2.3_r;
  initial.pos().x_[1] = 3.2_r;
  initial.pos().y_[1] = 2.7_r;
  initial.pos().z_[1] = 1.4_r;
  initial.pos().x_[2] = 5.6_r;
  initial.pos().y_[2] = 0.9_r;
  initial.pos().z_[2] = 4.8_r;

  EnergyTracker incremental{L, static_cast<real_t>(n)};
  incremental.initialize_structure_factors(initial);
  incremental.initialize_reciprocal_energy();
  incremental.initialize_real_energy(initial);

  Particles moved_particles{n};
  copy_positions(initial, moved_particles);

  const real_t old_x{moved_particles.pos().x_[moved]};
  const real_t old_y{moved_particles.pos().y_[moved]};
  const real_t old_z{moved_particles.pos().z_[moved]};

  moved_particles.pos().x_[moved] = old_x + 0.21_r;
  moved_particles.pos().y_[moved] = old_y - 0.17_r;
  moved_particles.pos().z_[moved] = old_z + 0.33_r;

  incremental.update_structure_factors(
    old_x, old_y, old_z,
    moved_particles.pos().x_[moved],
    moved_particles.pos().y_[moved],
    moved_particles.pos().z_[moved]
  );
  incremental.update_real_energy(moved, old_x, old_y, old_z, moved_particles);

  EnergyTracker rebuilt{L, static_cast<real_t>(n)};
  rebuilt.initialize_structure_factors(moved_particles);
  rebuilt.initialize_reciprocal_energy();
  rebuilt.initialize_real_energy(moved_particles);

  const real_t incremental_total{incremental.eval_total_energy(moved_particles)};
  const real_t rebuilt_total{rebuilt.eval_total_energy(moved_particles)};

  INFO("Incremental Ewald updates should agree with a full reinitialization after one move.");
  CAPTURE(moved, old_x, old_y, old_z, incremental_total, rebuilt_total);
  require_near(incremental_total, rebuilt_total, 1e-9_r);
}