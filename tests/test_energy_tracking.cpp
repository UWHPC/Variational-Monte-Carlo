#include "test_utilities.hpp"

#include <catch2/catch_message.hpp>

#include "energy_tracking/energy_tracking.hpp"

#include <cstddef>

TEST_CASE("EnergyTracker total energy changes by expected kinetic contribution", "[energy]") {
  constexpr std::size_t n{3U};
  constexpr double      L{8.0};
  const EnergyTracker   tracker{L, static_cast<double>(n)};

  Particles reference{n};
  reference.pos().x_[0] = 0.3;
  reference.pos().y_[0] = 1.5;
  reference.pos().z_[0] = 2.4;
  reference.pos().x_[1] = 3.7;
  reference.pos().y_[1] = 0.2;
  reference.pos().z_[1] = 5.1;
  reference.pos().x_[2] = 6.9;
  reference.pos().y_[2] = 7.2;
  reference.pos().z_[2] = 1.1;

  Particles with_derivatives{n};
  copy_positions(reference, with_derivatives);

  with_derivatives.grad_log_psi().x_[0] = 0.2;
  with_derivatives.grad_log_psi().y_[0] = -0.1;
  with_derivatives.grad_log_psi().z_[0] = 0.3;
  with_derivatives.lap_log_psi()[0] = 0.5;

  with_derivatives.grad_log_psi().x_[1] = -0.4;
  with_derivatives.grad_log_psi().y_[1] = 0.0;
  with_derivatives.grad_log_psi().z_[1] = 0.1;
  with_derivatives.lap_log_psi()[1] = -0.2;

  with_derivatives.grad_log_psi().x_[2] = 0.3;
  with_derivatives.grad_log_psi().y_[2] = -0.2;
  with_derivatives.grad_log_psi().z_[2] = -0.5;
  with_derivatives.lap_log_psi()[2] = 0.7;

  const double energy_without_derivatives{tracker.eval_total_energy(reference)};
  const double energy_with_derivatives{tracker.eval_total_energy(with_derivatives)};

  double expected_kinetic{};
  for (std::size_t i = 0; i < n; ++i) {
    const double gx{with_derivatives.grad_log_psi().x_[i]};
    const double gy{with_derivatives.grad_log_psi().y_[i]};
    const double gz{with_derivatives.grad_log_psi().z_[i]};
    const double lap{with_derivatives.lap_log_psi()[i]};
    expected_kinetic += -0.5 * (lap + gx * gx + gy * gy + gz * gz);
  }

  require_near(energy_with_derivatives - energy_without_derivatives, expected_kinetic);
}

TEST_CASE("EnergyTracker is invariant under box-periodic particle translations", "[energy]") {
  constexpr std::size_t n{3U};
  constexpr double      L{7.5};
  const EnergyTracker   tracker{L, static_cast<double>(n)};

  Particles particles{n};
  particles.pos().x_[0] = 1.1;
  particles.pos().y_[0] = 2.3;
  particles.pos().z_[0] = 3.7;
  particles.pos().x_[1] = 6.4;
  particles.pos().y_[1] = 5.6;
  particles.pos().z_[1] = 0.4;
  particles.pos().x_[2] = 2.9;
  particles.pos().y_[2] = 1.2;
  particles.pos().z_[2] = 7.0;

  particles.grad_log_psi().x_[0] = 0.12;
  particles.grad_log_psi().y_[0] = -0.08;
  particles.grad_log_psi().z_[0] = 0.05;
  particles.lap_log_psi()[0] = 0.2;

  particles.grad_log_psi().x_[1] = -0.15;
  particles.grad_log_psi().y_[1] = 0.11;
  particles.grad_log_psi().z_[1] = -0.03;
  particles.lap_log_psi()[1] = -0.1;

  particles.grad_log_psi().x_[2] = 0.07;
  particles.grad_log_psi().y_[2] = 0.09;
  particles.grad_log_psi().z_[2] = -0.04;
  particles.lap_log_psi()[2] = 0.05;

  Particles translated{n};
  copy_positions(particles, translated);
  copy_derivatives(particles, translated);

  translated.pos().x_[0] += L;
  translated.pos().y_[0] -= 2.0 * L;
  translated.pos().z_[1] += 3.0 * L;
  translated.pos().x_[2] -= L;

  const double baseline{tracker.eval_total_energy(particles)};
  const double shifted{tracker.eval_total_energy(translated)};
  require_near(shifted, baseline, 1e-9);
}

TEST_CASE("EnergyTracker handles degenerate positions and permutation symmetry", "[energy]") {
  constexpr std::size_t n{2U};
  constexpr double      L{6.0};
  const EnergyTracker   tracker{L, static_cast<double>(n)};

  Particles particles{n};
  particles.pos().x_[0] = 2.0;
  particles.pos().y_[0] = 2.0;
  particles.pos().z_[0] = 2.0;
  particles.pos().x_[1] = 2.0;
  particles.pos().y_[1] = 2.0;
  particles.pos().z_[1] = 2.0;

  particles.grad_log_psi().x_[0] = 0.1;
  particles.grad_log_psi().y_[0] = 0.2;
  particles.grad_log_psi().z_[0] = 0.3;
  particles.lap_log_psi()[0] = 0.4;
  particles.grad_log_psi().x_[1] = -0.3;
  particles.grad_log_psi().y_[1] = -0.1;
  particles.grad_log_psi().z_[1] = 0.2;
  particles.lap_log_psi()[1] = -0.5;

  const double degenerate_energy{tracker.eval_total_energy(particles)};
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

  const double permuted_energy{tracker.eval_total_energy(permuted)};
  require_near(permuted_energy, degenerate_energy, 1e-10);
}

TEST_CASE("EnergyTracker incremental reciprocal and real-energy updates match full recomputation",
          "[energy]") {
  constexpr std::size_t n{3U};
  constexpr double      L{8.5};
  constexpr std::size_t moved{1U};

  Particles initial{n};
  initial.pos().x_[0] = 0.7;
  initial.pos().y_[0] = 1.1;
  initial.pos().z_[0] = 2.3;
  initial.pos().x_[1] = 3.2;
  initial.pos().y_[1] = 2.7;
  initial.pos().z_[1] = 1.4;
  initial.pos().x_[2] = 5.6;
  initial.pos().y_[2] = 0.9;
  initial.pos().z_[2] = 4.8;

  EnergyTracker incremental{L, static_cast<double>(n)};
  incremental.initialize_structure_factors(initial);
  incremental.initialize_reciprocal_energy();
  incremental.initialize_real_energy(initial);

  Particles moved_particles{n};
  copy_positions(initial, moved_particles);

  const double old_x{moved_particles.pos().x_[moved]};
  const double old_y{moved_particles.pos().y_[moved]};
  const double old_z{moved_particles.pos().z_[moved]};

  moved_particles.pos().x_[moved] = old_x + 0.21;
  moved_particles.pos().y_[moved] = old_y - 0.17;
  moved_particles.pos().z_[moved] = old_z + 0.33;

  incremental.update_structure_factors(
    old_x, old_y, old_z,
    moved_particles.pos().x_[moved],
    moved_particles.pos().y_[moved],
    moved_particles.pos().z_[moved]
  );
  incremental.update_real_energy(moved, old_x, old_y, old_z, moved_particles);

  EnergyTracker rebuilt{L, static_cast<double>(n)};
  rebuilt.initialize_structure_factors(moved_particles);
  rebuilt.initialize_reciprocal_energy();
  rebuilt.initialize_real_energy(moved_particles);

  const double incremental_total{incremental.eval_total_energy(moved_particles)};
  const double rebuilt_total{rebuilt.eval_total_energy(moved_particles)};

  INFO("Incremental Ewald updates should agree with a full reinitialization after one move.");
  CAPTURE(moved, old_x, old_y, old_z, incremental_total, rebuilt_total);
  require_near(incremental_total, rebuilt_total, 1e-9);
}