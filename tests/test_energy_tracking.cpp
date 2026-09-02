#include "test_utilities.hpp"

#include <catch2/catch_message.hpp>

#include "energy_tracking/energy_tracking.hpp"

#include <cstddef>

#ifdef FP_64
constexpr fp_t INCREMENTAL_REBUILD_TOLERANCE{1e-9_fp};
#else
constexpr fp_t INCREMENTAL_REBUILD_TOLERANCE{5e-6_fp};
#endif

TEST_CASE("EnergyTracker total energy changes by expected kinetic contribution", "[energy]") {
  constexpr std::size_t n{3U};
  constexpr fp_t L{8.0_fp};
  const EnergyTracker tracker{L, n};

  Particles reference{n};
  reference.pos().x_[0] = 0.3_fp;
  reference.pos().y_[0] = 1.5_fp;
  reference.pos().z_[0] = 2.4_fp;
  reference.pos().x_[1] = 3.7_fp;
  reference.pos().y_[1] = 0.2_fp;
  reference.pos().z_[1] = 5.1_fp;
  reference.pos().x_[2] = 6.9_fp;
  reference.pos().y_[2] = 7.2_fp;
  reference.pos().z_[2] = 1.1_fp;

  Particles with_derivatives{n};
  copy_positions(reference, with_derivatives);

  with_derivatives.grad_log_psi().x_[0] = 0.2_fp;
  with_derivatives.grad_log_psi().y_[0] = -0.1_fp;
  with_derivatives.grad_log_psi().z_[0] = 0.3_fp;
  with_derivatives.lap_log_psi()[0] = 0.5_fp;

  with_derivatives.grad_log_psi().x_[1] = -0.4_fp;
  with_derivatives.grad_log_psi().y_[1] = 0.0_fp;
  with_derivatives.grad_log_psi().z_[1] = 0.1_fp;
  with_derivatives.lap_log_psi()[1] = -0.2_fp;

  with_derivatives.grad_log_psi().x_[2] = 0.3_fp;
  with_derivatives.grad_log_psi().y_[2] = -0.2_fp;
  with_derivatives.grad_log_psi().z_[2] = -0.5_fp;
  with_derivatives.lap_log_psi()[2] = 0.7_fp;

  const fp_t energy_without_derivatives{tracker.eval_total_energy(reference.view())};
  const fp_t energy_with_derivatives{tracker.eval_total_energy(with_derivatives.view())};

  fp_t expected_kinetic{};
  for (std::size_t i = 0; i < n; ++i) {
    const fp_t gx{with_derivatives.grad_log_psi().x_[i]};
    const fp_t gy{with_derivatives.grad_log_psi().y_[i]};
    const fp_t gz{with_derivatives.grad_log_psi().z_[i]};
    const fp_t lap{with_derivatives.lap_log_psi()[i]};
    expected_kinetic += -0.5_fp * (lap + gx * gx + gy * gy + gz * gz);
  }

  require_near(energy_with_derivatives - energy_without_derivatives, expected_kinetic);
}

TEST_CASE("EnergyTracker is invariant under box-periodic particle translations", "[energy]") {
  constexpr std::size_t n{3U};
  constexpr fp_t L{7.5_fp};
  const EnergyTracker tracker{L, n};

  Particles particles{n};
  particles.pos().x_[0] = 1.1_fp;
  particles.pos().y_[0] = 2.3_fp;
  particles.pos().z_[0] = 3.7_fp;
  particles.pos().x_[1] = 6.4_fp;
  particles.pos().y_[1] = 5.6_fp;
  particles.pos().z_[1] = 0.4_fp;
  particles.pos().x_[2] = 2.9_fp;
  particles.pos().y_[2] = 1.2_fp;
  particles.pos().z_[2] = 7.0_fp;

  particles.grad_log_psi().x_[0] = 0.12_fp;
  particles.grad_log_psi().y_[0] = -0.08_fp;
  particles.grad_log_psi().z_[0] = 0.05_fp;
  particles.lap_log_psi()[0] = 0.2_fp;

  particles.grad_log_psi().x_[1] = -0.15_fp;
  particles.grad_log_psi().y_[1] = 0.11_fp;
  particles.grad_log_psi().z_[1] = -0.03_fp;
  particles.lap_log_psi()[1] = -0.1_fp;

  particles.grad_log_psi().x_[2] = 0.07_fp;
  particles.grad_log_psi().y_[2] = 0.09_fp;
  particles.grad_log_psi().z_[2] = -0.04_fp;
  particles.lap_log_psi()[2] = 0.05_fp;

  Particles translated{n};
  copy_positions(particles, translated);
  copy_derivatives(particles, translated);

  translated.pos().x_[0] += L;
  translated.pos().y_[0] -= 2.0_fp * L;
  translated.pos().z_[1] += 3.0_fp * L;
  translated.pos().x_[2] -= L;

  const fp_t baseline{tracker.eval_total_energy(particles.view())};
  const fp_t shifted{tracker.eval_total_energy(translated.view())};
  require_near(shifted, baseline, 1e-9_fp);
}

TEST_CASE("EnergyTracker handles degenerate positions and permutation symmetry", "[energy]") {
  constexpr std::size_t n{2U};
  constexpr fp_t L{6.0_fp};
  const EnergyTracker tracker{L, n};

  Particles particles{n};
  particles.pos().x_[0] = 2.0_fp;
  particles.pos().y_[0] = 2.0_fp;
  particles.pos().z_[0] = 2.0_fp;
  particles.pos().x_[1] = 2.0_fp;
  particles.pos().y_[1] = 2.0_fp;
  particles.pos().z_[1] = 2.0_fp;

  particles.grad_log_psi().x_[0] = 0.1_fp;
  particles.grad_log_psi().y_[0] = 0.2_fp;
  particles.grad_log_psi().z_[0] = 0.3_fp;
  particles.lap_log_psi()[0] = 0.4_fp;
  particles.grad_log_psi().x_[1] = -0.3_fp;
  particles.grad_log_psi().y_[1] = -0.1_fp;
  particles.grad_log_psi().z_[1] = 0.2_fp;
  particles.lap_log_psi()[1] = -0.5_fp;

  const fp_t degenerate_energy{tracker.eval_total_energy(particles.view())};
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

  const fp_t permuted_energy{tracker.eval_total_energy(permuted.view())};
  require_near(permuted_energy, degenerate_energy, 1e-10_fp);
}

TEST_CASE("EnergyTracker incremental reciprocal and real-energy updates match full recomputation",
          "[energy]") {
  constexpr std::size_t n{3U};
  constexpr fp_t L{8.5_fp};
  constexpr std::size_t moved{1U};

  Particles initial{n};
  initial.pos().x_[0] = 0.7_fp;
  initial.pos().y_[0] = 1.1_fp;
  initial.pos().z_[0] = 2.3_fp;
  initial.pos().x_[1] = 3.2_fp;
  initial.pos().y_[1] = 2.7_fp;
  initial.pos().z_[1] = 1.4_fp;
  initial.pos().x_[2] = 5.6_fp;
  initial.pos().y_[2] = 0.9_fp;
  initial.pos().z_[2] = 4.8_fp;

  EnergyTracker incremental{L, n};
  incremental.initialize_structure_factors(initial.view());
  incremental.initialize_reciprocal_energy();
  incremental.initialize_real_energy(initial.view());

  Particles moved_particles{n};
  copy_positions(initial, moved_particles);

  auto pos{moved_particles.pos()};
  const xpu::array<fp_t, idx(Axis::NUM)> old_pos{
    pos[idx(Axis::X)][moved],
    pos[idx(Axis::Y)][moved],
    pos[idx(Axis::Z)][moved]
  };

  pos[idx(Axis::X)][moved] = old_pos[idx(Axis::X)] + 0.21_fp;
  pos[idx(Axis::Y)][moved] = old_pos[idx(Axis::Y)] - 0.17_fp;
  pos[idx(Axis::Z)][moved] = old_pos[idx(Axis::Z)] + 0.33_fp;

  const xpu::array<fp_t, idx(Axis::NUM)> new_pos{
    pos[idx(Axis::X)][moved],
    pos[idx(Axis::Y)][moved],
    pos[idx(Axis::Z)][moved]
  };

  incremental.update_structure_factors(
    old_pos, new_pos
  );
  incremental.update_real_energy(moved, old_pos, moved_particles.view());

  EnergyTracker rebuilt{L, n};
  rebuilt.initialize_structure_factors(moved_particles.view());
  rebuilt.initialize_reciprocal_energy();
  rebuilt.initialize_real_energy(moved_particles.view());

  const fp_t incremental_total{incremental.eval_total_energy(moved_particles.view())};
  const fp_t rebuilt_total{rebuilt.eval_total_energy(moved_particles.view())};

  INFO("Incremental Ewald updates should agree with a full reinitialization after one move.");
  CAPTURE(moved, old_pos, incremental_total, rebuilt_total);
  require_near(incremental_total, rebuilt_total, INCREMENTAL_REBUILD_TOLERANCE);
}
