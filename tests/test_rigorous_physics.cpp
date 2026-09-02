#include "test_utilities.hpp"

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include "energy_tracking/energy_tracking.hpp"
#include "jastrow_pade/jastrow_pade.hpp"
#include "particles/particles.hpp"
#include "slater_plane_wave/slater_plane_wave.hpp"
#include "wavefunction/wavefunction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <random>
#include <vector>

namespace {

#ifdef FP_64
constexpr fp_t RIGOROUS_PRECISION_SCALE{1.0_fp};
constexpr std::size_t UNREFRESHED_DRIFT_STEPS{64U};
constexpr fp_t MAINTAINED_RESIDUAL_TOLERANCE{1e-8_fp};
constexpr fp_t FRESH_RESIDUAL_TOLERANCE{1e-10_fp};
constexpr fp_t PROBE_RATIO_TOLERANCE{1e-8_fp};
constexpr fp_t TRANSLATION_LOG_PSI_TOLERANCE{1e-12_fp};
constexpr fp_t TRANSLATION_GRAD_TOLERANCE{1e-10_fp};
constexpr fp_t TRANSLATION_LAP_TOLERANCE{1e-7_fp};
#else
constexpr fp_t RIGOROUS_PRECISION_SCALE{1e6_fp};
constexpr std::size_t UNREFRESHED_DRIFT_STEPS{16U};
constexpr fp_t MAINTAINED_RESIDUAL_TOLERANCE{1e-1_fp};
constexpr fp_t FRESH_RESIDUAL_TOLERANCE{1e-2_fp};
constexpr fp_t PROBE_RATIO_TOLERANCE{1e-2_fp};
constexpr fp_t TRANSLATION_LOG_PSI_TOLERANCE{1e-4_fp};
constexpr fp_t TRANSLATION_GRAD_TOLERANCE{1e-2_fp};
constexpr fp_t TRANSLATION_LAP_TOLERANCE{1e-1_fp};
#endif

constexpr fp_t EWALD_RECIPROCAL_TOLERANCE{1.0e-6_fp};

fp_t phy_determinant3x3(const fp_t* matrix, std::size_t stride) {
  return
    matrix[0 * stride + 0] * (
      matrix[1 * stride + 1] * matrix[2 * stride + 2] -
      matrix[1 * stride + 2] * matrix[2 * stride + 1]
    ) -
    matrix[0 * stride + 1] * (
      matrix[1 * stride + 0] * matrix[2 * stride + 2] -
      matrix[1 * stride + 2] * matrix[2 * stride + 0]
    ) +
    matrix[0 * stride + 2] * (
      matrix[1 * stride + 0] * matrix[2 * stride + 1] -
      matrix[1 * stride + 1] * matrix[2 * stride + 0]
    );
}

void phy_copy_positions(const Particles& source, Particles& dest) {
  const std::size_t n{source.size()};

  for (std::size_t i = 0; i < n; ++i) {
    dest.pos().x_[i] = source.pos().x_[i];
    dest.pos().y_[i] = source.pos().y_[i];
    dest.pos().z_[i] = source.pos().z_[i];
  }
}

fp_t exact_real_potential(const Particles& particles, fp_t box_length) {
  const std::size_t n{particles.size()};
  const fp_t alpha{6.0_fp / box_length};

  fp_t total{};

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      const fp_t dx{minimum_image(particles.pos().x_[i] - particles.pos().x_[j], box_length)};
      const fp_t dy{minimum_image(particles.pos().y_[i] - particles.pos().y_[j], box_length)};
      const fp_t dz{minimum_image(particles.pos().z_[i] - particles.pos().z_[j], box_length)};

      const fp_t r{std::sqrt(dx * dx + dy * dy + dz * dz)};
      total += std::erfc(alpha * r) / r;
    }
  }

  return total;
}

fp_t exact_reciprocal_potential(const Particles& particles, fp_t box_length) {
  const std::size_t n{particles.size()};
  const fp_t alpha{6.0_fp / box_length};
  const fp_t two_pi_over_l{2.0_fp * std::numbers::pi_v<fp_t> / box_length};
  const fp_t four_alpha_sq{4.0_fp * alpha * alpha};
  const fp_t cutoff_factor{-std::log(EWALD_RECIPROCAL_TOLERANCE)};
  const fp_t g_max_mag_sq{four_alpha_sq * cutoff_factor};
  const int m_max{static_cast<int>(std::ceil(std::sqrt(g_max_mag_sq) / two_pi_over_l)) + 1};

  fp_t weighted_sum{};

  for (int mx = -m_max; mx <= m_max; ++mx) {
    for (int my = -m_max; my <= m_max; ++my) {
      for (int mz = -m_max; mz <= m_max; ++mz) {
        if (mx < 0) {
          continue;
        }
        if (mx == 0 && my < 0) {
          continue;
        }
        if (mx == 0 && my == 0 && mz <= 0) {
          continue;
        }

        const fp_t gx{two_pi_over_l * static_cast<fp_t>(mx)};
        const fp_t gy{two_pi_over_l * static_cast<fp_t>(my)};
        const fp_t gz{two_pi_over_l * static_cast<fp_t>(mz)};
        const fp_t g_sq{gx * gx + gy * gy + gz * gz};

        if (g_sq > g_max_mag_sq) {
          continue;
        }

        const fp_t weight{
          8.0_fp * std::numbers::pi_v<fp_t> * std::numbers::pi_v<fp_t> *
          std::exp(-g_sq / four_alpha_sq) / g_sq
        };

        fp_t s_real{};
        fp_t s_imag{};

        for (std::size_t i = 0; i < n; ++i) {
          const fp_t phase{
            gx * particles.pos().x_[i] +
            gy * particles.pos().y_[i] +
            gz * particles.pos().z_[i]
          };
          s_real += std::cos(phase);
          s_imag += std::sin(phase);
        }

        weighted_sum += weight * (s_real * s_real + s_imag * s_imag);
      }
    }
  }

  return weighted_sum / (2.0_fp * std::numbers::pi_v<fp_t> * box_length * box_length * box_length);
}

fp_t exact_total_potential(const Particles& particles, fp_t box_length) {
  const fp_t n{static_cast<fp_t>(particles.size())};
  const fp_t self_correction{-6.0_fp * n / (std::sqrt(std::numbers::pi_v<fp_t>) * box_length)};
  const fp_t background{-std::numbers::pi_v<fp_t> * n * n / (72.0_fp * box_length)};

  return
    exact_real_potential(particles, box_length) +
    exact_reciprocal_potential(particles, box_length) +
    self_correction +
    background;
}

} // namespace

TEST_CASE("Fully polarized Jastrow cusp matches the same-spin analytical form near contact",
          "[rigorous][jastrow]") {
  constexpr fp_t box_length{50.0_fp};
  const JastrowPade jastrow{box_length, 0.25_fp, 1.0_fp};
  Particles particles{2U};

  particles.pos().x_[0] = 0.0_fp;
  particles.pos().y_[0] = 0.0_fp;
  particles.pos().z_[0] = 0.0_fp;

  for (const fp_t r : {1.0e-5_fp, 2.0e-5_fp, 5.0e-5_fp, 1.0e-4_fp}) {
    particles.pos().x_[1] = r;
    particles.pos().y_[1] = 0.0_fp;
    particles.pos().z_[1] = 0.0_fp;

    const fp_t expected_value{0.25_fp * r / (1.0_fp + r)};
    const fp_t actual_value{jastrow.value(particles.view())};

    INFO("Checking Jastrow value against the exact same-spin Padé form near coalescence.");
    CAPTURE(r, actual_value, expected_value);
    REQUIRE(std::abs(actual_value - expected_value) <= (1e-14_fp * RIGOROUS_PRECISION_SCALE));

    std::vector<fp_t> grad_x(particles.p_stride(), 0.0_fp);
    std::vector<fp_t> grad_y(particles.p_stride(), 0.0_fp);
    std::vector<fp_t> grad_z(particles.p_stride(), 0.0_fp);
    std::vector<fp_t> lap(particles.p_stride(), 0.0_fp);

    jastrow.add_derivatives(particles.view(), grad_x.data(), grad_y.data(), grad_z.data(), lap.data());

    const fp_t first_derivative{0.25_fp / ((1.0_fp + r) * (1.0_fp + r))};
    const fp_t second_derivative{-0.5_fp / ((1.0_fp + r) * (1.0_fp + r) * (1.0_fp + r))};
    const fp_t expected_grad_particle_0{-first_derivative};
    const fp_t expected_grad_particle_1{first_derivative};
    const fp_t expected_laplacian{second_derivative + 2.0_fp * first_derivative / r};

    INFO("Checking same-spin Jastrow gradient and Laplacian near contact.");
    CAPTURE(r, first_derivative, second_derivative, expected_laplacian);
    REQUIRE(std::abs(grad_x[0] - expected_grad_particle_0) <= (1e-10_fp * RIGOROUS_PRECISION_SCALE));
    REQUIRE(std::abs(grad_x[1] - expected_grad_particle_1) <= (1e-10_fp * RIGOROUS_PRECISION_SCALE));
    REQUIRE(std::abs(grad_y[0]) <= (1e-14_fp * RIGOROUS_PRECISION_SCALE));
    REQUIRE(std::abs(grad_y[1]) <= (1e-14_fp * RIGOROUS_PRECISION_SCALE));
    REQUIRE(std::abs(grad_z[0]) <= (1e-14_fp * RIGOROUS_PRECISION_SCALE));
    REQUIRE(std::abs(grad_z[1]) <= (1e-14_fp * RIGOROUS_PRECISION_SCALE));
    REQUIRE(std::abs(lap[0] - expected_laplacian) <= (1e-7_fp * RIGOROUS_PRECISION_SCALE));
    REQUIRE(std::abs(lap[1] - expected_laplacian) <= (1e-7_fp * RIGOROUS_PRECISION_SCALE));
  }
}

TEST_CASE(
    "Slater determinant changes sign under particle exchange while log_abs_det stays invariant",
    "[rigorous][slater]") {
  constexpr std::size_t n{3U};
  constexpr fp_t box_length{9.0_fp};

  Particles original{n};
  original.pos().x_[0] = 0.7_fp;
  original.pos().y_[0] = 1.4_fp;
  original.pos().z_[0] = 2.1_fp;

  original.pos().x_[1] = 2.8_fp;
  original.pos().y_[1] = 0.9_fp;
  original.pos().z_[1] = 1.6_fp;

  original.pos().x_[2] = 4.3_fp;
  original.pos().y_[2] = 3.2_fp;
  original.pos().z_[2] = 0.5_fp;

  SlaterPlaneWave slater_original{original, box_length};
  const fp_t log_original{slater_original.log_abs_det(original.view())};
  REQUIRE(std::isfinite(log_original));

  const fp_t det_original{
    phy_determinant3x3(
      slater_original.determinant(),
      slater_original.matrix_row_stride()
    )
  };

  Particles swapped{n};
  phy_copy_positions(original, swapped);
  std::swap(swapped.pos().x_[0], swapped.pos().x_[1]);
  std::swap(swapped.pos().y_[0], swapped.pos().y_[1]);
  std::swap(swapped.pos().z_[0], swapped.pos().z_[1]);

  SlaterPlaneWave slater_swapped{swapped, box_length};
  const fp_t log_swapped{slater_swapped.log_abs_det(swapped.view())};
  REQUIRE(std::isfinite(log_swapped));

  const fp_t det_swapped{
    phy_determinant3x3(
      slater_swapped.determinant(),
      slater_swapped.matrix_row_stride()
    )
  };

  INFO("Swapping two particles must negate the Slater determinant but preserve log|det|.");
  CAPTURE(det_original, det_swapped, log_original, log_swapped);
  REQUIRE(std::abs(det_swapped + det_original) <= (1e-10_fp * RIGOROUS_PRECISION_SCALE));
  REQUIRE(std::abs(log_swapped - log_original) <= (1e-12_fp * RIGOROUS_PRECISION_SCALE));
}

TEST_CASE("SlaterPlaneWave is periodic under box translations of a single particle",
          "[rigorous][slater]") {
  constexpr std::size_t n{3U};
  constexpr fp_t box_length{10.0_fp};

  Particles particles{n};
  particles.pos().x_[0] = 1.1_fp;
  particles.pos().y_[0] = 2.2_fp;
  particles.pos().z_[0] = 3.3_fp;

  particles.pos().x_[1] = 4.4_fp;
  particles.pos().y_[1] = 5.5_fp;
  particles.pos().z_[1] = 6.6_fp;

  particles.pos().x_[2] = 7.7_fp;
  particles.pos().y_[2] = 1.8_fp;
  particles.pos().z_[2] = 2.9_fp;

  SlaterPlaneWave baseline{particles, box_length};
  const fp_t baseline_log_det{baseline.log_abs_det(particles.view())};
  REQUIRE(std::isfinite(baseline_log_det));

  Particles shifted{n};
  phy_copy_positions(particles, shifted);
  shifted.pos().x_[1] += box_length;
  shifted.pos().y_[1] -= box_length;
  shifted.pos().z_[1] += 2.0_fp * box_length;

  SlaterPlaneWave translated{shifted, box_length};
  const fp_t translated_log_det{translated.log_abs_det(shifted.view())};
  REQUIRE(std::isfinite(translated_log_det));

  INFO("Plane-wave Slater determinant must be invariant under integer box translations.");
  CAPTURE(baseline_log_det, translated_log_det);
  REQUIRE(std::abs(translated_log_det - baseline_log_det) <= (1e-12_fp * RIGOROUS_PRECISION_SCALE));

  const std::size_t S{baseline.matrix_row_stride()};
  for (std::size_t row = 0; row < n; ++row) {
    for (std::size_t col = 0; col < n; ++col) {
      const std::size_t idx{row * S + col};
      CAPTURE(row, col);
      REQUIRE(std::abs(baseline.determinant()[idx] - translated.determinant()[idx]) <= (1e-12_fp * RIGOROUS_PRECISION_SCALE));
    }
  }
}

TEST_CASE("Randomized determinant_ratio matches exact determinant ratio over many configurations",
          "[rigorous][slater]") {
  constexpr std::size_t n{3U};
  constexpr fp_t box_length{11.0_fp};
  constexpr std::size_t samples{40U};

  std::mt19937_64 rng{123456789ULL};
  std::uniform_real_distribution<fp_t> position_dist{0.0_fp, box_length};
  std::uniform_real_distribution<fp_t> move_dist{-0.20_fp, 0.20_fp};

  for (std::size_t sample = 0; sample < samples; ++sample) {
    Particles particles{n};

    for (std::size_t i = 0; i < n; ++i) {
      particles.pos().x_[i] = position_dist(rng);
      particles.pos().y_[i] = position_dist(rng);
      particles.pos().z_[i] = position_dist(rng);
    }

    SlaterPlaneWave slater{particles, box_length};
    const fp_t baseline_log_det{slater.log_abs_det(particles.view())};

    INFO("Baseline Slater rebuild must be finite for randomized determinant-ratio testing.");
    CAPTURE(sample, baseline_log_det);
    REQUIRE(std::isfinite(baseline_log_det));

    const fp_t det_old{
      phy_determinant3x3(
        slater.determinant(),
        slater.matrix_row_stride()
      )
    };
    REQUIRE(std::abs(det_old) > (1e-12_fp * RIGOROUS_PRECISION_SCALE));

    const std::size_t moved{sample % n};
    particles.pos().x_[moved] =
        wrap_coordinate(particles.pos().x_[moved] + move_dist(rng), box_length);
    particles.pos().y_[moved] =
        wrap_coordinate(particles.pos().y_[moved] + move_dist(rng), box_length);
    particles.pos().z_[moved] =
        wrap_coordinate(particles.pos().z_[moved] + move_dist(rng), box_length);

    slater.update_trig_cache(moved, particles.view());
    const fp_t* const new_row{slater.build_row(moved)};
    const fp_t ratio{slater.determinant_ratio(moved, new_row)};

    SlaterPlaneWave rebuilt{particles, box_length};
    const fp_t rebuilt_log_det{rebuilt.log_abs_det(particles.view())};

    INFO("Fresh rebuild must be finite for randomized determinant-ratio testing.");
    CAPTURE(sample, moved, rebuilt_log_det);
    REQUIRE(std::isfinite(rebuilt_log_det));

    const fp_t det_new{
      phy_determinant3x3(
        rebuilt.determinant(),
        rebuilt.matrix_row_stride()
      )
    };
    const fp_t exact_ratio{det_new / det_old};

    INFO("Fast determinant ratio must agree with exact determinant ratio over many random "
         "configurations.");
    CAPTURE(sample, moved, ratio, exact_ratio, det_old, det_new);
    REQUIRE(std::abs(ratio - exact_ratio) <= (1e-9_fp * RIGOROUS_PRECISION_SCALE));
  }
}

TEST_CASE("SlaterPlaneWave reports a singular determinant for duplicate particle rows",
          "[rigorous][slater]") {
  constexpr std::size_t n{3U};
  constexpr fp_t box_length{10.0_fp};

  Particles particles{n};
  particles.pos().x_[0] = 1.25_fp;
  particles.pos().y_[0] = 2.50_fp;
  particles.pos().z_[0] = 3.75_fp;

  particles.pos().x_[1] = particles.pos().x_[0];
  particles.pos().y_[1] = particles.pos().y_[0];
  particles.pos().z_[1] = particles.pos().z_[0];

  particles.pos().x_[2] = 4.10_fp;
  particles.pos().y_[2] = 0.90_fp;
  particles.pos().z_[2] = 1.70_fp;

  SlaterPlaneWave slater{particles, box_length};
  const fp_t log_abs_det{slater.log_abs_det(particles.view())};
  const fp_t det{phy_determinant3x3(slater.determinant(), slater.matrix_row_stride())};

  INFO("Duplicate particle positions should create duplicate Slater rows and a singular matrix.");
  CAPTURE(log_abs_det, det);
  REQUIRE_FALSE(std::isfinite(log_abs_det));
  REQUIRE(std::abs(det) <= (1e-12_fp * RIGOROUS_PRECISION_SCALE));
}

TEST_CASE("Repeated accepted Slater updates preserve predictive determinant ratios against fresh "
          "rebuilds",
          "[rigorous][slater]") {
  constexpr std::size_t n{7U};
  constexpr fp_t box_length{10.0_fp};
  constexpr std::size_t steps{UNREFRESHED_DRIFT_STEPS};

  Particles particles{n};
  set_stable_closed_shell_positions(particles);

  SlaterPlaneWave maintained{particles, box_length};
  const fp_t initial_log_det{maintained.log_abs_det(particles.view())};

  INFO("Initial Slater matrix for the multi-step drift test must be well-defined.");
  CAPTURE(initial_log_det);
  REQUIRE(std::isfinite(initial_log_det));

  for (std::size_t step = 0; step < steps; ++step) {
    const std::size_t moved{step % n};

    const fp_t dx{0.003_fp * static_cast<fp_t>((step % 5U) + 1U)};
    const fp_t dy{-0.0025_fp * static_cast<fp_t>((step % 7U) + 1U)};
    const fp_t dz{0.002_fp * static_cast<fp_t>((step % 3U) + 1U)};

    particles.pos().x_[moved] = wrap_coordinate(particles.pos().x_[moved] + dx, box_length);
    particles.pos().y_[moved] = wrap_coordinate(particles.pos().y_[moved] + dy, box_length);
    particles.pos().z_[moved] = wrap_coordinate(particles.pos().z_[moved] + dz, box_length);

    maintained.update_trig_cache(moved, particles.view());
    const fp_t* const accepted_row{maintained.build_row(moved)};
    const fp_t accepted_ratio{maintained.determinant_ratio(moved, accepted_row)};

    INFO("Accepted-step determinant ratio became non-finite or too close to singular.");
    CAPTURE(step, moved, accepted_ratio);
    REQUIRE(std::isfinite(accepted_ratio));
    REQUIRE(std::abs(accepted_ratio) > (1e-10_fp * RIGOROUS_PRECISION_SCALE));

    maintained.accept_move(moved, accepted_row, accepted_ratio);

    SlaterPlaneWave rebuilt{particles, box_length};
    const fp_t rebuilt_log_det{rebuilt.log_abs_det(particles.view())};

    INFO("Fresh rebuild produced a non-finite log determinant during the multi-step "
         "predictive-consistency test.");
    CAPTURE(step, moved, rebuilt_log_det);
    REQUIRE(std::isfinite(rebuilt_log_det));

    const fp_t maintained_residual{slater_identity_residual(maintained)};
    const fp_t rebuilt_residual{slater_identity_residual(rebuilt)};

    INFO("Both maintained and rebuilt Slater states must remain valid inverses of their "
         "determinant matrices.");
    CAPTURE(step, moved, maintained_residual, rebuilt_residual);
    REQUIRE(maintained_residual <= MAINTAINED_RESIDUAL_TOLERANCE);
    REQUIRE(rebuilt_residual <= FRESH_RESIDUAL_TOLERANCE);

    fp_t max_det_diff{};
    const std::size_t S_MAT{maintained.matrix_row_stride()};
    for (std::size_t row = 0; row < n; ++row) {
      for (std::size_t col = 0; col < n; ++col) {
        const std::size_t idx{row * S_MAT + col};
        const fp_t det_diff{
          std::abs(maintained.determinant()[idx] - rebuilt.determinant()[idx])
        };
        max_det_diff = std::max(max_det_diff, det_diff);
      }
    }

    INFO("Maintained and rebuilt determinant matrices should match after each accepted update.");
    CAPTURE(step, moved, max_det_diff);
    REQUIRE(max_det_diff <= (1e-11_fp * RIGOROUS_PRECISION_SCALE));

    const std::size_t probe{(moved + 1U) % n};
    Particles probe_particles{n};
    phy_copy_positions(particles, probe_particles);

    const fp_t probe_dx{0.0017_fp * static_cast<fp_t>((step % 4U) + 1U)};
    const fp_t probe_dy{-0.0011_fp * static_cast<fp_t>((step % 6U) + 1U)};
    const fp_t probe_dz{0.0013_fp * static_cast<fp_t>((step % 5U) + 1U)};

    probe_particles.pos().x_[probe] =
        wrap_coordinate(probe_particles.pos().x_[probe] + probe_dx, box_length);
    probe_particles.pos().y_[probe] =
        wrap_coordinate(probe_particles.pos().y_[probe] + probe_dy, box_length);
    probe_particles.pos().z_[probe] =
        wrap_coordinate(probe_particles.pos().z_[probe] + probe_dz, box_length);

    maintained.update_trig_cache(probe, probe_particles.view());
    rebuilt.update_trig_cache(probe, probe_particles.view());

    const fp_t* const maintained_probe_row{maintained.build_row(probe)};
    const fp_t* const rebuilt_probe_row{rebuilt.build_row(probe)};

    const fp_t maintained_probe_ratio{maintained.determinant_ratio(probe, maintained_probe_row)};
    const fp_t rebuilt_probe_ratio{rebuilt.determinant_ratio(probe, rebuilt_probe_row)};

    maintained.restore_trig_row(probe);
    rebuilt.restore_trig_row(probe);

    INFO("Maintained and rebuilt Slater states should predict the same future determinant ratio.");
    CAPTURE(step, moved, probe, maintained_probe_ratio, rebuilt_probe_ratio);
    REQUIRE(std::isfinite(maintained_probe_ratio));
    REQUIRE(std::isfinite(rebuilt_probe_ratio));
    REQUIRE(std::abs(maintained_probe_ratio - rebuilt_probe_ratio) <= PROBE_RATIO_TOLERANCE);
  }
}

TEST_CASE("WaveFunction log_psi and derivatives are invariant under integer box translations",
          "[rigorous][wavefunction]") {
  constexpr std::size_t n{7U};
  constexpr fp_t box_length{9.5_fp};

  Particles particles{n};
  set_stable_closed_shell_positions(particles);
  for (std::size_t i = 0; i < n; ++i) {
    particles.pos().x_[i] += 0.8_fp * static_cast<fp_t>((3U * i) % 5U);
    particles.pos().y_[i] += 0.6_fp * static_cast<fp_t>((2U * i) % 7U);
    particles.pos().z_[i] += 0.5_fp * static_cast<fp_t>((5U * i) % 3U);
  }

  WaveFunction baseline{particles, box_length};
  const fp_t baseline_log_psi{baseline.evaluate_log_psi(particles.view())};
  baseline.evaluate_derivatives(particles.view());

  std::vector<fp_t> baseline_grad_x(n);
  std::vector<fp_t> baseline_grad_y(n);
  std::vector<fp_t> baseline_grad_z(n);
  std::vector<fp_t> baseline_lap(n);

  for (std::size_t i = 0; i < n; ++i) {
    baseline_grad_x[i] = particles.grad_log_psi().x_[i];
    baseline_grad_y[i] = particles.grad_log_psi().y_[i];
    baseline_grad_z[i] = particles.grad_log_psi().z_[i];
    baseline_lap[i] = particles.lap_log_psi()[i];
  }

  Particles shifted{n};
  phy_copy_positions(particles, shifted);

  shifted.pos().x_[3] += box_length;
  shifted.pos().y_[3] -= 2.0_fp * box_length;
  shifted.pos().z_[3] += 3.0_fp * box_length;

  // After shifting:
  for (std::size_t i = 0; i < n; ++i) {
    shifted.pos().x_[i] -= box_length * std::floor(shifted.pos().x_[i] / box_length);
    shifted.pos().y_[i] -= box_length * std::floor(shifted.pos().y_[i] / box_length);
    shifted.pos().z_[i] -= box_length * std::floor(shifted.pos().z_[i] / box_length);
  }

  WaveFunction translated{shifted, box_length};
  const fp_t translated_log_psi{translated.evaluate_log_psi(shifted.view())};
  translated.evaluate_derivatives(shifted.view());

  INFO("WaveFunction must be periodic under integer box translations.");
  CAPTURE(baseline_log_psi, translated_log_psi);
  REQUIRE(std::abs(baseline_log_psi - translated_log_psi) <= TRANSLATION_LOG_PSI_TOLERANCE);

  for (std::size_t i = 0; i < n; ++i) {
    CAPTURE(i);
    REQUIRE(std::abs(baseline_grad_x[i] - shifted.grad_log_psi().x_[i]) <= TRANSLATION_GRAD_TOLERANCE);
    REQUIRE(std::abs(baseline_grad_y[i] - shifted.grad_log_psi().y_[i]) <= TRANSLATION_GRAD_TOLERANCE);
    REQUIRE(std::abs(baseline_grad_z[i] - shifted.grad_log_psi().z_[i]) <= TRANSLATION_GRAD_TOLERANCE);
    REQUIRE(std::abs(baseline_lap[i] - shifted.lap_log_psi()[i]) <= TRANSLATION_LAP_TOLERANCE);
  }
}

TEST_CASE("EnergyTracker matches an exact Ewald reference on many random small configurations",
          "[rigorous][energy]") {
  constexpr std::size_t n{3U};
  constexpr fp_t box_length{8.25_fp};
  constexpr std::size_t samples{25U};

  std::mt19937_64 rng{987654321ULL};
  std::uniform_real_distribution<fp_t> dist{0.0_fp, box_length};

  for (std::size_t sample = 0; sample < samples; ++sample) {
    Particles particles{n};

    for (std::size_t i = 0; i < n; ++i) {
      particles.pos().x_[i] = dist(rng);
      particles.pos().y_[i] = dist(rng);
      particles.pos().z_[i] = dist(rng);
    }

    EnergyTracker tracker{box_length, n};
    tracker.initialize_structure_factors(particles.view());
    tracker.initialize_reciprocal_energy();
    tracker.initialize_real_energy(particles.view());

    const fp_t tracker_total{tracker.eval_total_energy(particles.view())};
    const fp_t exact_total{exact_total_potential(particles, box_length)};

    INFO("EnergyTracker total potential should match the exact Ewald reference on many random "
         "small configurations.");
    CAPTURE(sample, tracker_total, exact_total);
    REQUIRE(std::abs(tracker_total - exact_total) <= (8e-6_fp * RIGOROUS_PRECISION_SCALE));
  }
}

TEST_CASE("EnergyTracker remains close to the exact Ewald reference across many cached updates",
          "[rigorous][energy]") {
  constexpr std::size_t n{3U};
  constexpr fp_t box_length{8.75_fp};
  constexpr std::size_t steps{48U};

  Particles particles{n};
  particles.pos().x_[0] = 0.55_fp;
  particles.pos().y_[0] = 1.20_fp;
  particles.pos().z_[0] = 2.10_fp;

  particles.pos().x_[1] = 3.15_fp;
  particles.pos().y_[1] = 2.45_fp;
  particles.pos().z_[1] = 1.05_fp;

  particles.pos().x_[2] = 5.60_fp;
  particles.pos().y_[2] = 0.85_fp;
  particles.pos().z_[2] = 4.45_fp;

  EnergyTracker tracker{box_length, n};
  tracker.initialize_structure_factors(particles.view());
  tracker.initialize_reciprocal_energy();
  tracker.initialize_real_energy(particles.view());

  for (std::size_t step = 0; step < steps; ++step) {
    const std::size_t moved{step % n};

    auto pos{particles.pos()};
    const xpu::array<fp_t, idx(Axis::NUM)> old_pos{
      pos[idx(Axis::X)][moved],
      pos[idx(Axis::Y)][moved],
      pos[idx(Axis::Z)][moved]
    };

    const fp_t dx{0.017_fp * static_cast<fp_t>((step % 4U) + 1U)};
    const fp_t dy{-0.012_fp * static_cast<fp_t>((step % 5U) + 1U)};
    const fp_t dz{0.010_fp * static_cast<fp_t>((step % 3U) + 1U)};

    pos[idx(Axis::X)][moved] = wrap_coordinate(
      old_pos[idx(Axis::X)] + dx, box_length
    );
    pos[idx(Axis::Y)][moved] = wrap_coordinate(
      old_pos[idx(Axis::Y)] + dy, box_length
    );
    pos[idx(Axis::Z)][moved] = wrap_coordinate(
      old_pos[idx(Axis::Z)] + dz, box_length
    );

    const xpu::array<fp_t, idx(Axis::NUM)> new_pos{
      pos[idx(Axis::X)][moved],
      pos[idx(Axis::Y)][moved],
      pos[idx(Axis::Z)][moved]
    };

    tracker.update_structure_factors(
      old_pos, new_pos
    );
    tracker.update_real_energy(moved, old_pos, particles.view());

    const fp_t tracker_total{tracker.eval_total_energy(particles.view())};
    const fp_t exact_total{exact_total_potential(particles, box_length)};

    INFO("Cached Ewald updates drifted too far from an exact recomputation.");
    CAPTURE(step, moved, tracker_total, exact_total);
    REQUIRE(std::abs(tracker_total - exact_total) <= (8e-6_fp * RIGOROUS_PRECISION_SCALE));
  }
}
