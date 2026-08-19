/// @file test_validation_known_energies.cpp
/// @brief Validates VMC energies against published Perdew-Zunger total energies
///        for the fully polarized 3D homogeneous electron gas.
///
/// Published reference: Perdew-Zunger thermodynamic-limit total energy per particle
/// for the fully polarized (ζ=1) HEG, constructed from:
///   - Exact HF kinetic + exchange (analytical)
///   - PZ81 correlation fit to Ceperley-Alder DMC data
///   Refs: Ceperley & Alder, PRL 45, 566 (1980)
///         Perdew & Zunger, PRB 23, 5048 (1981)
///
/// Our VMC runs at finite N (7, 19) with a simple Pade-Jastrow (no backflow,
/// unoptimized b=1). We compare E/N against the published thermodynamic-limit
/// values with tolerances that account for:
///   - Finite-size effects (dominant at small N, high density)
///   - Unoptimized Jastrow b parameter
///   - No backflow (VMC ceiling above DMC)
///
/// The key validation: our numbers must be in the right ballpark and converge
/// toward published values as N increases and r_s increases (where finite-size
/// effects shrink).

#include "test_utilities.hpp"

#include "simulation/simulation.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <numbers>

namespace {

/// Published Perdew-Zunger total energy per particle for fully polarized HEG (Ha).
/// E/N = T_HF + E_x + e_c, all at thermodynamic limit.
///
/// Hardcoded from PZ81 parametrization:
///   r_s=1:   +1.1592 Ha
///   r_s=2:   +0.1366 Ha
///   r_s=5:   -0.0537 Ha
///   r_s=10:  -0.0457 Ha
///   r_s=20:  -0.0279 Ha
fp_t published_pz_total_energy(fp_t r_s) {
  // T_HF/N = (2.21 * 2^(2/3)) / (2 * r_s^2) Ha
  constexpr fp_t T_COEFF{2.21_fp * 1.587401051968199_fp * 0.5_fp};
  const fp_t t{T_COEFF / (r_s * r_s)};

  // E_x/N = -(0.9163 * 2^(1/3)) / (2 * r_s) Ha
  constexpr fp_t EX_COEFF{0.9163_fp * 1.2599210498948732_fp * 0.5_fp};
  const fp_t ex{-EX_COEFF / r_s};

  // PZ81 correlation for ζ=1, in Ry -> Ha (divide by 2)
  fp_t ec_ry{};
  if (r_s >= 1.0_fp) {
    constexpr fp_t GAMMA{-0.0843_fp};
    constexpr fp_t BETA1{1.0529_fp};
    constexpr fp_t BETA2{0.3334_fp};
    ec_ry = GAMMA / (1.0_fp + BETA1 * std::sqrt(r_s) + BETA2 * r_s);
  } else {
    constexpr fp_t A{0.01555_fp};
    constexpr fp_t B{-0.0269_fp};
    constexpr fp_t C{0.0007_fp};
    constexpr fp_t D{-0.0048_fp};
    const fp_t ln_rs{std::log(r_s)};
    ec_ry = A * ln_rs + B + C * r_s * ln_rs + D * r_s;
  }

  return t + ex + 0.5_fp * ec_ry;
}

SimResult run_published_vmc(
  std::size_t N,
  fp_t r_s,
  std::size_t warmup_steps,
  std::size_t measure_steps,
  fp_t step_size,
  uint64_t seed,
  std::size_t block_size
) {
  const fp_t L{box_length_from_rs(r_s, N)};
  Config config{make_config(N, L, warmup_steps, measure_steps, step_size, seed, block_size)};

  auto writer{std::make_unique<RecordingOutputWriter>()};
  RecordingOutputWriter* const sink{writer.get()};

  Simulation sim{config, std::move(writer)};
  sim.run();

  SimResult result{};
  if (sink->done.has_value()) {
    result.mean_energy = sink->done->final_mean_energy;
    result.standard_error = sink->done->final_standard_error.value_or(0.0_fp);
    result.acceptance_rate = sink->done->final_acceptance_rate;
  }
  return result;
}

} // namespace

// Validation against Perdew-Zunger published values.
//
// At low r_s (high density), finite-size effects on kinetic energy
// dominate and the VMC value can differ significantly from the
// thermodynamic limit. At higher r_s the agreement tightens.
//
// Tolerance strategy:
//   r_s=5,  N=19: within 0.01  Ha/electron (finite-size ~15%)
//   r_s=10, N=19: within 0.005 Ha/electron (finite-size ~3%)
//   r_s=5,  N=7 vs N=19: convergence check

TEST_CASE("Published values: E/N at r_s=5, N=19 matches PZ81 within 0.01 Ha", "[published]") {
  constexpr std::size_t N{19U};
  constexpr fp_t R_S{5.0_fp};
  const fp_t E_PZ{published_pz_total_energy(R_S)};

  const auto result{run_published_vmc(N, R_S, 5000U, 80000U, 0.6_fp, 5019U, 2000U)};
  const fp_t E_VMC{result.mean_energy / static_cast<fp_t>(N)};

  INFO("E_VMC/N    = " << E_VMC << " Ha");
  INFO("E_PZ/N     = " << E_PZ << " Ha (Perdew-Zunger, thermo limit)");
  INFO("Difference = " << (E_VMC - E_PZ) << " Ha");

  REQUIRE(std::isfinite(E_VMC));
  REQUIRE(std::abs(E_VMC - E_PZ) < 0.01_fp);
}

TEST_CASE("Published values: E/N at r_s=10, N=19 matches PZ81 within 0.005 Ha", "[published]") {
  constexpr std::size_t N{19U};
  constexpr fp_t R_S{10.0_fp};
  const fp_t E_PZ{published_pz_total_energy(R_S)};

  const auto result{run_published_vmc(N, R_S, 5000U, 80000U, 1.0_fp, 10019U, 2000U)};
  const fp_t E_VMC{result.mean_energy / static_cast<fp_t>(N)};

  INFO("E_VMC/N    = " << E_VMC << " Ha");
  INFO("E_PZ/N     = " << E_PZ << " Ha (Perdew-Zunger, thermo limit)");
  INFO("Difference = " << (E_VMC - E_PZ) << " Ha");

  REQUIRE(std::isfinite(E_VMC));
  REQUIRE(std::abs(E_VMC - E_PZ) < 0.005_fp);
}

TEST_CASE("Published values: E/N at r_s=10, N=7 matches PZ81 within 0.005 Ha", "[published]") {
  constexpr std::size_t N{7U};
  constexpr fp_t R_S{10.0_fp};
  const fp_t E_PZ{published_pz_total_energy(R_S)};

  const auto result{run_published_vmc(N, R_S, 5000U, 80000U, 1.2_fp, 10007U, 2000U)};
  const fp_t E_VMC{result.mean_energy / static_cast<fp_t>(N)};

  INFO("E_VMC/N    = " << E_VMC << " Ha");
  INFO("E_PZ/N     = " << E_PZ << " Ha (Perdew-Zunger, thermo limit)");
  INFO("Difference = " << (E_VMC - E_PZ) << " Ha");

  REQUIRE(std::isfinite(E_VMC));
  REQUIRE(std::abs(E_VMC - E_PZ) < 0.005_fp);
}

// At r_s=2 finite-size effects are larger, use wider tolerance
TEST_CASE("Published values: E/N at r_s=2, N=19 matches PZ81 within 0.07 Ha", "[published]") {
  constexpr std::size_t N{19U};
  constexpr fp_t R_S{2.0_fp};
  const fp_t E_PZ{published_pz_total_energy(R_S)};

  const auto result{run_published_vmc(N, R_S, 5000U, 80000U, 0.3_fp, 2019U, 2000U)};
  const fp_t E_VMC{result.mean_energy / static_cast<fp_t>(N)};

  INFO("E_VMC/N    = " << E_VMC << " Ha");
  INFO("E_PZ/N     = " << E_PZ << " Ha (Perdew-Zunger, thermo limit)");
  INFO("Difference = " << (E_VMC - E_PZ) << " Ha");

  REQUIRE(std::isfinite(E_VMC));
  REQUIRE(std::abs(E_VMC - E_PZ) < 0.07_fp);
}

// Convergence toward published values with increasing N
TEST_CASE("Published values: N=19 closer to PZ81 than N=7 at r_s=5", "[published]") {
  constexpr fp_t R_S{5.0_fp};
  const fp_t E_PZ{published_pz_total_energy(R_S)};

  const auto r7{run_published_vmc(7U, R_S, 5000U, 60000U, 0.8_fp, 7500U, 1500U)};
  const auto r19{run_published_vmc(19U, R_S, 5000U, 60000U, 0.6_fp, 19500U, 1500U)};

  const fp_t E7{r7.mean_energy / 7.0_fp};
  const fp_t E19{r19.mean_energy / 19.0_fp};

  const fp_t ERR7{std::abs(E7 - E_PZ)};
  const fp_t ERR19{std::abs(E19 - E_PZ)};

  INFO("E/N (N=7)  = " << E7 << " Ha, |err| = " << ERR7);
  INFO("E/N (N=19) = " << E19 << " Ha, |err| = " << ERR19);
  INFO("E_PZ       = " << E_PZ << " Ha");

  // N=19 should be closer to thermodynamic limit than N=7
  REQUIRE(ERR19 <= ERR7);
}

// Physical correctness checks
TEST_CASE("Published values: correct sign at each density", "[published]") {
  // r_s=1: kinetic dominates, E/N > 0
  {
    const auto r{run_published_vmc(19U, 1.0_fp, 5000U, 40000U, 0.2_fp, 1100U, 1000U)};
    REQUIRE(r.mean_energy / 19.0_fp > 0.0_fp);
  }
  // r_s=5: correlation dominates, E/N < 0
  {
    const auto r{run_published_vmc(19U, 5.0_fp, 5000U, 40000U, 0.6_fp, 5100U, 1000U)};
    REQUIRE(r.mean_energy / 19.0_fp < 0.0_fp);
  }
  // r_s=10: E/N < 0
  {
    const auto r{run_published_vmc(19U, 10.0_fp, 5000U, 40000U, 1.0_fp, 10100U, 1000U)};
    REQUIRE(r.mean_energy / 19.0_fp < 0.0_fp);
  }
}

TEST_CASE("Published values: E/N decreases from r_s=2 to r_s=5", "[published]") {
  const auto r2{run_published_vmc(19U, 2.0_fp, 5000U, 60000U, 0.3_fp, 2200U, 1500U)};
  const auto r5{run_published_vmc(19U, 5.0_fp, 5000U, 60000U, 0.6_fp, 5200U, 1500U)};
  REQUIRE(r5.mean_energy / 19.0_fp < r2.mean_energy / 19.0_fp);
}

TEST_CASE("Published values: deterministic reproducibility", "[published]") {
  const auto r1{run_published_vmc(7U, 5.0_fp, 2000U, 20000U, 0.8_fp, 42U, 500U)};
  const auto r2{run_published_vmc(7U, 5.0_fp, 2000U, 20000U, 0.8_fp, 42U, 500U)};
  REQUIRE(r1.mean_energy == r2.mean_energy);
}

// Print comparison table
TEST_CASE("Published values: print comparison table", "[published][.print]") {
  struct Run {
    fp_t r_s;
    std::size_t N;
    fp_t step;
    uint64_t seed;
  };
  const Run runs[] = {
      {1.0_fp, 7U, 0.2_fp, 1007U},   {1.0_fp, 19U, 0.2_fp, 1019U},   {2.0_fp, 7U, 0.4_fp, 2007U},
      {2.0_fp, 19U, 0.3_fp, 2019U},  {5.0_fp, 7U, 0.8_fp, 5007U},    {5.0_fp, 19U, 0.6_fp, 5019U},
      {10.0_fp, 7U, 1.2_fp, 10007U}, {10.0_fp, 19U, 1.0_fp, 10019U},
  };

  std::cout << "\n=== VMC vs Perdew-Zunger Published Values (Ha/electron) ===\n";
  std::cout << "  r_s    N    E_VMC/N      E_PZ/N      Diff        % err\n";
  std::cout << "  ----  ---  ----------   ----------   ---------   ------\n";

  for (const auto& run : runs) {
    const fp_t e_pz{published_pz_total_energy(run.r_s)};
    const auto result{run_published_vmc(run.N, run.r_s, 5000U, 80000U, run.step, run.seed, 2000U)};
    const fp_t e_vmc{result.mean_energy / static_cast<fp_t>(run.N)};
    const fp_t diff{e_vmc - e_pz};
    const fp_t pct{(std::abs(e_pz) > 1e-10_fp) ? 100.0_fp * diff / e_pz : 0.0_fp};

    std::printf(
      "  %4.1f  %3zu  %+10.6f   %+10.6f   %+9.6f   %+6.1f%%\n",
      run.r_s, run.N, e_vmc, e_pz, diff, pct
    );
  }
  std::cout << std::endl;
}