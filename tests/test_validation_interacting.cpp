#include "test_utilities.hpp"

#include "simulation/simulation.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>

namespace {

SimResult run_vmc(
  std::size_t N,
  real_t L,
  std::size_t warmup,
  std::size_t measure,
  real_t step_size,
  uint64_t seed,
  std::size_t block_size
) {
  const Config config{make_config(N, L, warmup, measure, step_size, seed, block_size)};

  auto writer{std::make_unique<RecordingOutputWriter>()};
  RecordingOutputWriter* const sink{writer.get()};

  Simulation sim{config, std::move(writer)};
  sim.run();

  SimResult result{};
  if (sink->done.has_value()) {
    result.mean_energy = sink->done->final_mean_energy;
    result.standard_error = sink->done->final_standard_error.value_or(0.0_r);
    result.acceptance_rate = sink->done->final_acceptance_rate;
  }
  return result;
}

} // namespace

// At r_s=5 the interaction energy is strongly negative, VMC total must be below non-interacting KE
TEST_CASE("Interacting VMC energy is below non-interacting KE at r_s=5", "[interacting]") {
  constexpr std::size_t N{19U};
  constexpr real_t R_S{5.0_r};
  const real_t L{box_length_from_rs(R_S, N)};

  Particles p{N};
  const SlaterPlaneWave slater{p, L};
  const real_t T_EXACT{exact_kinetic_energy(slater)};

  const auto result{run_vmc(N, L, 5000U, 40000U, 0.6_r, 2048U, 1000U)};

  REQUIRE(std::isfinite(result.mean_energy));
  REQUIRE(result.mean_energy < T_EXACT);
}

// At r_s=5, exchange-correlation dominates -> E/N must be negative
TEST_CASE("Energy per particle is negative at r_s=5", "[interacting]") {
  constexpr std::size_t N{19U};
  constexpr real_t R_S{5.0_r};
  const real_t L{box_length_from_rs(R_S, N)};

  const auto result{run_vmc(N, L, 5000U, 40000U, 0.6_r, 4096U, 1000U)};
  const real_t E_PER_N{result.mean_energy / static_cast<real_t>(N)};

  REQUIRE(E_PER_N < 0.0_r);
}

// E/N at r_s=5 should be more negative than at r_s=2 (approaching energy minimum)
TEST_CASE("Energy per particle decreases from r_s=2 to r_s=5", "[interacting]") {
  constexpr std::size_t N{19U};

  const real_t L_2{box_length_from_rs(2.0_r, N)};
  const real_t L_5{box_length_from_rs(5.0_r, N)};

  const auto result_2{run_vmc(N, L_2, 5000U, 40000U, 0.3_r, 1024U, 1000U)};
  const auto result_5{run_vmc(N, L_5, 5000U, 40000U, 0.6_r, 2048U, 1000U)};

  const real_t E_PER_N_2{result_2.mean_energy / static_cast<real_t>(N)};
  const real_t E_PER_N_5{result_5.mean_energy / static_cast<real_t>(N)};

  REQUIRE(E_PER_N_5 < E_PER_N_2);
}

// E/N for N=7 and N=19 should agree within ~0.01 Ha/electron at r_s=5
TEST_CASE("Finite-size convergence: N=7 and N=19 agree at r_s=5", "[interacting]") {
  constexpr real_t R_S{5.0_r};

  const real_t L_7{box_length_from_rs(R_S, 7U)};
  const real_t L_19{box_length_from_rs(R_S, 19U)};

  const auto result_7{run_vmc(7U, L_7, 5000U, 40000U, 0.8_r, 456U, 1000U)};
  const auto result_19{run_vmc(19U, L_19, 5000U, 40000U, 0.6_r, 2048U, 1000U)};

  const real_t E_PER_N_7{result_7.mean_energy / 7.0_r};
  const real_t E_PER_N_19{result_19.mean_energy / 19.0_r};

  const real_t DELTA{std::abs(E_PER_N_7 - E_PER_N_19)};
  REQUIRE(DELTA < 0.01_r);
}

// 40k measure steps with block_size=1000 -> 40 blocks, should produce reliable SE
TEST_CASE("Blocking analysis produces finite standard error", "[interacting]") {
  constexpr std::size_t N{7U};
  constexpr real_t R_S{2.0_r};
  const real_t L{box_length_from_rs(R_S, N)};

  const auto result{run_vmc(N, L, 2000U, 40000U, 0.4_r, 123U, 1000U)};

  REQUIRE(result.standard_error > 0.0_r);
  REQUIRE(std::isfinite(result.standard_error));

  const real_t SE_PER_N{result.standard_error / static_cast<real_t>(N)};
  REQUIRE(SE_PER_N < 0.01_r);
}

// Warmup targets ~50% acceptance, measurement should land between 30-70%
TEST_CASE("Warmup produces reasonable acceptance rate", "[interacting]") {
  constexpr std::size_t N{19U};
  constexpr real_t R_S{2.0_r};
  const real_t L{box_length_from_rs(R_S, N)};

  const auto result{run_vmc(N, L, 5000U, 20000U, 0.3_r, 7777U, 1000U)};

  REQUIRE(result.acceptance_rate > 0.30_r);
  REQUIRE(result.acceptance_rate < 0.70_r);
}

// At r_s=1 (high density), KE dominates -> E/N positive, bounded [0.3, 2.0] Ha
TEST_CASE("Energy per particle at r_s=1 is positive and bounded", "[interacting]") {
  constexpr std::size_t N{19U};
  constexpr real_t R_S{1.0_r};
  const real_t L{box_length_from_rs(R_S, N)};

  const auto result{run_vmc(N, L, 5000U, 40000U, 0.2_r, 789U, 1000U)};
  const real_t E_PER_N{result.mean_energy / static_cast<real_t>(N)};

  REQUIRE(E_PER_N > 0.0_r);
  REQUIRE(E_PER_N > 0.3_r);
  REQUIRE(E_PER_N < 2.0_r);
}

// SE per particle should be small fraction of |E/N|
TEST_CASE("Standard error is small relative to energy scale", "[interacting]") {
  constexpr std::size_t N{19U};
  constexpr real_t R_S{5.0_r};
  const real_t L{box_length_from_rs(R_S, N)};

  const auto result{run_vmc(N, L, 5000U, 40000U, 0.6_r, 31415U, 1000U)};

  const real_t SE_PER_N{result.standard_error / static_cast<real_t>(N)};
  const real_t ABS_E_PER_N{std::abs(result.mean_energy / static_cast<real_t>(N))};

  REQUIRE(SE_PER_N < 0.1_r * ABS_E_PER_N);
}

// Deterministic RNG seeding must produce identical results
TEST_CASE("Same seed produces identical energy", "[interacting]") {
  constexpr std::size_t N{7U};
  constexpr real_t R_S{2.0_r};
  const real_t L{box_length_from_rs(R_S, N)};

  const auto r1{run_vmc(N, L, 1000U, 10000U, 0.4_r, 42U, 500U)};
  const auto r2{run_vmc(N, L, 1000U, 10000U, 0.4_r, 42U, 500U)};

  REQUIRE(r1.mean_energy == r2.mean_energy);
}

// Different seeds should agree within 3 combined standard errors
TEST_CASE("Independent seeds produce consistent energies", "[interacting]") {
  constexpr std::size_t N{7U};
  constexpr real_t R_S{5.0_r};
  const real_t L{box_length_from_rs(R_S, N)};

  const auto r1{run_vmc(N, L, 3000U, 40000U, 0.8_r, 111U, 1000U)};
  const auto r2{run_vmc(N, L, 3000U, 40000U, 0.8_r, 222U, 1000U)};

  const real_t COMBINED_SE{
    std::sqrt(r1.standard_error * r1.standard_error + r2.standard_error * r2.standard_error)
  };
  const real_t DELTA{std::abs(r1.mean_energy - r2.mean_energy)};

  REQUIRE(DELTA < 3.0_r * COMBINED_SE);
}