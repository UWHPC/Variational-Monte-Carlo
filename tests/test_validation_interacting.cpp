#include "test_utilities.hpp"

#include "simulation/simulation.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>

namespace {

SimResult run_vmc(std::size_t N, double L, std::size_t warmup, std::size_t measure, double step_size, uint64_t seed,
                  std::size_t block_size) {
    const Config config{make_config(N, L, warmup, measure, step_size, seed, block_size)};

    auto writer{std::make_unique<RecordingOutputWriter>()};
    RecordingOutputWriter* const sink{writer.get()};

    Simulation sim{config, std::move(writer)};
    sim.run();

    SimResult result{};
    if (sink->done.has_value()) {
        result.mean_energy = sink->done->final_mean_energy;
        result.standard_error = sink->done->final_standard_error.value_or(0.0);
        result.acceptance_rate = sink->done->final_acceptance_rate;
    }
    return result;
}

} // namespace

// At r_s=5 the interaction energy is strongly negative, VMC total must be below non-interacting KE
TEST_CASE("Interacting VMC energy is below non-interacting KE at r_s=5", "[interacting]") {
    constexpr std::size_t N{19U};
    constexpr double R_S{5.0};
    const double L{box_length_from_rs(R_S, N)};

    Particles p{N};
    const SlaterPlaneWave slater{p, L};
    const double T_EXACT{exact_kinetic_energy(slater)};

    const auto result{run_vmc(N, L, 5000U, 40000U, 0.6, 2048U, 1000U)};

    REQUIRE(std::isfinite(result.mean_energy));
    REQUIRE(result.mean_energy < T_EXACT);
}

// At r_s=5, exchange-correlation dominates -> E/N must be negative
TEST_CASE("Energy per particle is negative at r_s=5", "[interacting]") {
    constexpr std::size_t N{19U};
    constexpr double R_S{5.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto result{run_vmc(N, L, 5000U, 40000U, 0.6, 4096U, 1000U)};
    const double E_PER_N{result.mean_energy / static_cast<double>(N)};

    REQUIRE(E_PER_N < 0.0);
}

// E/N at r_s=5 should be more negative than at r_s=2 (approaching energy minimum)
TEST_CASE("Energy per particle decreases from r_s=2 to r_s=5", "[interacting]") {
    constexpr std::size_t N{19U};

    const double L_2{box_length_from_rs(2.0, N)};
    const double L_5{box_length_from_rs(5.0, N)};

    const auto result_2{run_vmc(N, L_2, 5000U, 40000U, 0.3, 1024U, 1000U)};
    const auto result_5{run_vmc(N, L_5, 5000U, 40000U, 0.6, 2048U, 1000U)};

    const double E_PER_N_2{result_2.mean_energy / static_cast<double>(N)};
    const double E_PER_N_5{result_5.mean_energy / static_cast<double>(N)};

    REQUIRE(E_PER_N_5 < E_PER_N_2);
}

// E/N for N=7 and N=19 should agree within ~0.01 Ha/electron at r_s=5
TEST_CASE("Finite-size convergence: N=7 and N=19 agree at r_s=5", "[interacting]") {
    constexpr double R_S{5.0};

    const double L_7{box_length_from_rs(R_S, 7U)};
    const double L_19{box_length_from_rs(R_S, 19U)};

    const auto result_7{run_vmc(7U, L_7, 5000U, 40000U, 0.8, 456U, 1000U)};
    const auto result_19{run_vmc(19U, L_19, 5000U, 40000U, 0.6, 2048U, 1000U)};

    const double E_PER_N_7{result_7.mean_energy / 7.0};
    const double E_PER_N_19{result_19.mean_energy / 19.0};

    const double DELTA{std::abs(E_PER_N_7 - E_PER_N_19)};
    REQUIRE(DELTA < 0.01);
}

// 40k measure steps with block_size=1000 -> 40 blocks, should produce reliable SE
TEST_CASE("Blocking analysis produces finite standard error", "[interacting]") {
    constexpr std::size_t N{7U};
    constexpr double R_S{2.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto result{run_vmc(N, L, 2000U, 40000U, 0.4, 123U, 1000U)};

    REQUIRE(result.standard_error > 0.0);
    REQUIRE(std::isfinite(result.standard_error));

    const double SE_PER_N{result.standard_error / static_cast<double>(N)};
    REQUIRE(SE_PER_N < 0.01);
}

// Warmup targets ~50% acceptance, measurement should land between 30-70%
TEST_CASE("Warmup produces reasonable acceptance rate", "[interacting]") {
    constexpr std::size_t N{19U};
    constexpr double R_S{2.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto result{run_vmc(N, L, 5000U, 20000U, 0.3, 7777U, 1000U)};

    REQUIRE(result.acceptance_rate > 0.30);
    REQUIRE(result.acceptance_rate < 0.70);
}

// At r_s=1 (high density), KE dominates -> E/N positive, bounded [0.3, 2.0] Ha
TEST_CASE("Energy per particle at r_s=1 is positive and bounded", "[interacting]") {
    constexpr std::size_t N{19U};
    constexpr double R_S{1.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto result{run_vmc(N, L, 5000U, 40000U, 0.2, 789U, 1000U)};
    const double E_PER_N{result.mean_energy / static_cast<double>(N)};

    REQUIRE(E_PER_N > 0.0);
    REQUIRE(E_PER_N > 0.3);
    REQUIRE(E_PER_N < 2.0);
}

// SE per particle should be small fraction of |E/N|
TEST_CASE("Standard error is small relative to energy scale", "[interacting]") {
    constexpr std::size_t N{19U};
    constexpr double R_S{5.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto result{run_vmc(N, L, 5000U, 40000U, 0.6, 31415U, 1000U)};

    const double SE_PER_N{result.standard_error / static_cast<double>(N)};
    const double ABS_E_PER_N{std::abs(result.mean_energy / static_cast<double>(N))};

    REQUIRE(SE_PER_N < 0.1 * ABS_E_PER_N);
}

// Deterministic RNG seeding must produce identical results
TEST_CASE("Same seed produces identical energy", "[interacting]") {
    constexpr std::size_t N{7U};
    constexpr double R_S{2.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto r1{run_vmc(N, L, 1000U, 10000U, 0.4, 42U, 500U)};
    const auto r2{run_vmc(N, L, 1000U, 10000U, 0.4, 42U, 500U)};

    REQUIRE(r1.mean_energy == r2.mean_energy);
}

// Different seeds should agree within 3 combined standard errors
TEST_CASE("Independent seeds produce consistent energies", "[interacting]") {
    constexpr std::size_t N{7U};
    constexpr double R_S{5.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto r1{run_vmc(N, L, 3000U, 40000U, 0.8, 111U, 1000U)};
    const auto r2{run_vmc(N, L, 3000U, 40000U, 0.8, 222U, 1000U)};

    const double COMBINED_SE{std::sqrt(r1.standard_error * r1.standard_error + r2.standard_error * r2.standard_error)};
    const double DELTA{std::abs(r1.mean_energy - r2.mean_energy)};

    REQUIRE(DELTA < 3.0 * COMBINED_SE);
}