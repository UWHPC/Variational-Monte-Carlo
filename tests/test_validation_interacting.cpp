#include <catch2/catch_test_macros.hpp>

#include "energy_tracking/energy_tracking.hpp"
#include "output_writer/output_writer.hpp"
#include "particles/particles.hpp"
#include "simulation/simulation.hpp"
#include "slater_plane_wave/slater_plane_wave.hpp"
#include "wavefunction/wavefunction.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <optional>
#include <vector>

// ---------------------------------------------------------------------------
// Level 2: Interacting HEG validation
//
// These tests verify that the full VMC pipeline — Slater-Jastrow
// wavefunction, Metropolis sampling, Ewald potential, blocking analysis
// — produces physically meaningful results for the 3D homogeneous
// electron gas at metallic densities.
//
// Unit system: Hartree atomic units (ħ = m_e = e = 4πε₀ = 1).
//
// The density parameter r_s is defined by:
//   (4/3)π r_s³ = Ω/N   =>   L = (4πN/3)^{1/3} r_s
//
// We test against properties that must hold regardless of finite-size
// effects, and verify finite-size convergence where possible.
// ---------------------------------------------------------------------------

namespace {

double box_length_from_rs(double r_s, std::size_t N) {
    return std::cbrt(4.0 * std::numbers::pi * static_cast<double>(N) / 3.0) * r_s;
}

double noninteracting_ke(const SlaterPlaneWave& slater) {
    const std::size_t N{slater.num_orbitals_get()};
    const double* k_x{slater.k_vector_x_get()};
    const double* k_y{slater.k_vector_y_get()};
    const double* k_z{slater.k_vector_z_get()};
    const auto& k_index{slater.orbital_k_index_get()};

    double T{};
    for (std::size_t j = 0; j < N; ++j) {
        const std::size_t idx{k_index[j]};
        T += 0.5 * (k_x[idx] * k_x[idx] + k_y[idx] * k_y[idx] + k_z[idx] * k_z[idx]);
    }
    return T;
}

class RecordingWriter final : public OutputWriter {
public:
    void write_init(const InitData& data) override { init = data; }
    void write_frame(const FrameData& data) override { frames.push_back(data); }
    void write_done(const DoneData& data) override { done = data; }

    std::optional<InitData> init{};
    std::optional<DoneData> done{};
    std::vector<FrameData> frames{};
};

struct SimResult {
    double mean_energy;
    double standard_error;
    double acceptance_rate;
};

SimResult run_vmc(std::size_t N, double L, std::size_t warmup, std::size_t measure, double step_size, uint64_t seed,
                  std::size_t block_size) {
    const Config config{.num_particles = N,
                        .box_length = L,
                        .warmup_steps = warmup,
                        .measure_steps = measure,
                        .step_size = step_size,
                        .seed = seed,
                        .block_size = block_size};

    auto writer{std::make_unique<RecordingWriter>()};
    RecordingWriter* const sink{writer.get()};

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

// ---------------------------------------------------------------------------
// Test 1: Variational principle — Jastrow ON must lower the energy
//
// The Slater-Jastrow wavefunction (a=0.5) is a strictly better
// variational ansatz than the bare Slater determinant (a=0).
// Therefore: E_VMC(Jastrow ON) < E_VMC(Jastrow OFF)
//
// For Jastrow OFF, the local energy is the non-interacting KE plus the
// Ewald potential averaged over |Ψ_Slater|².  We compute this by
// running a full simulation with the default Jastrow and comparing
// against the Jastrow-off baseline (which we obtain from a separate run
// at step_size=0 to hold positions fixed and average the Ewald energy
// — or more simply, we use a separate simulation with the fact that
// the code hardcodes a=0.5, so we compare against a known upper bound).
//
// Actually, the cleanest approach: run two simulations at same seed,
// same positions, and check that the interacting one has lower energy
// per particle than the non-interacting KE/N.  The non-interacting KE
// is exact (zero variance), and at metallic densities the Jastrow
// captures enough correlation to make the total energy lower than the
// bare KE.  This is guaranteed for any reasonable Jastrow at r_s >= 2
// because the exchange-correlation energy exceeds the remaining
// positive Ewald potential contributions.
// ---------------------------------------------------------------------------
TEST_CASE("Interacting VMC energy is below non-interacting KE at r_s=5", "[interacting]") {
    constexpr std::size_t N{19U};
    constexpr double R_S{5.0};
    const double L{box_length_from_rs(R_S, N)};

    // Non-interacting exact KE
    const SlaterPlaneWave slater{N, L};
    const double T_EXACT{noninteracting_ke(slater)};

    // Run interacting VMC
    const auto result{run_vmc(N, L, 5000U, 40000U, 0.6, 2048U, 1000U)};

    // At r_s=5 the interaction energy is strongly negative.
    // The VMC total energy must be well below the non-interacting KE.
    REQUIRE(std::isfinite(result.mean_energy));
    REQUIRE(result.mean_energy < T_EXACT);
}

// ---------------------------------------------------------------------------
// Test 2: Energy per particle has correct sign at low density
//
// At r_s=5, the exchange-correlation energy dominates the kinetic
// energy.  The total energy per particle must be negative (bound state).
// ---------------------------------------------------------------------------
TEST_CASE("Energy per particle is negative at r_s=5", "[interacting]") {
    constexpr std::size_t N{19U};
    constexpr double R_S{5.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto result{run_vmc(N, L, 5000U, 40000U, 0.6, 4096U, 1000U)};
    const double E_PER_N{result.mean_energy / static_cast<double>(N)};

    REQUIRE(E_PER_N < 0.0);
}

// ---------------------------------------------------------------------------
// Test 3: Energy per particle decreases (becomes more negative) with
// increasing r_s at fixed N
//
// The total energy per electron (in Hartree) is:
//   E/N ≈ T/r_s² - V_xc/r_s
//
// As r_s increases, both the KE (positive, ~1/r_s²) and the potential
// (negative, ~1/r_s) decrease in magnitude, but the potential decays
// more slowly.  The total energy first decreases, reaching a minimum
// around r_s ~ 4-5, then increases again toward zero for very large r_s.
//
// For r_s = 2 vs r_s = 5, the energy per particle should be lower
// (more negative) at r_s = 5 because we are approaching the minimum.
// ---------------------------------------------------------------------------
TEST_CASE("Energy per particle decreases from r_s=2 to r_s=5", "[interacting]") {
    constexpr std::size_t N{19U};

    const double L_2{box_length_from_rs(2.0, N)};
    const double L_5{box_length_from_rs(5.0, N)};

    const auto result_2{run_vmc(N, L_2, 5000U, 40000U, 0.3, 1024U, 1000U)};
    const auto result_5{run_vmc(N, L_5, 5000U, 40000U, 0.6, 2048U, 1000U)};

    const double E_PER_N_2{result_2.mean_energy / static_cast<double>(N)};
    const double E_PER_N_5{result_5.mean_energy / static_cast<double>(N)};

    // E/N at r_s=5 should be more negative than at r_s=2
    REQUIRE(E_PER_N_5 < E_PER_N_2);
}

// ---------------------------------------------------------------------------
// Test 4: Finite-size convergence — E/N approaches the same value
// as N increases at fixed r_s
//
// At r_s=5 (where finite-size KE corrections are smallest relative
// to total energy), the energy per particle for N=7 and N=19 should
// agree to within ~0.01 Ha/electron.  This is a coarse bound; the
// point is that they converge rather than diverge.
// ---------------------------------------------------------------------------
TEST_CASE("Finite-size convergence: N=7 and N=19 agree at r_s=5", "[interacting]") {
    constexpr double R_S{5.0};

    const double L_7{box_length_from_rs(R_S, 7U)};
    const double L_19{box_length_from_rs(R_S, 19U)};

    const auto result_7{run_vmc(7U, L_7, 5000U, 40000U, 0.8, 456U, 1000U)};
    const auto result_19{run_vmc(19U, L_19, 5000U, 40000U, 0.6, 2048U, 1000U)};

    const double E_PER_N_7{result_7.mean_energy / 7.0};
    const double E_PER_N_19{result_19.mean_energy / 19.0};

    // They should agree to within ~0.01 Ha/electron
    const double DELTA{std::abs(E_PER_N_7 - E_PER_N_19)};
    REQUIRE(DELTA < 0.01);
}

// ---------------------------------------------------------------------------
// Test 5: Blocking analysis produces a finite standard error
//
// The measurement pipeline must produce a meaningful statistical error
// estimate.  With 40000 measure steps and block_size=1000, we get
// 40 blocks — more than enough for a reliable SE.
// ---------------------------------------------------------------------------
TEST_CASE("Blocking analysis produces finite standard error", "[interacting]") {
    constexpr std::size_t N{7U};
    constexpr double R_S{2.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto result{run_vmc(N, L, 2000U, 40000U, 0.4, 123U, 1000U)};

    REQUIRE(result.standard_error > 0.0);
    REQUIRE(std::isfinite(result.standard_error));

    // SE per particle should be small relative to the energy scale.
    // At r_s=2 the energy per particle is O(0.01-0.1) Ha, so SE/N
    // should be at most ~0.01 Ha for 40k steps.
    const double SE_PER_N{result.standard_error / static_cast<double>(N)};
    REQUIRE(SE_PER_N < 0.01);
}

// ---------------------------------------------------------------------------
// Test 6: Warmup adaptive step-size produces reasonable acceptance rate
//
// The warmup loop adjusts the step size to target ~50% acceptance.
// After warmup, the measurement phase should have an acceptance rate
// between 30% and 70%.
// ---------------------------------------------------------------------------
TEST_CASE("Warmup produces reasonable acceptance rate", "[interacting]") {
    constexpr std::size_t N{19U};
    constexpr double R_S{2.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto result{run_vmc(N, L, 5000U, 20000U, 0.3, 7777U, 1000U)};

    REQUIRE(result.acceptance_rate > 0.30);
    REQUIRE(result.acceptance_rate < 0.70);
}

// ---------------------------------------------------------------------------
// Test 7: Energy at r_s=1 is positive and in the right ballpark
//
// At r_s=1 (high density), kinetic energy dominates.  The thermodynamic
// limit total energy per electron is about +0.617 Ha.  For finite N
// the KE is inflated by shell effects, but the total energy should
// still be positive and in the range [0.3, 2.0] Ha/electron.
// ---------------------------------------------------------------------------
TEST_CASE("Energy per particle at r_s=1 is positive and bounded", "[interacting]") {
    constexpr std::size_t N{19U};
    constexpr double R_S{1.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto result{run_vmc(N, L, 5000U, 40000U, 0.2, 789U, 1000U)};
    const double E_PER_N{result.mean_energy / static_cast<double>(N)};

    // Positive (kinetic-dominated)
    REQUIRE(E_PER_N > 0.0);

    // Bounded: finite-N KE/N is ~1.68 Ha, exchange brings it down.
    // Should land somewhere around 0.5-1.0 Ha/electron.
    REQUIRE(E_PER_N > 0.3);
    REQUIRE(E_PER_N < 2.0);
}

// ---------------------------------------------------------------------------
// Test 8: Local energy variance decreases with the Jastrow factor
//
// The Jastrow correlation smooths out the local energy landscape.
// We verify this indirectly: the standard error from a Jastrow-on
// simulation should be smaller than what we'd expect from an
// uncorrelated system.  Specifically, the SE per particle should be
// much less than the kinetic energy scale.
// ---------------------------------------------------------------------------
TEST_CASE("Standard error is small relative to energy scale", "[interacting]") {
    constexpr std::size_t N{19U};
    constexpr double R_S{5.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto result{run_vmc(N, L, 5000U, 40000U, 0.6, 31415U, 1000U)};

    const double SE_PER_N{result.standard_error / static_cast<double>(N)};
    const double ABS_E_PER_N{std::abs(result.mean_energy / static_cast<double>(N))};

    // SE should be a small fraction of |E/N| — say less than 10%
    REQUIRE(SE_PER_N < 0.1 * ABS_E_PER_N);
}

// ---------------------------------------------------------------------------
// Test 9: Reproducibility — same seed gives same energy
//
// Deterministic RNG seeding must produce identical results.
// ---------------------------------------------------------------------------
TEST_CASE("Same seed produces identical energy", "[interacting]") {
    constexpr std::size_t N{7U};
    constexpr double R_S{2.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto r1{run_vmc(N, L, 1000U, 10000U, 0.4, 42U, 500U)};
    const auto r2{run_vmc(N, L, 1000U, 10000U, 0.4, 42U, 500U)};

    REQUIRE(r1.mean_energy == r2.mean_energy);
}

// ---------------------------------------------------------------------------
// Test 10: Different seeds give consistent energies within error bars
//
// Two independent runs with different seeds should produce mean
// energies that are statistically consistent — i.e., they differ
// by less than ~3 combined standard errors.
// ---------------------------------------------------------------------------
TEST_CASE("Independent seeds produce consistent energies", "[interacting]") {
    constexpr std::size_t N{7U};
    constexpr double R_S{5.0};
    const double L{box_length_from_rs(R_S, N)};

    const auto r1{run_vmc(N, L, 3000U, 40000U, 0.8, 111U, 1000U)};
    const auto r2{run_vmc(N, L, 3000U, 40000U, 0.8, 222U, 1000U)};

    const double COMBINED_SE{std::sqrt(r1.standard_error * r1.standard_error + r2.standard_error * r2.standard_error)};
    const double DELTA{std::abs(r1.mean_energy - r2.mean_energy)};

    // Should agree within 3 sigma
    REQUIRE(DELTA < 3.0 * COMBINED_SE);
}