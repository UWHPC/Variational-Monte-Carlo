#include "test_utilities.hpp"

#include "optimizer/jastrow_optimizer.hpp"
#include "simulation/simulation.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <numbers>

namespace {

double energy_at_b(std::size_t N, double r_s, double b, uint64_t seed) {
    const double L{box_length_from_rs(r_s, N)};
    Config config{make_config(N, L, 3000U, 40000U, L / 10.0, seed, 1000U)};
    config.jastrow_b = b;

    Simulation sim{config};
    const auto summary{sim.run()};
    return summary.mean_energy;
}

} // namespace

// The optimizer should find a b value, and it should be positive and finite
TEST_CASE("Optimizer: produces valid b parameter", "[optimizer]") {
    constexpr std::size_t N{7U};
    constexpr double R_S{5.0};
    const double L{box_length_from_rs(R_S, N)};

    Config config{};
    config.num_particles = N;
    config.box_length = L;
    config.master_seed = 12345U;
    config.warmup_sweeps = 100U;
    config.measure_sweeps = 100U;
    config.warmup_steps = N * config.warmup_sweeps;
    config.measure_steps = N * config.measure_sweeps;
    config.step_size = L / 10.0;
    config.block_size = 100U;

    const auto result{JastrowOptimizer::optimize(config)};

    INFO("optimal_b = " << result.optimal_b);
    INFO("energy    = " << result.energy);

    REQUIRE(std::isfinite(result.optimal_b));
    REQUIRE(result.optimal_b > 0.0);
    REQUIRE(std::isfinite(result.energy));
}

// The optimized b should give lower (or equal) energy than b=1.0
TEST_CASE("Optimizer: optimized b lowers energy vs default b=1", "[optimizer]") {
    constexpr std::size_t N{7U};
    constexpr double R_S{5.0};
    const double L{box_length_from_rs(R_S, N)};

    Config config{};
    config.num_particles = N;
    config.box_length = L;
    config.master_seed = 54321U;
    config.warmup_sweeps = 100U;
    config.measure_sweeps = 100U;
    config.warmup_steps = N * config.warmup_sweeps;
    config.measure_steps = N * config.measure_sweeps;
    config.step_size = L / 10.0;
    config.block_size = 100U;

    const auto result{JastrowOptimizer::optimize(config)};

    // Run a longer simulation at the optimized b and at the default b=1
    const double E_OPTIMIZED{energy_at_b(N, R_S, result.optimal_b, 9999U)};
    const double E_DEFAULT{energy_at_b(N, R_S, 1.0, 9999U)};

    INFO("E(b_opt=" << result.optimal_b << ") = " << E_OPTIMIZED);
    INFO("E(b=1.0) = " << E_DEFAULT);

    REQUIRE(E_OPTIMIZED <= E_DEFAULT);
}

// At r_s=10 the optimizer should also work
TEST_CASE("Optimizer: works at r_s=10", "[optimizer]") {
    constexpr std::size_t N{7U};
    constexpr double R_S{10.0};
    const double L{box_length_from_rs(R_S, N)};

    Config config{};
    config.num_particles = N;
    config.box_length = L;
    config.master_seed = 77777U;
    config.warmup_sweeps = 100U;
    config.measure_sweeps = 100U;
    config.warmup_steps = N * config.warmup_sweeps;
    config.measure_steps = N * config.measure_sweeps;
    config.step_size = L / 10.0;
    config.block_size = 100U;

    const auto result{JastrowOptimizer::optimize(config)};

    INFO("r_s=10, optimal_b = " << result.optimal_b);
    INFO("energy = " << result.energy);

    REQUIRE(std::isfinite(result.optimal_b));
    REQUIRE(result.optimal_b > 0.0);
    REQUIRE(result.optimal_b < 20.0);
}
