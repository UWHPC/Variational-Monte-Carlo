#include "jastrow_optimizer.hpp"
#include "../simulation/simulation.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <iostream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <vector>

double JastrowOptimizer::compute_rs(const Config& cfg) {
    const double N{static_cast<double>(cfg.num_particles)};
    const double volume{cfg.box_length * cfg.box_length * cfg.box_length};
    const double density{N / volume};
    return std::cbrt(3.0 / (4.0 * std::numbers::pi * density));
}


JastrowOptimizer::Result JastrowOptimizer::optimize(const Config& base_config, bool verbose) {
    const double R_S{compute_rs(base_config)};
    const std::size_t N{base_config.num_particles};

    // b controls the inverse Jastrow range: effective range ~ 1/b.
    // b_min = 1/r_s: Jastrow range ≤ one mean interparticle spacing.
    //         Beyond this, correlations are long-range and should be handled
    //         by the Ewald sum, not the Jastrow factor.
    // b_max = 5.0: at b=5 the Jastrow range is 0.2 bohr, capturing only
    //         the electron-electron cusp with negligible correlation beyond.
    //         The energy is flat well before this point.
    const double b_min{1.0 / R_S};
    const double b_max{std::max(5.0, b_min + 0.5)};
    const std::size_t grid_points{2U * base_config.num_threads};

    // Adapt scan statistics based on system size.
    const std::size_t scan_warmup{std::max<std::size_t>(50U, 500U / std::max<std::size_t>(N / 7U, 1U))};
    const std::size_t scan_measure{std::max<std::size_t>(50U, 500U / std::max<std::size_t>(N / 7U, 1U))};

    if (verbose) {
        std::cout << "[Optimizer] r_s=" << std::fixed << std::setprecision(2) << R_S
                  << ", N=" << N
                  << ", scanning b in [" << b_min << ", " << b_max << "]"
                  << " (" << grid_points << " points, "
                  << scan_warmup << "+" << scan_measure << " sweeps)\n";
    }

    // Build b values
    std::vector<double> b_values(grid_points);
    for (std::size_t i = 0; i < grid_points; ++i) {
        b_values[i] = b_min + (b_max - b_min) * static_cast<double>(i)
                      / static_cast<double>(grid_points - 1);
    }

    // Phase 1: parallel grid scan
    std::vector<std::future<EvalResult>> futures;
    futures.reserve(grid_points);

    for (std::size_t i = 0; i < grid_points; ++i) {
        const double b{b_values[i]};
        const std::size_t warmup{scan_warmup};
        const std::size_t measure{scan_measure};
        futures.push_back(std::async(std::launch::async,
            [&base_config, b, warmup, measure]() {
                return evaluate(base_config, b, warmup, measure);
            }));
    }

    std::vector<EvalResult> results;
    results.reserve(grid_points);

    for (auto& f : futures) {
        results.push_back(f.get());
        if (verbose) {
            const auto& r{results.back()};
            std::cout << "  b=" << std::fixed << std::setprecision(3) << r.b
                      << "  E=" << std::setprecision(6) << r.energy << "\n";
        }
    }

    // Find the best point, but guard against boundary artifacts.
    // If the lowest energy is at b_min and it's significantly lower than the
    // second-best, the long-range Jastrow bias is dominating. In that case,
    // skip the boundary point and take the best of the interior points.
    // "Significantly lower" = more than 2x the spread of the interior points.
    double best_b{1.0};
    double best_energy{std::numeric_limits<double>::max()};
    double second_best_energy{std::numeric_limits<double>::max()};

    for (const auto& r : results) {
        if (r.energy < best_energy) {
            second_best_energy = best_energy;
            best_energy = r.energy;
            best_b = r.b;
        } else if (r.energy < second_best_energy) {
            second_best_energy = r.energy;
        }
    }

    // Check if the winner is the boundary point and an outlier
    const bool at_boundary{std::abs(best_b - b_min) < 1e-10};
    if (at_boundary && results.size() > 2) {
        // Compute the energy spread of all non-boundary points
        double interior_min{std::numeric_limits<double>::max()};
        double interior_max{std::numeric_limits<double>::lowest()};
        for (const auto& r : results) {
            if (std::abs(r.b - b_min) < 1e-10) continue;
            interior_min = std::min(interior_min, r.energy);
            interior_max = std::max(interior_max, r.energy);
        }
        const double interior_spread{interior_max - interior_min};
        const double boundary_gap{interior_min - best_energy};

        // If the boundary point is more than 2x the interior spread below
        // the interior minimum, it's an artifact; use the interior best instead.
        if (boundary_gap > 2.0 * interior_spread) {
            if (verbose) {
                std::cout << "[Optimizer] Boundary b=" << std::setprecision(3) << best_b
                          << " rejected (artifact), using interior best\n";
            }
            best_b = 1.0;
            best_energy = std::numeric_limits<double>::max();
            for (const auto& r : results) {
                if (std::abs(r.b - b_min) < 1e-10) continue;
                if (r.energy < best_energy) {
                    best_energy = r.energy;
                    best_b = r.b;
                }
            }
        }
    }

    if (verbose) {
        std::cout << "[Optimizer] Result: b=" << std::setprecision(4) << best_b << "\n";
    }

    return Result{
        .optimal_b = best_b,
        .energy = best_energy,
        .standard_error = 0.0
    };
}

JastrowOptimizer::EvalResult JastrowOptimizer::evaluate(const Config& base_config, double b,
                                                         std::size_t warmup_sweeps,
                                                         std::size_t measure_sweeps) {
    Config cfg{};
    cfg.num_particles = base_config.num_particles;
    cfg.box_length = base_config.box_length;
    cfg.jastrow_a = base_config.jastrow_a;
    cfg.jastrow_b = b;
    cfg.master_seed = base_config.master_seed;
    cfg.is_master_thread = false;
    cfg.num_threads = 1;

    cfg.warmup_sweeps = warmup_sweeps;
    cfg.measure_sweeps = measure_sweeps;
    cfg.block_size = std::max<std::size_t>(50U, measure_sweeps / 5U);
    cfg.warmup_steps = cfg.num_particles * warmup_sweeps;
    cfg.measure_steps = cfg.num_particles * measure_sweeps;
    cfg.step_size = base_config.box_length / 10.0;

    Simulation sim{cfg};
    const auto summary{sim.run()};

    return EvalResult{
        .b = b,
        .energy = summary.mean_energy,
        .standard_error = summary.standard_error.value_or(0.0)
    };
}