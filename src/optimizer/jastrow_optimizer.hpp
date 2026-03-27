#pragma once

#include "../config/config.hpp"

class JastrowOptimizer {
public:
    struct Result {
        double optimal_b;
        double energy;
        double standard_error;
    };

    // Optimize the Jastrow b parameter for the given config.
    // Phase 1: parallel grid scan to locate the basin.
    // Phase 2: serial golden-section refinement to pin down the minimum.
    // Uses variance-penalized energy to avoid unphysical long-range Jastrow states.
    [[nodiscard]] static Result optimize(const Config& base_config, bool verbose = false);

private:
    struct EvalResult {
        double b;
        double energy;
        double standard_error;
    };

    [[nodiscard]] static EvalResult evaluate(const Config& base_config, double b,
                                              std::size_t warmup_sweeps, std::size_t measure_sweeps);
    [[nodiscard]] static double compute_rs(const Config& cfg);
};