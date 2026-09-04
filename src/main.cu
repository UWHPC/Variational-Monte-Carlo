#include "config/config.hpp"
#include "optimizer/jastrow_optimizer.hpp"
#include "simulation/simulation.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <print>
#include <string_view>

namespace {

void print_startup() {
  std::print(
    "\n<--- Variational Monte Carlo Simulation --->\n\n"
    "<--- Optimizing Jastrow b parameter --->\n"
  );
}

void print_config(const Config& config) {
  std::print(
    "\n<--- Config Settings --->\n"
    "Number of CPU threads: {}\n"
    "Number of walkers: {}\n"
    "Number of particles: {}\n"
    "Number of warmup sweeps: {}\n"
    "Number of measure sweeps: {}\n"
    "Length of box: {}\n"
    "Samples per block: {}\n"
    "Master seed: {}\n"
    "Jastrow a: {}\n"
    "Jastrow b: {} (optimized)\n\n",
    config.num_threads,
    config.num_walkers,
    config.num_particles,
    config.warmup_sweeps,
    config.measure_sweeps,
    config.box_length,
    config.block_size,
    config.master_seed,
    config.jastrow_a,
    config.jastrow_b
  );
}

void print_summary(
  const Simulation::MeasurementSummary& summary,
  const std::chrono::duration<double>& elapsed
) {
  std::print("<--- Final Measurements --->\nFinal Energy: {:.6}", summary.mean_energy);

  if (summary.standard_error.has_value()) {
    std::print(" +/- {:.6}", *summary.standard_error);
  } else {
    std::print(" +/- N/A (insufficient blocks)");
  }

  constexpr auto percent_scale{100.0_fp};
  std::print(
    "\nElapsed: {} s\nAcceptance Rate: {}%\n\n",
    elapsed.count(),
    summary.acceptance_rate * percent_scale
  );
}

void print_error(const std::string_view message) {
  std::print(stderr, "Exception: {}\n", message);
}

}

int main() {
  try {
    auto config{Config::from_file("config.cfg")};

    print_startup();

    constexpr auto verbose{true};
    const auto optimization{JastrowOptimizer::optimize(config, verbose)};
    config.jastrow_b = optimization.optimal_b;

    print_config(config);

    const auto start{std::chrono::steady_clock::now()};

    Simulation simulation{config};
    const auto summary{simulation.run()};

    const auto end{std::chrono::steady_clock::now()};
    const auto elapsed{std::chrono::duration<double>{end - start}};

    print_summary(summary, elapsed);

    return EXIT_SUCCESS;
  } catch (const std::exception& exception) {
    print_error(exception.what());
    return EXIT_FAILURE;
  }
}
