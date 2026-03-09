#include "config/config.hpp"
#include "output_writer/output_writer.hpp"
#include "simulation/simulation.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {
    // Temporary local config while CLI parsing is disabled.
    // Closed shells: N = 1, 7, 19, 27, 33, 57, 81, 93, 123, 147, 171, 179, 203, 251, 257, 305, 341, 365, 389, 437, 461,
    // 485
    const Config config{.num_particles = 203U,
                        .box_length = 6.75,
                        .warmup_steps = 2000U,
                        .measure_steps = 20000U,
                        .step_size = 0.05,
                        .seed = 12345U,
                        .block_size = 1000U};

    // const Config config{parse_args(argc, argv)};
    std::ofstream out_file{"data/run.jsonl"};
    if (!out_file) {
        std::cerr << "Failed to open output file: data/run.jsonl\n";
        return 1;
    }
    std::unique_ptr<OutputWriter> writer{make_output_writer(OutputFormat::JSON, out_file)};

    Simulation sim{config, std::move(writer)};

    auto start{std::chrono::steady_clock::now()};
    sim.run();
    auto end{std::chrono::steady_clock::now()};

    std::chrono::duration<double> elapsed{end - start};
    std::cout << "Elapsed: " << elapsed.count() << " s" << std::endl;

    return 0;
}
