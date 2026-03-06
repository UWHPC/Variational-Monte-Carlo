#include "config/config.hpp"
#include "output_writer/output_writer.hpp"
#include "simulation/simulation.hpp"

#include <fstream>
#include <iostream>
#include <memory>

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {
    // Temporary local config while CLI parsing is disabled.
    const Config config{.num_particles = 50U,
                        .box_length = 6.75,
                        .warmup_steps = 5000U,
                        .measure_steps = 50000U,
                        .step_size = 0.20,
                        .seed = 12345U,
                        .block_size = 200U};

    // const Config config{parse_args(argc, argv)};
    std::ofstream out_file{"data/run.jsonl"};
    if (!out_file) {
        std::cerr << "Failed to open output file: data/run.jsonl\n";
        return 1;
    }

    std::unique_ptr<OutputWriter> writer{make_output_writer(OutputFormat::JSON, out_file)};
    Simulation sim{config, std::move(writer)};
    sim.run();
    return 0;
}
