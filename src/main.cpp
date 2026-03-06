#include "config/config.hpp"
#include "output_writer/output_writer.hpp"
#include "simulation/simulation.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {
    // auto writer = make_output_writer(OutputFormat::JSON, out_file);

    // Temporary local config while CLI parsing is disabled.
    const Config config{.num_particles = 16U,
                        .box_length = 4.5,
                        .warmup_steps = 50U,
                        .measure_steps = 200U,
                        .step_size = 0.25,
                        .seed = 12345U,
                        .block_size = 20U};

    // const Config config{parse_args(argc, argv)};
    Simulation sim{config};
    sim.run();
    return 0;
}
