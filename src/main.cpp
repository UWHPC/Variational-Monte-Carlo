#include "config/config.hpp"
#include "output_writer/output_writer.hpp"
#include "simulation/simulation.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <omp.h>

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {
    // Temporary local config while CLI parsing is disabled.

    // // Closed shells (up to 10k):
    // N = 1, 7, 19, 27, 33, 57, 81, 93, 123, 147, 171, 179,
    // 203, 251, 257, 305, 341, 365, 389, 437, 461, 485, 515, 587,
    // 619, 691, 739, 751, 799, 847, 895, 925, 949, 1021, 1045, 1141,
    // 1189, 1213, 1237, 1309, 1357, 1365, 1419, 1503, 1551, 1575, 1647, 1743,
    // 1791, 1839, 1863, 1935, 2007, 2103, 2109, 2205, 2301, 2325, 2373, 2469,
    // 2517, 2553, 2601, 2721, 2777, 2801, 2897, 2945, 2969, 3071, 3119, 3191,
    // 3239, 3287, 3407, 3431, 3575, 3695, 3743, 3791, 3887, 3911, 3959, 4067,
    // 4139, 4169, 4337, 4385, 4457, 4553, 4625, 4697, 4729, 4801, 4945, 5041,
    // 5137, 5185, 5257, 5377, 5449, 5497, 5575, 5695, 5743, 5887, 6031, 6043,
    // 6187, 6235, 6355, 6403, 6451, 6619, 6667, 6763, 6859, 6931, 6979, 7075,
    // 7123, 7153, 7249, 7441, 7497, 7521, 7689, 7809, 7881, 8025, 8121, 8217,
    // 8289, 8385, 8409, 8601, 8709, 8733, 8829, 8925, 9045, 9093, 9171, 9315,
    // 9435, 9459, 9627, 9771, 9795, 9843, 9939, 10059
    const Config config{.num_particles = 1791U,
                        .box_length = 9.0,
                        .warmup_steps = 500U,
                        .measure_steps = 2000U,
                        .step_size = 0.05,
                        .seed = 12345U,
                        .block_size = 500U};

    // const Config config{parse_args(argc, argv)};
    // std::ofstream out_file{"data/run.jsonl"};
    // if (!out_file) {
    //     std::cerr << "Failed to open output file: data/run.jsonl\n";
    //     return 1;
    // }
    // std::unique_ptr<OutputWriter> writer{make_output_writer(OutputFormat::JSON, out_file)};

    Simulation sim{config};

    auto start{std::chrono::steady_clock::now()};
    sim.run();
    auto end{std::chrono::steady_clock::now()};

    std::chrono::duration<double> elapsed{end - start};
    std::cout << "Elapsed: " << elapsed.count() << " s" << std::endl;

    return 0;
}
