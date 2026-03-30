#include "config/config.hpp"
#include "optimizer/jastrow_optimizer.hpp"
#include "output_writer/output_writer.hpp"
#include "simulation/simulation.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <future>
#include <thread>
#include <iomanip>
#include <memory>

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {
    Config master_config{Config::from_file("config.cfg")};

    try {
        std::size_t num_threads{master_config.num_threads};

        // Optimize Jastrow b parameter before production run
        std::cout << "\n<--- Variational Monte Carlo Simulation --->\n\n"
                  << "<--- Optimizing Jastrow b parameter --->\n";

        const auto opt_result{JastrowOptimizer::optimize(master_config, true)};
        master_config.jastrow_b = opt_result.optimal_b;

        std::cout << "\n<--- Config Settings --->\n"
                  << "Number of threads: " << num_threads << "\n"
                  << "Number of particles: " << master_config.num_particles << "\n"
                  << "Number of warmup sweeps: " << master_config.warmup_sweeps << "\n"
                  << "Number of measure sweeps: " << master_config.measure_sweeps << "\n"
                  << "Length of box: " << master_config.box_length << "\n"
                  << "Samples per block: " << master_config.block_size << "\n"
                  << "Master seed: " << master_config.master_seed << "\n"
                  << "Jastrow a: " << master_config.jastrow_a << "\n"
                  << "Jastrow b: " << master_config.jastrow_b << " (optimized)\n"
                  << std::endl;

        std::mt19937_64 master_rng(master_config.master_seed);
        std::uniform_int_distribution<uint64_t> seedDist;
        std::vector<std::future<Simulation::MeasurementSummary>> futures;

        std::ofstream bin_out{"output/vmc.bin", std::ios::binary | std::ios::trunc};
        if (!bin_out) {
            std::cerr << "Warning: could not open output/vmc.bin for writing\n";
        }

        auto start{std::chrono::steady_clock::now()};

        for (std::size_t thread{}; thread < num_threads; ++thread) {
            Config thread_config = master_config;
            thread_config.master_seed = seedDist(master_rng);
            thread_config.is_master_thread = (thread == 0);

            if (thread == 0 && bin_out) {
                futures.push_back(std::async(std::launch::async,
                    [thread_config, &bin_out]() {
                        auto writer{std::make_unique<BinOutputWriter>(bin_out)};
                        Simulation sim{thread_config, std::move(writer)};
                        return sim.run();
                    }));
            } else {
                futures.push_back(std::async(std::launch::async, [thread_config]() {
                    Simulation sim{thread_config};
                    return sim.run();
                }));
            }
        }
        
        double global_energy_sum{};
        double global_variance_sum{};
        double global_acceptance_sum{};
        std::size_t threads_with_se{};

        for (auto& f : futures) {
            Simulation::MeasurementSummary summary{f.get()};
            global_energy_sum += summary.mean_energy;
            if (summary.standard_error.has_value()) {
                global_variance_sum += (*summary.standard_error) * (*summary.standard_error);
                ++threads_with_se;
            }
            global_acceptance_sum += summary.acceptance_rate;
        }

        auto end{std::chrono::steady_clock::now()};
        std::chrono::duration<double> elapsed{end - start};

        const double inv_num_threads{1.0 / static_cast<double>(num_threads)};
        const double final_mean{global_energy_sum * inv_num_threads};
        const double final_acceptance_rate{global_acceptance_sum * inv_num_threads * 100.0};

        std::cout << "<--- Final Measurements --->" << std::endl
                  << "Final Aggregated Energy: " << std::setprecision(6) << final_mean;

        if (threads_with_se > 0U) {
            const double final_se{std::sqrt(global_variance_sum) * inv_num_threads};
            std::cout << " +/- " << final_se;
        } else {
            std::cout << " +/- N/A (insufficient blocks)";
        }

        std::cout << std::endl
                  << "Elapsed: " << elapsed.count() << " s" << std::endl
                  << "Acceptance Rate: " << final_acceptance_rate << "%\n" << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cout << "Exception: " << e.what() << std::endl;
        std::exit(-1);
    }  

    return 0;
}