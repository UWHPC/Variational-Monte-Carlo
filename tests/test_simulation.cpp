#include "test_utilities.hpp"

#include "simulation/simulation.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

TEST_CASE("Simulation API is present", "[simulation]") {
    REQUIRE(std::is_constructible_v<Simulation, Config>);
    REQUIRE(std::is_same_v<decltype(std::declval<Simulation&>().run()), Simulation::MeasurementSummary>);
}

TEST_CASE("Simulation emits consistent init/frame/done records through OutputWriter", "[simulation]") {
    const Config config{make_config(7U, 6.0, 0U, 4U, 0.35, 12345U, 2U)};

    auto writer{std::make_unique<RecordingOutputWriter>()};
    RecordingOutputWriter* const sink{writer.get()};

    Simulation simulation{config, std::move(writer)};
    simulation.run();

    REQUIRE(sink->saw_init);
    REQUIRE(sink->saw_done);
    REQUIRE(sink->init.has_value());
    REQUIRE(sink->done.has_value());
    REQUIRE(sink->frames.size() == config.measure_steps);

    const InitData& init{*sink->init};
    REQUIRE(init.run_id == "vmc-seed-" + std::to_string(config.master_seed));
    REQUIRE(init.num_particles == config.num_particles);
    REQUIRE(init.box_length == config.box_length);
    REQUIRE(init.warmup_steps == config.warmup_steps);
    REQUIRE(init.measure_steps == config.measure_steps);
    REQUIRE(init.step_size == config.step_size);
    REQUIRE(init.seed == config.master_seed);
    REQUIRE(init.block_size == config.block_size);

    for (std::size_t i = 0; i < sink->frames.size(); ++i) {
        const FrameData& frame{sink->frames[i]};
        const std::size_t expected_step{i + 1U};

        REQUIRE(frame.step == expected_step);
        REQUIRE(frame.proposed == expected_step);
        REQUIRE(frame.accepted <= frame.proposed);
        require_near(frame.acceptance_rate,
                                static_cast<double>(frame.accepted) / static_cast<double>(frame.proposed));
        REQUIRE(frame.positions.size() == 3U * config.num_particles);
    }

    REQUIRE_FALSE(sink->frames.front().standard_error.has_value());
    REQUIRE(sink->frames.back().standard_error.has_value());

    const DoneData& done{*sink->done};
    REQUIRE(done.total_proposed == config.measure_steps);
    REQUIRE(done.total_accepted == sink->frames.back().accepted);
    require_near(done.final_acceptance_rate,
                            static_cast<double>(done.total_accepted) / static_cast<double>(done.total_proposed));
    REQUIRE(done.final_standard_error.has_value());
}

TEST_CASE("Simulation produces blocked energy summary when no writer is attached", "[simulation]") {
    const Config config{make_config(1U, 5.0, 0U, 4U, 0.1, 7U, 2U)};

    Simulation simulation{config};
    const auto summary{simulation.run()};

    // 4 measure steps, block_size=2 -> 2 complete blocks -> blocking is ready
    REQUIRE(summary.standard_error.has_value());
    REQUIRE(std::isfinite(summary.mean_energy));
    REQUIRE(summary.acceptance_rate >= 0.0);
    REQUIRE(summary.acceptance_rate <= 1.0);
}

TEST_CASE("Simulation skips blocked energy summary when insufficient blocks are available", "[simulation]") {
    const Config config{make_config(1U, 5.0, 0U, 3U, 0.1, 9U, 10U)};

    Simulation simulation{config};
    const auto summary{simulation.run()};

    // 3 measure steps, block_size=10 -> 0 complete blocks -> no SE available
    REQUIRE_FALSE(summary.standard_error.has_value());
    REQUIRE(std::isfinite(summary.mean_energy));
    REQUIRE(summary.acceptance_rate >= 0.0);
    REQUIRE(summary.acceptance_rate <= 1.0);
}

TEST_CASE("Simulation accepts zero-step proposals and preserves local energy across frames", "[simulation]") {
    const Config config{make_config(7U, 6.0, 0U, 6U, 0.0, 314159U, 3U)};

    auto writer{std::make_unique<RecordingOutputWriter>()};
    RecordingOutputWriter* const sink{writer.get()};

    Simulation simulation{config, std::move(writer)};
    simulation.run();

    REQUIRE(sink->frames.size() == config.measure_steps);
    const double first_energy{sink->frames.front().local_energy};

    for (std::size_t i = 0; i < sink->frames.size(); ++i) {
        const FrameData& frame{sink->frames[i]};
        REQUIRE(frame.accepted == i + 1U);
        REQUIRE(frame.proposed == i + 1U);
        require_near(frame.local_energy, first_energy);
    }
}

TEST_CASE("Simulation performs rejected moves for a large proposal step size", "[simulation]") {
    const Config config{make_config(7U, 6.0, 0U, 100U, 2.0, 20250308U, 10U)};

    auto writer{std::make_unique<RecordingOutputWriter>()};
    RecordingOutputWriter* const sink{writer.get()};

    Simulation simulation{config, std::move(writer)};
    simulation.run();

    REQUIRE_FALSE(sink->frames.empty());
    const std::size_t accepted_final{sink->frames.back().accepted};
    REQUIRE(accepted_final < config.measure_steps);
}

TEST_CASE("Simulation warmup path executes with adaptive proposal updates", "[simulation]") {
    const Config config{make_config(1U, 4.0, 5U, 1U, 0.2, 42U, 1U)};

    auto writer{std::make_unique<RecordingOutputWriter>()};
    RecordingOutputWriter* const sink{writer.get()};

    Simulation simulation{config, std::move(writer)};
    simulation.run();

    REQUIRE(sink->frames.size() == config.measure_steps);
    REQUIRE(sink->done.has_value());
    REQUIRE(sink->done->total_proposed == config.measure_steps);
}