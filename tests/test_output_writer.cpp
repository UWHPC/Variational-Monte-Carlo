#include <catch2/catch_test_macros.hpp>

#include "output_writer/output_writer.hpp"

#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

TEST_CASE("JsonOutputWriter writes init/frame/done records and escapes JSON strings", "[output_writer]") {
    std::ostringstream out{};
    JsonOutputWriter writer{out};

    const InitData init{.run_id = std::string{"run\"\\\b\f\n\r\t\x01"},
                        .num_particles = 7U,
                        .box_length = 6.5,
                        .warmup_steps = 10U,
                        .measure_steps = 20U,
                        .step_size = 0.2,
                        .seed = 42U,
                        .block_size = 5U};
    writer.write_init(init);

    writer.write_frame(FrameData{.step = 1U,
                                 .accepted = 0U,
                                 .proposed = 1U,
                                 .acceptance_rate = 0.0,
                                 .local_energy = -1.25,
                                 .mean_energy = -1.25,
                                 .standard_error = std::nullopt,
                                 .positions = std::vector<double>{1.0, 2.0, 3.0}});

    writer.write_frame(FrameData{.step = 2U,
                                 .accepted = 1U,
                                 .proposed = 2U,
                                 .acceptance_rate = 0.5,
                                 .local_energy = -1.0,
                                 .mean_energy = -1.125,
                                 .standard_error = 0.25,
                                 .positions = std::vector<double>{4.0, 5.0, 6.0}});

    writer.write_done(DoneData{.total_accepted = 1U,
                               .total_proposed = 2U,
                               .final_acceptance_rate = 0.5,
                               .final_mean_energy = -1.125,
                               .final_standard_error = 0.125});
    writer.write_done(DoneData{.total_accepted = 1U,
                               .total_proposed = 2U,
                               .final_acceptance_rate = 0.5,
                               .final_mean_energy = -1.125,
                               .final_standard_error = std::nullopt});

    const std::string text{out.str()};

    REQUIRE(text.find("\"type\":\"init\"") != std::string::npos);
    REQUIRE(text.find("\"runId\":\"run\\\"\\\\\\b\\f\\n\\r\\t\\u0001\"") != std::string::npos);
    REQUIRE(text.find("\"numParticles\":7") != std::string::npos);
    REQUIRE(text.find("\"seed\":42") != std::string::npos);

    REQUIRE(text.find("\"type\":\"frame\"") != std::string::npos);
    REQUIRE(text.find("\"standardErrorAvailable\":false") != std::string::npos);
    REQUIRE(text.find("\"standardErrorAvailable\":true") != std::string::npos);
    REQUIRE(text.find("\"positions\":[1,2,3]") != std::string::npos);
    REQUIRE(text.find("\"positions\":[4,5,6]") != std::string::npos);

    REQUIRE(text.find("\"type\":\"done\"") != std::string::npos);
    REQUIRE(text.find("\"finalStandardErrorAvailable\":false") != std::string::npos);
    REQUIRE(text.find("\"finalStandardErrorAvailable\":true") != std::string::npos);
}

TEST_CASE("CsvOutputWriter methods throw until implemented", "[output_writer]") {
    std::ostringstream out{};
    CsvOutputWriter writer{out};

    REQUIRE_THROWS_AS(
        writer.write_init(InitData{.run_id = "run",
                                   .num_particles = 1U,
                                   .box_length = 1.0,
                                   .warmup_steps = 0U,
                                   .measure_steps = 1U,
                                   .step_size = 0.1,
                                   .seed = 0U,
                                   .block_size = 1U}),
        std::logic_error);

    REQUIRE_THROWS_AS(writer.write_frame(FrameData{.step = 1U,
                                                   .accepted = 1U,
                                                   .proposed = 1U,
                                                   .acceptance_rate = 1.0,
                                                   .local_energy = 0.0,
                                                   .mean_energy = 0.0,
                                                   .standard_error = std::nullopt,
                                                   .positions = {0.0, 0.0, 0.0}}),
                      std::logic_error);

    REQUIRE_THROWS_AS(writer.write_done(DoneData{.total_accepted = 1U,
                                                 .total_proposed = 1U,
                                                 .final_acceptance_rate = 1.0,
                                                 .final_mean_energy = 0.0,
                                                 .final_standard_error = std::nullopt}),
                      std::logic_error);
}

TEST_CASE("make_output_writer selects concrete writer type and validates format", "[output_writer]") {
    std::ostringstream out{};

    std::unique_ptr<OutputWriter> json{make_output_writer(OutputFormat::JSON, out)};
    REQUIRE(dynamic_cast<JsonOutputWriter*>(json.get()) != nullptr);

    std::unique_ptr<OutputWriter> csv{make_output_writer(OutputFormat::CSV, out)};
    REQUIRE(dynamic_cast<CsvOutputWriter*>(csv.get()) != nullptr);

    REQUIRE_THROWS_AS(make_output_writer(static_cast<OutputFormat>(999), out), std::invalid_argument);
}
