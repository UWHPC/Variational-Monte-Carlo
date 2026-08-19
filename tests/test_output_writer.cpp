#include <catch2/catch_test_macros.hpp>

#include "output_writer/output_writer.hpp"

#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

TEST_CASE("CsvOutputWriter methods throw until implemented", "[output_writer]") {
  std::ostringstream out{};
  CsvOutputWriter writer{out};

  REQUIRE_THROWS_AS(
    writer.write_init(InitData{
      .run_id = "run",
      .num_particles = 1U,
      .box_length = 1.0_fp,
      .warmup_steps = 0U,
      .measure_steps = 1U,
      .step_size = 0.1_fp,
      .seed = 0U,
      .block_size = 1U
    }),
    std::logic_error
  );

  REQUIRE_THROWS_AS(
    writer.write_frame(FrameData{
      .step = 1U,
      .accepted = 1U,
      .proposed = 1U,
      .acceptance_rate = 1.0_fp,
      .local_energy = 0.0_fp,
      .mean_energy = 0.0_fp,
      .standard_error = std::nullopt,
      .positions = {0.0_fp, 0.0_fp, 0.0_fp}
    }),
    std::logic_error
  );

  REQUIRE_THROWS_AS(
    writer.write_done(DoneData{
      .total_accepted = 1U,
      .total_proposed = 1U,
      .final_acceptance_rate = 1.0_fp,
      .final_mean_energy = 0.0_fp,
      .final_standard_error = std::nullopt
    }),
    std::logic_error
  );
}

TEST_CASE("make_output_writer selects concrete writer type and validates format",
          "[output_writer]") {
  std::ostringstream out{};

  std::unique_ptr<OutputWriter> csv{make_output_writer(OutputFormat::CSV, out)};
  REQUIRE(dynamic_cast<CsvOutputWriter*>(csv.get()) != nullptr);

  REQUIRE_THROWS_AS(make_output_writer(static_cast<OutputFormat>(999), out), std::invalid_argument);
}
