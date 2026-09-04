#pragma once

#include "../support/checks.hpp"

#include "config/config.hpp"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace test {

class ConfigFile {
private:
  std::filesystem::path path_;

public:
  explicit ConfigFile(const std::string& content)
  : path_{std::filesystem::temp_directory_path() / "vmc-test-config.cfg"} {
    std::ofstream output{path_};
    output << content;
  }

  ~ConfigFile() {
    std::error_code error{};
    std::filesystem::remove(path_, error);
  }

  ConfigFile(const ConfigFile&) = delete;
  ConfigFile& operator=(const ConfigFile&) = delete;

  [[nodiscard]] Config parse() const {
    return Config::from_file(path_.string());
  }
};

} // namespace test

TEST_CASE("Configuration defaults form a valid simulation", "[config]") {
  const Config config{};

  REQUIRE(config.num_particles == 7uz);
  REQUIRE(config.num_walkers == 1uz);
  REQUIRE(config.warmup_sweeps == 100uz);
  REQUIRE(config.measure_sweeps == 100uz);
  REQUIRE(config.block_size == 100uz);
  REQUIRE(config.master_seed == 42uz);
  require_near(config.box_length, 10.0_fp);
  require_near(config.step_size, 1.0_fp);
  REQUIRE(config.warmup_steps == 700uz);
  REQUIRE(config.measure_steps == 700uz);
}

TEST_CASE("Configuration parses documented values and derives total work", "[config]") {
  const test::ConfigFile file{
    "Num_Particles = 19\n"
    "Num_Walkers = 4\n"
    "Warmup_Sweeps = 20\n"
    "Measure_Sweeps = 50\n"
    "Box_Length = 8.5\n"
    "Block_Size = 5\n"
    "Master_Seed = 999\n"
    "Jastrow_A = 0.3\n"
    "Jastrow_B = 0.7\n"
  };
  const Config config{file.parse()};

  REQUIRE(config.num_particles == 19uz);
  REQUIRE(config.num_walkers == 4uz);
  REQUIRE(config.warmup_steps == 1520uz);
  REQUIRE(config.measure_steps == 3800uz);
  REQUIRE(config.block_size == 5uz);
  REQUIRE(config.master_seed == 999uz);
  require_near(config.box_length, 8.5_fp);
  require_near(config.step_size, 0.85_fp);
  require_near(config.jastrow_a, 0.3_fp);
  require_near(config.jastrow_b, 0.7_fp);
}

TEST_CASE("Configuration ignores comments, whitespace, and unknown keys", "[config]") {
  const test::ConfigFile file{
    "# comment\n"
    "  Num_Particles = 3  # inline comment\n"
    "Unknown_Key = ignored\n"
    "\n"
    "Measure_Sweeps = 4\n"
  };
  const Config config{file.parse()};

  REQUIRE(config.num_particles == 3uz);
  REQUIRE(config.measure_sweeps == 4uz);
  REQUIRE(config.num_walkers == 1uz);
}

TEST_CASE("Configuration rejects invalid physical and sampling values", "[config]") {
  for (const auto& content : {
    "Num_Particles = 0\n",
    "Num_Walkers = 0\n",
    "Measure_Sweeps = 0\n",
    "Block_Size = 0\n",
    "Box_Length = 0\n",
    "Box_Length = -1\n"
  }) {
    const test::ConfigFile file{content};
    REQUIRE_THROWS_AS(file.parse(), std::invalid_argument);
  }
}
