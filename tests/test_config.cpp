#include "test_utilities.hpp"

#include "config/config.hpp"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>

namespace {

/// Writes a temporary config file and returns a Config parsed from it.
Config parse_config_string(const std::string& content) {
    const std::string PATH{"test_config_tmp.cfg"};
    {
        std::ofstream out{PATH};
        out << content;
    }
    return Config::from_file(PATH);
}

} // namespace

TEST_CASE("Config default constructor produces expected defaults and derived values", "[config]") {
    const Config cfg{};

    REQUIRE(cfg.num_threads == 1U);
    REQUIRE(cfg.num_particles == 7U);
    REQUIRE(cfg.warmup_sweeps == 100U);
    REQUIRE(cfg.measure_sweeps == 100U);
    REQUIRE(cfg.box_length == 10.0);
    REQUIRE(cfg.block_size == 100U);
    REQUIRE(cfg.master_seed == 42U);
    REQUIRE(cfg.is_master_thread == false);

    // Derived fields:
    REQUIRE(cfg.warmup_steps == 7U * 100U);
    REQUIRE(cfg.measure_steps == 7U * 100U);
    require_near(cfg.step_size, 1.0); // 10.0 / 10.0
}

TEST_CASE("Config::from_file parses all recognized keys", "[config]") {
    const Config cfg{parse_config_string(
        "Num_Threads = 4\n"
        "Num_Particles = 19\n"
        "Warmup_Sweeps = 200\n"
        "Measure_Sweeps = 500\n"
        "Box_Length = 8.5\n"
        "Block_Size = 50\n"
        "Master_Seed = 99999\n"
        "Is_Master_Thread = true\n")};

    REQUIRE(cfg.num_threads == 4U);
    REQUIRE(cfg.num_particles == 19U);
    REQUIRE(cfg.warmup_sweeps == 200U);
    REQUIRE(cfg.measure_sweeps == 500U);
    require_near(cfg.box_length, 8.5);
    REQUIRE(cfg.block_size == 50U);
    REQUIRE(cfg.master_seed == 99999U);
    REQUIRE(cfg.is_master_thread == true);

    // Derived:
    REQUIRE(cfg.warmup_steps == 19U * 200U);
    REQUIRE(cfg.measure_steps == 19U * 500U);
    require_near(cfg.step_size, 8.5 / 10.0);
}

TEST_CASE("Config::from_file ignores comments and blank lines", "[config]") {
    const Config cfg{parse_config_string(
        "# This is a comment\n"
        "\n"
        "Num_Particles = 33  # inline comment\n"
        "   \n"
        "Box_Length = 12.0\n")};

    REQUIRE(cfg.num_particles == 33U);
    require_near(cfg.box_length, 12.0);
    // Other fields remain at defaults:
    REQUIRE(cfg.warmup_sweeps == 100U);
}

TEST_CASE("Config::from_file handles whitespace around keys and values", "[config]") {
    const Config cfg{parse_config_string(
        "  Num_Particles  =  27  \n"
        "  Box_Length  =  7.5  \n")};

    REQUIRE(cfg.num_particles == 27U);
    require_near(cfg.box_length, 7.5);
}

TEST_CASE("Config::from_file warns on unknown keys and preserves known values", "[config]") {
    const std::string output{capture_stdout([] {
        parse_config_string(
            "Num_Particles = 7\n"
            "Unknown_Key = 42\n"
            "Another_Bad = hello\n");
    })};

    REQUIRE(output.find("Warning") != std::string::npos);
    REQUIRE(output.find("unknown keys ignored") != std::string::npos);
    REQUIRE(output.find("Unknown_Key") != std::string::npos);
}

TEST_CASE("Config::from_file uses defaults when file does not exist", "[config]") {
    const std::string output{capture_stdout([] {
        const Config cfg{Config::from_file("nonexistent_file.cfg")};
        // Should still produce valid defaults:
        REQUIRE(cfg.num_particles == 7U);
        REQUIRE(cfg.warmup_sweeps == 100U);
    })};

    REQUIRE(output.find("File not found") != std::string::npos);
}

TEST_CASE("Config::from_file skips lines without an equals sign", "[config]") {
    const Config cfg{parse_config_string(
        "Num_Particles = 19\n"
        "this line has no equals\n"
        "Box_Length = 5.0\n")};

    REQUIRE(cfg.num_particles == 19U);
    require_near(cfg.box_length, 5.0);
}

TEST_CASE("Config::from_file parses Is_Master_Thread boolean variants", "[config]") {
    {
        const Config cfg{parse_config_string("Is_Master_Thread = true\n")};
        REQUIRE(cfg.is_master_thread == true);
    }
    {
        const Config cfg{parse_config_string("Is_Master_Thread = TRUE\n")};
        REQUIRE(cfg.is_master_thread == true);
    }
    {
        const Config cfg{parse_config_string("Is_Master_Thread = 1\n")};
        REQUIRE(cfg.is_master_thread == true);
    }
    {
        const Config cfg{parse_config_string("Is_Master_Thread = false\n")};
        REQUIRE(cfg.is_master_thread == false);
    }
    {
        const Config cfg{parse_config_string("Is_Master_Thread = 0\n")};
        REQUIRE(cfg.is_master_thread == false);
    }
}

TEST_CASE("Config derived fields recompute correctly for various particle counts", "[config]") {
    {
        const Config cfg{parse_config_string(
            "Num_Particles = 1\n"
            "Warmup_Sweeps = 50\n"
            "Measure_Sweeps = 200\n"
            "Box_Length = 20.0\n")};

        REQUIRE(cfg.warmup_steps == 1U * 50U);
        REQUIRE(cfg.measure_steps == 1U * 200U);
        require_near(cfg.step_size, 2.0);
    }
    {
        const Config cfg{parse_config_string(
            "Num_Particles = 100\n"
            "Warmup_Sweeps = 10\n"
            "Measure_Sweeps = 30\n"
            "Box_Length = 5.0\n")};

        REQUIRE(cfg.warmup_steps == 100U * 10U);
        REQUIRE(cfg.measure_steps == 100U * 30U);
        require_near(cfg.step_size, 0.5);
    }
}

TEST_CASE("Config::from_file handles empty config file gracefully", "[config]") {
    const Config cfg{parse_config_string("")};

    // All defaults:
    REQUIRE(cfg.num_particles == 7U);
    REQUIRE(cfg.warmup_sweeps == 100U);
    REQUIRE(cfg.measure_sweeps == 100U);
    require_near(cfg.box_length, 10.0);
    REQUIRE(cfg.block_size == 100U);
    REQUIRE(cfg.master_seed == 42U);
}

TEST_CASE("Config::from_file partial override leaves other fields at defaults", "[config]") {
    const Config cfg{parse_config_string("Block_Size = 7\n")};

    REQUIRE(cfg.block_size == 7U);
    // Everything else is default:
    REQUIRE(cfg.num_particles == 7U);
    REQUIRE(cfg.warmup_sweeps == 100U);
    REQUIRE(cfg.measure_sweeps == 100U);
    require_near(cfg.box_length, 10.0);
    REQUIRE(cfg.master_seed == 42U);
}