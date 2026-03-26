#include "test_utilities.hpp"

#include "config/config.hpp"

#include <initializer_list>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

class ArgvBuilder {
public:
    explicit ArgvBuilder(std::initializer_list<std::string_view> args) {
        storage_.reserve(args.size());
        for (const std::string_view arg : args) {
            storage_.emplace_back(arg);
        }

        argv_.reserve(storage_.size());
        for (std::string& arg : storage_) {
            argv_.push_back(arg.data());
        }
    }

    [[nodiscard]] int argc() const { return static_cast<int>(argv_.size()); }
    [[nodiscard]] char** argv() { return argv_.data(); }

private:
    std::vector<std::string> storage_{};
    std::vector<char*> argv_{};
};

Config parse_config(std::initializer_list<std::string_view> args) {
    ArgvBuilder argv{args};
    return parse_args(argv.argc(), argv.argv());
}

std::string invalidArgumentMessage(std::initializer_list<std::string_view> args) {
    ArgvBuilder argv{args};
    try {
        (void)parse_args(argv.argc(), argv.argv());
    } catch (const std::invalid_argument& ex) {
        return ex.what();
    }

    FAIL("Expected parse_args to throw std::invalid_argument");
    return {};
}

} // namespace

TEST_CASE("parse_args parses canonical options with mixed token styles", "[config]") {
    const Config cfg{
        parse_config({"vmc", "--num_particles", "16", "--box_length=4.5", "--warmup_steps", "50", "--measure_steps=200",
                     "--step_size", "0.25", "--seed", "12345", "--block_size", "20"})};

    REQUIRE(cfg.num_particles == 16U);
    REQUIRE(cfg.box_length == 4.5);
    REQUIRE(cfg.warmup_steps == 50U);
    REQUIRE(cfg.measure_steps == 200U);
    REQUIRE(cfg.step_size == 0.25);
    REQUIRE(cfg.seed == 12345U);
    REQUIRE(cfg.block_size == 20U);
}

TEST_CASE("parse_args parses alias options", "[config]") {
    const Config cfg{
        parse_config({"vmc", "--num-particles", "8", "--box-length", "9", "--warmup-steps", "10", "--measure-steps",
                     "100", "--step-size", "0.1", "--seed", "42", "--block-size", "5"})};

    REQUIRE(cfg.num_particles == 8U);
    REQUIRE(cfg.box_length == 9.0);
    REQUIRE(cfg.warmup_steps == 10U);
    REQUIRE(cfg.measure_steps == 100U);
    REQUIRE(cfg.step_size == 0.1);
    REQUIRE(cfg.seed == 42U);
    REQUIRE(cfg.block_size == 5U);
}

TEST_CASE("parse_args rejects missing required options", "[config]") {
    const std::string message{
        invalidArgumentMessage({"vmc", "--num_particles", "16", "--box_length", "4", "--warmup_steps", "10",
                                "--measure_steps", "100", "--step_size", "0.25", "--seed", "1"})};

    REQUIRE(message.find("Missing required options") != std::string::npos);
    REQUIRE(message.find("--block_size") != std::string::npos);
}

TEST_CASE("parse_args rejects duplicate aliases for the same option", "[config]") {
    const std::string message{invalidArgumentMessage(
        {"vmc", "--num_particles", "16", "--num-particles", "16", "--box_length", "4", "--warmup_steps", "10",
         "--measure_steps", "100", "--step_size", "0.25", "--seed", "1", "--block_size", "4"})};

    REQUIRE(message.find("Duplicate option") != std::string::npos);
    REQUIRE(message.find("--num_particles") != std::string::npos);
}

TEST_CASE("parse_args rejects unknown options", "[config]") {
    const std::string message{invalidArgumentMessage({"vmc", "--num_particles", "16", "--box_length", "4",
                                                      "--warmup_steps", "10", "--measure_steps", "100", "--step_size",
                                                      "0.25", "--seed", "1", "--block_size", "4", "--badFlag", "1"})};

    REQUIRE(message.find("Unknown option") != std::string::npos);
    REQUIRE(message.find("--badFlag") != std::string::npos);
}

TEST_CASE("parse_args rejects invalid numeric values", "[config]") {
    const std::string message{
        invalidArgumentMessage({"vmc", "--num_particles", "16", "--box_length", "4", "--warmup_steps", "10",
                                "--measure_steps", "100", "--step_size", "oops", "--seed", "1", "--block_size", "4"})};

    REQUIRE(message.find("Invalid floating-point value") != std::string::npos);
    REQUIRE(message.find("--step_size") != std::string::npos);
}

TEST_CASE("parse_args rejects unsupported short options", "[config]") {
    const std::string message{
        invalidArgumentMessage({"vmc", "-n", "16", "--box_length", "4", "--warmup_steps", "10", "--measure_steps",
                                "100", "--step_size", "0.25", "--seed", "1", "--block_size", "4"})};

    REQUIRE(message.find("Unsupported short option") != std::string::npos);
    REQUIRE(message.find("-n") != std::string::npos);
}

TEST_CASE("parse_args throws HelpRequested for --help and -h", "[config]") {
    {
        ArgvBuilder argv{"vmc", "--help"};
        REQUIRE_THROWS_AS(parse_args(argv.argc(), argv.argv()), HelpRequested);
    }

    {
        ArgvBuilder argv{"vmc", "-h"};
        REQUIRE_THROWS_AS(parse_args(argv.argc(), argv.argv()), HelpRequested);
    }
}

TEST_CASE("parse_args rejects missing option values", "[config]") {
    {
        const std::string message{invalidArgumentMessage({"vmc", "--num_particles"})};
        REQUIRE(message.find("Missing value for option --num_particles") != std::string::npos);
    }

    {
        const std::string message{invalidArgumentMessage({"vmc", "--num_particles=", "--box_length", "4"})};
        REQUIRE(message.find("Missing value for option --num_particles") != std::string::npos);
    }
}

TEST_CASE("parse_args rejects unexpected positional arguments and empty option names", "[config]") {
    {
        const std::string message{invalidArgumentMessage(
            {"vmc", "--num_particles", "16", "oops", "--box_length", "4", "--warmup_steps", "10"})};
        REQUIRE(message.find("Unexpected positional argument") != std::string::npos);
        REQUIRE(message.find("oops") != std::string::npos);
    }

    {
        const std::string message{invalidArgumentMessage({"vmc", "--", "1"})};
        REQUIRE(message.find("Invalid empty option name") != std::string::npos);
    }
}

TEST_CASE("parse_args rejects empty invocation and invalid integer values", "[config]") {
    {
        const std::string message{invalidArgumentMessage({"vmc"})};
        REQUIRE(message.find("No arguments provided") != std::string::npos);
    }

    {
        const std::string message{invalidArgumentMessage({"vmc", "--num_particles", "1.5", "--box_length", "4",
                                                          "--warmup_steps", "10", "--measure_steps", "100",
                                                          "--step_size", "0.25", "--seed", "1", "--block_size", "4"})};
        REQUIRE(message.find("Invalid integer value") != std::string::npos);
        REQUIRE(message.find("--num_particles") != std::string::npos);
    }
}

TEST_CASE("parse_args enforces runtime numeric constraints", "[config]") {
    {
        const std::string message{invalidArgumentMessage({"vmc", "--num_particles", "0", "--box_length", "4",
                                                          "--warmup_steps", "10", "--measure_steps", "100",
                                                          "--step_size", "0.25", "--seed", "1", "--block_size", "4"})};
        REQUIRE(message.find("--num_particles must be >= 1") != std::string::npos);
    }

    {
        const std::string message{invalidArgumentMessage({"vmc", "--num_particles", "4", "--box_length", "0",
                                                          "--warmup_steps", "10", "--measure_steps", "100",
                                                          "--step_size", "0.25", "--seed", "1", "--block_size", "4"})};
        REQUIRE(message.find("--box_length must be finite and > 0") != std::string::npos);
    }

    {
        const std::string message{invalidArgumentMessage({"vmc", "--num_particles", "4", "--box_length", "4",
                                                          "--warmup_steps", "10", "--measure_steps", "0", "--step_size",
                                                          "0.25", "--seed", "1", "--block_size", "4"})};
        REQUIRE(message.find("--measure_steps must be >= 1") != std::string::npos);
    }

    {
        const std::string message{
            invalidArgumentMessage({"vmc", "--num_particles", "4", "--box_length", "4", "--warmup_steps", "10",
                                    "--measure_steps", "100", "--step_size", "0", "--seed", "1", "--block_size", "4"})};
        REQUIRE(message.find("--step_size must be finite and > 0") != std::string::npos);
    }

    {
        const std::string message{invalidArgumentMessage({"vmc", "--num_particles", "4", "--box_length", "4",
                                                          "--warmup_steps", "10", "--measure_steps", "100",
                                                          "--step_size", "0.25", "--seed", "1", "--block_size", "0"})};
        REQUIRE(message.find("--block_size must be >= 1") != std::string::npos);
    }
}

TEST_CASE("printUsage and printConfig include all expected fields", "[config]") {
    const std::string usage{capture_stdout([] { print_usage("vmc-test"); })};
    REQUIRE(usage.find("Usage:") != std::string::npos);
    REQUIRE(usage.find("vmc-test") != std::string::npos);
    REQUIRE(usage.find("--num_particles") != std::string::npos);
    REQUIRE(usage.find("--block_size") != std::string::npos);
    REQUIRE(usage.find("Aliases:") != std::string::npos);

    const Config cfg{.num_particles = 12U,
                     .box_length = 9.5,
                     .warmup_steps = 30U,
                     .measure_steps = 500U,
                     .step_size = 0.2,
                     .seed = 99U,
                     .block_size = 25U};
    const std::string configText{capture_stdout([&cfg] { print_config(cfg); })};
    REQUIRE(configText.find("num_particles: 12") != std::string::npos);
    REQUIRE(configText.find("box_length: 9.5") != std::string::npos);
    REQUIRE(configText.find("warmup_steps: 30") != std::string::npos);
    REQUIRE(configText.find("measure_steps: 500") != std::string::npos);
    REQUIRE(configText.find("step_size: 0.2") != std::string::npos);
    REQUIRE(configText.find("seed: 99") != std::string::npos);
    REQUIRE(configText.find("block_size: 25") != std::string::npos);
}