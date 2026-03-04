#include <catch2/catch_test_macros.hpp>

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

Config parseConfig(std::initializer_list<std::string_view> args) {
    ArgvBuilder argv{args};
    return parseArgs(argv.argc(), argv.argv());
}

std::string invalidArgumentMessage(std::initializer_list<std::string_view> args) {
    ArgvBuilder argv{args};
    try {
        (void)parseArgs(argv.argc(), argv.argv());
    } catch (const std::invalid_argument& ex) {
        return ex.what();
    }

    FAIL("Expected parseArgs to throw std::invalid_argument");
    return {};
}

} // namespace

TEST_CASE("parseArgs parses canonical options with mixed token styles", "[config]") {
    const Config cfg{parseConfig({"vmc", "--numParticles", "16", "--boxLength=4.5", "--warmupSteps", "50",
                                  "--measureSteps=200", "--stepSize", "0.25", "--seed", "12345", "--blockSize", "20"})};

    REQUIRE(cfg.numParticles == 16U);
    REQUIRE(cfg.boxLength == 4.5);
    REQUIRE(cfg.warmupSteps == 50U);
    REQUIRE(cfg.measureSteps == 200U);
    REQUIRE(cfg.stepSize == 0.25);
    REQUIRE(cfg.seed == 12345U);
    REQUIRE(cfg.blockSize == 20U);
}

TEST_CASE("parseArgs parses alias options", "[config]") {
    const Config cfg{
        parseConfig({"vmc", "--num-particles", "8", "--box-length", "9", "--warmup-steps", "10", "--measure-steps",
                     "100", "--step-size", "0.1", "--seed", "42", "--block-size", "5"})};

    REQUIRE(cfg.numParticles == 8U);
    REQUIRE(cfg.boxLength == 9.0);
    REQUIRE(cfg.warmupSteps == 10U);
    REQUIRE(cfg.measureSteps == 100U);
    REQUIRE(cfg.stepSize == 0.1);
    REQUIRE(cfg.seed == 42U);
    REQUIRE(cfg.blockSize == 5U);
}

TEST_CASE("parseArgs rejects missing required options", "[config]") {
    const std::string message{
        invalidArgumentMessage({"vmc", "--numParticles", "16", "--boxLength", "4", "--warmupSteps", "10",
                                "--measureSteps", "100", "--stepSize", "0.25", "--seed", "1"})};

    REQUIRE(message.find("Missing required options") != std::string::npos);
    REQUIRE(message.find("--blockSize") != std::string::npos);
}

TEST_CASE("parseArgs rejects duplicate aliases for the same option", "[config]") {
    const std::string message{invalidArgumentMessage({"vmc", "--numParticles", "16", "--num-particles", "16",
                                                      "--boxLength", "4", "--warmupSteps", "10", "--measureSteps",
                                                      "100", "--stepSize", "0.25", "--seed", "1", "--blockSize", "4"})};

    REQUIRE(message.find("Duplicate option") != std::string::npos);
    REQUIRE(message.find("--numParticles") != std::string::npos);
}

TEST_CASE("parseArgs rejects unknown options", "[config]") {
    const std::string message{invalidArgumentMessage({"vmc", "--numParticles", "16", "--boxLength", "4",
                                                      "--warmupSteps", "10", "--measureSteps", "100", "--stepSize",
                                                      "0.25", "--seed", "1", "--blockSize", "4", "--badFlag", "1"})};

    REQUIRE(message.find("Unknown option") != std::string::npos);
    REQUIRE(message.find("--badFlag") != std::string::npos);
}

TEST_CASE("parseArgs rejects invalid numeric values", "[config]") {
    const std::string message{
        invalidArgumentMessage({"vmc", "--numParticles", "16", "--boxLength", "4", "--warmupSteps", "10",
                                "--measureSteps", "100", "--stepSize", "oops", "--seed", "1", "--blockSize", "4"})};

    REQUIRE(message.find("Invalid floating-point value") != std::string::npos);
    REQUIRE(message.find("--stepSize") != std::string::npos);
}

TEST_CASE("parseArgs rejects unsupported short options", "[config]") {
    const std::string message{
        invalidArgumentMessage({"vmc", "-n", "16", "--boxLength", "4", "--warmupSteps", "10", "--measureSteps", "100",
                                "--stepSize", "0.25", "--seed", "1", "--blockSize", "4"})};

    REQUIRE(message.find("Unsupported short option") != std::string::npos);
    REQUIRE(message.find("-n") != std::string::npos);
}

TEST_CASE("parseArgs throws HelpRequested for --help and -h", "[config]") {
    {
        ArgvBuilder argv{"vmc", "--help"};
        REQUIRE_THROWS_AS(parseArgs(argv.argc(), argv.argv()), HelpRequested);
    }

    {
        ArgvBuilder argv{"vmc", "-h"};
        REQUIRE_THROWS_AS(parseArgs(argv.argc(), argv.argv()), HelpRequested);
    }
}
