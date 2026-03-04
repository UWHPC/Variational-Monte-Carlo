#pragma once

#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

struct Config {
    std::size_t numParticles; // Number of particles
    double boxLength;         // Length of box (grid)
    std::size_t warmupSteps;  // Warm up
    std::size_t measureSteps; // M
    double stepSize;          // proposal half-width s
    uint64_t seed;            // Random seed
    std::size_t blockSize;    // Size of block
};

namespace {

struct HelpRequested final : std::exception {
    [[nodiscard]] const char* what() const noexcept override { return "Help requested"; }
};

enum class OptionId {
    numParticles,
    boxLength,
    warmupSteps,
    measureSteps,
    stepSize,
    seed,
    blockSize,
};

struct RawConfigValues {
    std::optional<std::string_view> numParticles;
    std::optional<std::string_view> boxLength;
    std::optional<std::string_view> warmupSteps;
    std::optional<std::string_view> measureSteps;
    std::optional<std::string_view> stepSize;
    std::optional<std::string_view> seed;
    std::optional<std::string_view> blockSize;
};

[[nodiscard]] std::optional<OptionId> parseOptionName(std::string_view name) {
    if (name == "numParticles" || name == "num-particles") {
        return OptionId::numParticles;
    }
    if (name == "boxLength" || name == "box-length") {
        return OptionId::boxLength;
    }
    if (name == "warmupSteps" || name == "warmup-steps") {
        return OptionId::warmupSteps;
    }
    if (name == "measureSteps" || name == "measure-steps") {
        return OptionId::measureSteps;
    }
    if (name == "stepSize" || name == "step-size") {
        return OptionId::stepSize;
    }
    if (name == "seed") {
        return OptionId::seed;
    }
    if (name == "blockSize" || name == "block-size") {
        return OptionId::blockSize;
    }
    return std::nullopt;
}

[[nodiscard]] std::string_view canonicalOptionName(OptionId id) {
    switch (id) {
    case OptionId::numParticles:
        return "numParticles";
    case OptionId::boxLength:
        return "boxLength";
    case OptionId::warmupSteps:
        return "warmupSteps";
    case OptionId::measureSteps:
        return "measureSteps";
    case OptionId::stepSize:
        return "stepSize";
    case OptionId::seed:
        return "seed";
    case OptionId::blockSize:
        return "blockSize";
    }
    throw std::logic_error{"Unknown option id"};
}

[[nodiscard]] std::optional<std::string_view> getOptionValue(const RawConfigValues& raw, OptionId id) {
    switch (id) {
    case OptionId::numParticles:
        return raw.numParticles;
    case OptionId::boxLength:
        return raw.boxLength;
    case OptionId::warmupSteps:
        return raw.warmupSteps;
    case OptionId::measureSteps:
        return raw.measureSteps;
    case OptionId::stepSize:
        return raw.stepSize;
    case OptionId::seed:
        return raw.seed;
    case OptionId::blockSize:
        return raw.blockSize;
    }
    throw std::logic_error{"Unknown option id"};
}

void setOptionValue(RawConfigValues& raw, OptionId id, std::string_view value) {
    switch (id) {
    case OptionId::numParticles:
        raw.numParticles = value;
        return;
    case OptionId::boxLength:
        raw.boxLength = value;
        return;
    case OptionId::warmupSteps:
        raw.warmupSteps = value;
        return;
    case OptionId::measureSteps:
        raw.measureSteps = value;
        return;
    case OptionId::stepSize:
        raw.stepSize = value;
        return;
    case OptionId::seed:
        raw.seed = value;
        return;
    case OptionId::blockSize:
        raw.blockSize = value;
        return;
    }
    throw std::logic_error{"Unknown option id"};
}

template <typename T> [[nodiscard]] T parseInteger(std::string_view text, std::string_view optionName) {
    T value{};
    const char* const begin{text.data()};
    const char* const end{text.data() + text.size()};
    const auto [ptr, ec]{std::from_chars(begin, end, value)};
    if (ec != std::errc{} || ptr != end) {
        throw std::invalid_argument{"Invalid integer value for --" + std::string{optionName} + ": '" +
                                    std::string{text} + "'"};
    }
    return value;
}

[[nodiscard]] double parseDouble(std::string_view text, std::string_view optionName) {
    double value{};
    const char* const begin{text.data()};
    const char* const end{text.data() + text.size()};
    const auto [ptr, ec]{std::from_chars(begin, end, value)};
    if (ec != std::errc{} || ptr != end) {
        throw std::invalid_argument{"Invalid floating-point value for --" + std::string{optionName} + ": '" +
                                    std::string{text} + "'"};
    }
    return value;
}

void validateConfig(const Config& config) {
    if (config.numParticles < 1U) {
        throw std::invalid_argument{"--numParticles must be >= 1"};
    }
    if (!std::isfinite(config.boxLength) || config.boxLength <= 0.0) {
        throw std::invalid_argument{"--boxLength must be finite and > 0"};
    }
    if (!std::isfinite(config.stepSize) || config.stepSize <= 0.0) {
        throw std::invalid_argument{"--stepSize must be finite and > 0"};
    }
    if (config.measureSteps < 1U) {
        throw std::invalid_argument{"--measureSteps must be >= 1"};
    }
    if (config.blockSize < 1U) {
        throw std::invalid_argument{"--blockSize must be >= 1"};
    }
}

void printUsage(const char* programName) {
    std::cout
        << "Usage:\n"
        << "  " << programName
        << " --numParticles N --boxLength L --warmupSteps W --measureSteps M --stepSize S --seed R --blockSize B\n\n"
        << "Aliases:\n"
        << "  --num-particles, --box-length, --warmup-steps, --measure-steps, --step-size, --block-size\n";
}

[[nodiscard]] Config parseArgs(int argc, char** argv) {
    if (argc <= 1) {
        throw std::invalid_argument{"No arguments provided"};
    }

    RawConfigValues raw{};
    for (int idx{1}; idx < argc; ++idx) {
        std::string_view token{argv[idx]};

        if (token.empty() || token[0] != '-') {
            throw std::invalid_argument{"Unexpected positional argument: '" + std::string{token} + "'"};
        }
        if (token == "--help" || token == "-h") {
            throw HelpRequested{};
        }
        if (!token.starts_with("--")) {
            throw std::invalid_argument{"Unsupported short option: '" + std::string{token} + "'"};
        }

        token.remove_prefix(2);
        if (token.empty()) {
            throw std::invalid_argument{"Invalid empty option name"};
        }

        std::string_view optionName{};
        std::string_view optionValue{};
        const std::size_t equalsPos{token.find('=')};
        if (equalsPos == std::string_view::npos) {
            optionName = token;
            if (idx + 1 >= argc) {
                throw std::invalid_argument{"Missing value for option --" + std::string{optionName}};
            }
            ++idx;
            optionValue = argv[idx];
        } else {
            optionName = token.substr(0, equalsPos);
            optionValue = token.substr(equalsPos + 1);
        }

        const auto parsedOption{parseOptionName(optionName)};
        if (!parsedOption.has_value()) {
            throw std::invalid_argument{"Unknown option: --" + std::string{optionName}};
        }
        if (optionValue.empty()) {
            throw std::invalid_argument{"Missing value for option --" +
                                        std::string{canonicalOptionName(*parsedOption)}};
        }
        if (getOptionValue(raw, *parsedOption).has_value()) {
            throw std::invalid_argument{"Duplicate option: --" + std::string{canonicalOptionName(*parsedOption)}};
        }

        setOptionValue(raw, *parsedOption, optionValue);
    }

    constexpr std::array<OptionId, 7> allOptions{
        OptionId::numParticles, OptionId::boxLength, OptionId::warmupSteps, OptionId::measureSteps,
        OptionId::stepSize,     OptionId::seed,      OptionId::blockSize,
    };

    std::string missing{};
    for (const OptionId option : allOptions) {
        if (!getOptionValue(raw, option).has_value()) {
            if (!missing.empty()) {
                missing += ", ";
            }
            missing += "--";
            missing += canonicalOptionName(option);
        }
    }
    if (!missing.empty()) {
        throw std::invalid_argument{"Missing required options: " + missing};
    }

    const Config config{
        .numParticles = parseInteger<std::size_t>(*raw.numParticles, "numParticles"),
        .boxLength = parseDouble(*raw.boxLength, "boxLength"),
        .warmupSteps = parseInteger<std::size_t>(*raw.warmupSteps, "warmupSteps"),
        .measureSteps = parseInteger<std::size_t>(*raw.measureSteps, "measureSteps"),
        .stepSize = parseDouble(*raw.stepSize, "stepSize"),
        .seed = parseInteger<std::uint64_t>(*raw.seed, "seed"),
        .blockSize = parseInteger<std::size_t>(*raw.blockSize, "blockSize"),
    };
    validateConfig(config);
    return config;
}

void printConfig(const Config& config) {
    std::cout << "numParticles: " << config.numParticles << '\n'
              << "boxLength: " << config.boxLength << '\n'
              << "warmupSteps: " << config.warmupSteps << '\n'
              << "measureSteps: " << config.measureSteps << '\n'
              << "stepSize: " << config.stepSize << '\n'
              << "seed: " << config.seed << '\n'
              << "blockSize: " << config.blockSize << '\n';
}

} // Namespace
