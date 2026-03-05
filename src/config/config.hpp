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
    std::size_t num_particles; // Number of particles
    double box_length;         // Length of box (grid)
    std::size_t warmup_steps;  // Warm up
    std::size_t measure_steps; // M
    double step_size;          // proposal half-width s
    uint64_t seed;             // Random seed
    std::size_t block_size;    // Size of block
};

namespace {

struct HelpRequested final : std::exception {
    [[nodiscard]] const char* what() const noexcept override { return "Help requested"; }
};

enum class OptionId {
    NUM_PARTICLES,
    BOX_LENGTH,
    WARMUP_STEPS,
    MEASURE_STEPS,
    STEP_SIZE,
    SEED,
    BLOCK_SIZE,
};

struct RawConfigValues {
    std::optional<std::string_view> num_particles;
    std::optional<std::string_view> box_length;
    std::optional<std::string_view> warmup_steps;
    std::optional<std::string_view> measure_steps;
    std::optional<std::string_view> step_size;
    std::optional<std::string_view> seed;
    std::optional<std::string_view> block_size;
};

[[nodiscard]] std::optional<OptionId> parse_option_name(std::string_view name) {
    if (name == "num_particles" || name == "num-particles") {
        return OptionId::NUM_PARTICLES;
    }
    if (name == "box_length" || name == "box-length") {
        return OptionId::BOX_LENGTH;
    }
    if (name == "warmup_steps" || name == "warmup-steps") {
        return OptionId::WARMUP_STEPS;
    }
    if (name == "measure_steps" || name == "measure-steps") {
        return OptionId::MEASURE_STEPS;
    }
    if (name == "step_size" || name == "step-size") {
        return OptionId::STEP_SIZE;
    }
    if (name == "seed") {
        return OptionId::SEED;
    }
    if (name == "block_size" || name == "block-size") {
        return OptionId::BLOCK_SIZE;
    }
    return std::nullopt;
}

[[nodiscard]] std::string_view canonical_option_name(OptionId id) {
    switch (id) {
    case OptionId::NUM_PARTICLES:
        return "num_particles";
    case OptionId::BOX_LENGTH:
        return "box_length";
    case OptionId::WARMUP_STEPS:
        return "warmup_steps";
    case OptionId::MEASURE_STEPS:
        return "measure_steps";
    case OptionId::STEP_SIZE:
        return "step_size";
    case OptionId::SEED:
        return "seed";
    case OptionId::BLOCK_SIZE:
        return "block_size";
    }
    throw std::logic_error{"Unknown option id"};
}

[[nodiscard]] std::optional<std::string_view> get_option_value(const RawConfigValues& raw, OptionId id) {
    switch (id) {
    case OptionId::NUM_PARTICLES:
        return raw.num_particles;
    case OptionId::BOX_LENGTH:
        return raw.box_length;
    case OptionId::WARMUP_STEPS:
        return raw.warmup_steps;
    case OptionId::MEASURE_STEPS:
        return raw.measure_steps;
    case OptionId::STEP_SIZE:
        return raw.step_size;
    case OptionId::SEED:
        return raw.seed;
    case OptionId::BLOCK_SIZE:
        return raw.block_size;
    }
    throw std::logic_error{"Unknown option id"};
}

void set_option_value(RawConfigValues& raw, OptionId id, std::string_view value) {
    switch (id) {
    case OptionId::NUM_PARTICLES:
        raw.num_particles = value;
        return;
    case OptionId::BOX_LENGTH:
        raw.box_length = value;
        return;
    case OptionId::WARMUP_STEPS:
        raw.warmup_steps = value;
        return;
    case OptionId::MEASURE_STEPS:
        raw.measure_steps = value;
        return;
    case OptionId::STEP_SIZE:
        raw.step_size = value;
        return;
    case OptionId::SEED:
        raw.seed = value;
        return;
    case OptionId::BLOCK_SIZE:
        raw.block_size = value;
        return;
    }
    throw std::logic_error{"Unknown option id"};
}

template <typename T> [[nodiscard]] T parse_integer(std::string_view text, std::string_view option_name) {
    T value{};
    const char* const begin{text.data()};
    const char* const end{text.data() + text.size()};
    const auto [ptr, ec]{std::from_chars(begin, end, value)};
    if (ec != std::errc{} || ptr != end) {
        throw std::invalid_argument{"Invalid integer value for --" + std::string{option_name} + ": '" +
                                    std::string{text} + "'"};
    }
    return value;
}

[[nodiscard]] double parse_double(std::string_view text, std::string_view option_name) {
    double value{};
    const char* const begin{text.data()};
    const char* const end{text.data() + text.size()};
    const auto [ptr, ec]{std::from_chars(begin, end, value)};
    if (ec != std::errc{} || ptr != end) {
        throw std::invalid_argument{"Invalid floating-point value for --" + std::string{option_name} + ": '" +
                                    std::string{text} + "'"};
    }
    return value;
}

void validate_config(const Config& config) {
    if (config.num_particles < 1U) {
        throw std::invalid_argument{"--num_particles must be >= 1"};
    }
    if (!std::isfinite(config.box_length) || config.box_length <= 0.0) {
        throw std::invalid_argument{"--box_length must be finite and > 0"};
    }
    if (!std::isfinite(config.step_size) || config.step_size <= 0.0) {
        throw std::invalid_argument{"--step_size must be finite and > 0"};
    }
    if (config.measure_steps < 1U) {
        throw std::invalid_argument{"--measure_steps must be >= 1"};
    }
    if (config.block_size < 1U) {
        throw std::invalid_argument{"--block_size must be >= 1"};
    }
}

void print_usage(const char* program_name) {
    std::cout << "Usage:\n"
              << "  " << program_name
              << " --num_particles N --box_length L --warmup_steps W --measure_steps M --step_size S --seed R "
                 "--block_size B\n\n"
              << "Aliases:\n"
              << "  --num-particles, --box-length, --warmup-steps, --measure-steps, --step-size, --block-size\n";
}

[[nodiscard]] Config parse_args(int argc, char** argv) {
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

        std::string_view option_name{};
        std::string_view option_value{};
        const std::size_t equals_pos{token.find('=')};
        if (equals_pos == std::string_view::npos) {
            option_name = token;
            if (idx + 1 >= argc) {
                throw std::invalid_argument{"Missing value for option --" + std::string{option_name}};
            }
            ++idx;
            option_value = argv[idx];
        } else {
            option_name = token.substr(0, equals_pos);
            option_value = token.substr(equals_pos + 1);
        }

        const auto parsed_option{parse_option_name(option_name)};
        if (!parsed_option.has_value()) {
            throw std::invalid_argument{"Unknown option: --" + std::string{option_name}};
        }
        if (option_value.empty()) {
            throw std::invalid_argument{"Missing value for option --" +
                                        std::string{canonical_option_name(*parsed_option)}};
        }
        if (get_option_value(raw, *parsed_option).has_value()) {
            throw std::invalid_argument{"Duplicate option: --" + std::string{canonical_option_name(*parsed_option)}};
        }

        set_option_value(raw, *parsed_option, option_value);
    }

    constexpr std::array<OptionId, 7> all_options{
        OptionId::NUM_PARTICLES, OptionId::BOX_LENGTH, OptionId::WARMUP_STEPS, OptionId::MEASURE_STEPS,
        OptionId::STEP_SIZE,     OptionId::SEED,       OptionId::BLOCK_SIZE,
    };

    std::string missing{};
    for (const OptionId option : all_options) {
        if (!get_option_value(raw, option).has_value()) {
            if (!missing.empty()) {
                missing += ", ";
            }
            missing += "--";
            missing += canonical_option_name(option);
        }
    }
    if (!missing.empty()) {
        throw std::invalid_argument{"Missing required options: " + missing};
    }

    const Config config{
        .num_particles = parse_integer<std::size_t>(*raw.num_particles, "num_particles"),
        .box_length = parse_double(*raw.box_length, "box_length"),
        .warmup_steps = parse_integer<std::size_t>(*raw.warmup_steps, "warmup_steps"),
        .measure_steps = parse_integer<std::size_t>(*raw.measure_steps, "measure_steps"),
        .step_size = parse_double(*raw.step_size, "step_size"),
        .seed = parse_integer<std::uint64_t>(*raw.seed, "seed"),
        .block_size = parse_integer<std::size_t>(*raw.block_size, "block_size"),
    };
    validate_config(config);
    return config;
}

[[maybe_unused]] void print_config(const Config& config) {
    std::cout << "num_particles: " << config.num_particles << '\n'
              << "box_length: " << config.box_length << '\n'
              << "warmup_steps: " << config.warmup_steps << '\n'
              << "measure_steps: " << config.measure_steps << '\n'
              << "step_size: " << config.step_size << '\n'
              << "seed: " << config.seed << '\n'
              << "block_size: " << config.block_size << '\n';
}

} // Namespace
