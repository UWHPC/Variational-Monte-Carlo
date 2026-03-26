#pragma once

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>

class Config {
public:
    std::size_t num_threads{1};
    std::size_t num_particles{7U};
    std::size_t warmup_sweeps{100U};
    std::size_t measure_sweeps{100U};
    double box_length{10.0};
    std::size_t block_size{100U};
    uint64_t master_seed{42U};
    bool is_master_thread{};

    std::size_t warmup_steps;
    std::size_t measure_steps;
    double step_size;

    Config() { compute_derived(); }

    [[nodiscard]] static Config from_file(std::string const &path) {
        Config cfg{};
        std::ifstream file{path};
        if (!file.is_open()) {
            std::cout << "[Config] File not found: " << path
                      << "; using defaults.\n";
            return cfg;
        }
        auto map{parse(file)};
        read(map, "Num_Threads", cfg.num_threads);
        read(map, "Num_Particles", cfg.num_particles);
        read(map, "Warmup_Sweeps", cfg.warmup_sweeps);
        read(map, "Measure_Sweeps", cfg.measure_sweeps);
        read(map, "Box_Length", cfg.box_length);
        read(map, "Block_Size", cfg.block_size);
        read(map, "Master_Seed", cfg.master_seed);
        read(map, "Is_Master_Thread", cfg.is_master_thread);
        if (!map.empty()) {
            std::cout << "[Config] Warning! unknown keys ignored:";
            for (auto const &[key, val] : map) {
                std::cout << " " << key;
            }
            std::cout << "\n";
        }
        cfg.compute_derived();
        return cfg;
    }

private:
    void compute_derived() {
        if (num_particles < 1U) {
            throw std::invalid_argument("[Config] Num_Particles must be >= 1");
        }
        if (block_size < 1U) {
            throw std::invalid_argument("[Config] Block_Size must be >= 1");
        }
        if (box_length <= 0.0 || !std::isfinite(box_length)) {
            throw std::invalid_argument("[Config] Box_Length must be finite and > 0");
        }
        if (measure_sweeps < 1U) {
            throw std::invalid_argument("[Config] Measure_Sweeps must be >= 1");
        }
        this->warmup_steps = num_particles * warmup_sweeps;
        this->measure_steps = num_particles * measure_sweeps;
        this->step_size = box_length / 10.0;
    }

    using Map = std::unordered_map<std::string, std::string>;

    [[nodiscard]] static Map parse(std::ifstream &file) {
        Map map{};
        std::string line{};
        while (std::getline(file, line)) {
            if (auto pos{line.find('#')}; pos != std::string::npos) {
                line.erase(pos);
            }
            if (line.find_first_not_of(" \t\r") == std::string::npos) {
                continue;
            }
            auto eq{line.find('=')};
            if (eq == std::string::npos) {
                continue;
            }
            std::string key{trim(line.substr(0, eq))};
            std::string val{trim(line.substr(eq + 1))};
            if (!key.empty() && !val.empty()) {
                map[key] = val;
            }
        }
        return map;
    }

    [[nodiscard]] static std::string trim(std::string s) {
        auto const start{s.find_first_not_of(" \t\r")};
        if (start == std::string::npos) return {};
        auto const end{s.find_last_not_of(" \t\r")};
        return s.substr(start, end - start + 1);
    }

    static void read(Map &map, std::string const &key, double &out) {
        if (auto it{map.find(key)}; it != map.end()) {
            out = std::stod(it->second);
            map.erase(it);
        }
    }

    template <typename T>
        requires (std::is_integral_v<T> && !std::is_same_v<T, bool>)
    static void read(Map &map, std::string const &key, T &out) {
        if (auto it{map.find(key)}; it != map.end()) {
            if constexpr (std::is_signed_v<T>) {
                out = static_cast<T>(std::stoll(it->second));
            } else {
                out = static_cast<T>(std::stoull(it->second));
            }
            map.erase(it);
        }
    }

    static void read(Map &map, std::string const &key, bool &out) {
        if (auto it{map.find(key)}; it != map.end()) {
            std::string val{it->second};
            for (auto &ch : val) ch = static_cast<char>(std::tolower(ch));
            out = (val == "true" || val == "1");
            map.erase(it);
        }
    }
};