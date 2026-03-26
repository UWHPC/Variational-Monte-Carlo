#pragma once

#include <catch2/catch_test_macros.hpp>

#include "config/config.hpp"
#include "output_writer/output_writer.hpp"
#include "particles/particles.hpp"
#include "slater_plane_wave/slater_plane_wave.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <numbers>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

inline void require_near(double actual, double expected, double tolerance = 1e-12) {
    REQUIRE(std::abs(actual - expected) <= tolerance);
}

inline void copy_positions(const Particles& source, Particles& dest) {
    const std::size_t N{source.num_particles_get()};
    std::copy_n(source.pos_x_get(), N, dest.pos_x_get());
    std::copy_n(source.pos_y_get(), N, dest.pos_y_get());
    std::copy_n(source.pos_z_get(), N, dest.pos_z_get());
}

inline void copy_derivatives(const Particles& source, Particles& dest) {
    const std::size_t N{source.num_particles_get()};
    std::copy_n(source.grad_log_psi_x_get(), N, dest.grad_log_psi_x_get());
    std::copy_n(source.grad_log_psi_y_get(), N, dest.grad_log_psi_y_get());
    std::copy_n(source.grad_log_psi_z_get(), N, dest.grad_log_psi_z_get());
    std::copy_n(source.lap_log_psi_get(), N, dest.lap_log_psi_get());
}

inline Particles copy_particle_positions(const Particles& source) {
    Particles copy{source.num_particles_get()};
    copy_positions(source, copy);
    return copy;
}

inline std::size_t matrix_index(std::size_t row, std::size_t col, std::size_t n) {
    return row * n + col;
}

inline double determinant_3x3(const double* matrix) {
    return matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7]) -
           matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6]) +
           matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6]);
}

inline double slater_identity_residual(const SlaterPlaneWave& slater) {
    const std::size_t N{slater.num_orbitals_get()};
    double max_residual{};

    for (std::size_t row = 0; row < N; ++row) {
        for (std::size_t col = 0; col < N; ++col) {
            double value{};

            for (std::size_t k = 0; k < N; ++k) {
                value += slater.determinant_get()[row * N + k] *
                         slater.inv_determinant_get()[col * N + k];
            }

            const double expected{row == col ? 1.0 : 0.0};
            const double residual{std::abs(value - expected)};
            max_residual = std::max(max_residual, residual);
        }
    }

    return max_residual;
}

inline double wrap_coordinate(double value, double box_length) {
    return value - box_length * std::floor(value / box_length);
}

inline double minimum_image(double dx, double box_length) {
    const double HALF_LENGTH{0.5 * box_length};

    if (dx <= -HALF_LENGTH) {
        dx += box_length;
    } else if (dx > HALF_LENGTH) {
        dx -= box_length;
    }

    return dx;
}

inline double exact_kinetic_energy(const SlaterPlaneWave& slater) {
    const std::size_t N{slater.num_orbitals_get()};
    const double* k_x{slater.k_vector_x_get()};
    const double* k_y{slater.k_vector_y_get()};
    const double* k_z{slater.k_vector_z_get()};
    const auto& K_INDEX{slater.orbital_k_index_get()};

    double T_exact{};
    for (std::size_t j = 0; j < N; ++j) {
        const std::size_t IDX{K_INDEX[j]};
        T_exact += 0.5 * (k_x[IDX] * k_x[IDX] + k_y[IDX] * k_y[IDX] + k_z[IDX] * k_z[IDX]);
    }
    return T_exact;
}

inline double local_kinetic_energy(const Particles& particles) {
    const std::size_t N{particles.num_particles_get()};
    double T_local{};
    for (std::size_t i = 0; i < N; ++i) {
        const double GX{particles.grad_log_psi_x_get()[i]};
        const double GY{particles.grad_log_psi_y_get()[i]};
        const double GZ{particles.grad_log_psi_z_get()[i]};
        const double LAP{particles.lap_log_psi_get()[i]};
        T_local += -0.5 * (LAP + GX * GX + GY * GY + GZ * GZ);
    }
    return T_local;
}

inline double box_length_from_rs(double r_s, std::size_t N) {
    return std::cbrt(4.0 * std::numbers::pi * static_cast<double>(N) / 3.0) * r_s;
}

inline void set_stable_closed_shell_positions(Particles& particles) {
    const std::size_t N{particles.num_particles_get()};

    for (std::size_t i = 0; i < N; ++i) {
        particles.pos_x_get()[i] = 1.0 + static_cast<double>(i) * 1.1;
        particles.pos_y_get()[i] = 0.5 + static_cast<double>(i) * 0.7;
        particles.pos_z_get()[i] = 0.3 + static_cast<double>(i) * 1.3;
    }
}

template <typename Fn>
std::string capture_stdout(Fn&& fn) {
    std::ostringstream output{};
    std::streambuf* const OLD_BUFFER{std::cout.rdbuf(output.rdbuf())};
    try {
        std::forward<Fn>(fn)();
    } catch (...) {
        std::cout.rdbuf(OLD_BUFFER);
        throw;
    }
    std::cout.rdbuf(OLD_BUFFER);
    return output.str();
}

/// Builds a Config with explicit control over derived fields.
/// warmup_steps and measure_steps are per-particle steps (not sweeps).
/// step_size overrides the default box_length/10 derivation.
inline Config make_config(std::size_t num_particles, double box_length,
                          std::size_t warmup_steps, std::size_t measure_steps,
                          double step_size, uint64_t master_seed,
                          std::size_t block_size, std::size_t num_threads = 1U,
                          bool is_master_thread = false) {
    Config cfg{};
    cfg.num_threads = num_threads;
    cfg.num_particles = num_particles;
    cfg.box_length = box_length;
    cfg.block_size = block_size;
    cfg.master_seed = master_seed;
    cfg.is_master_thread = is_master_thread;

    // Set derived fields directly (bypass compute_derived sweep logic):
    cfg.warmup_sweeps = (num_particles > 0U) ? (warmup_steps / num_particles) : 0U;
    cfg.measure_sweeps = (num_particles > 0U) ? (measure_steps / num_particles) : 0U;
    cfg.warmup_steps = warmup_steps;
    cfg.measure_steps = measure_steps;
    cfg.step_size = step_size;
    return cfg;
}

class RecordingOutputWriter final : public OutputWriter {
public:
    void write_init(const InitData& data) override {
        init = data;
        saw_init = true;
    }

    void write_frame(const FrameData& data) override { frames.push_back(data); }

    void write_done(const DoneData& data) override {
        done = data;
        saw_done = true;
    }

    bool saw_init{false};
    bool saw_done{false};
    std::optional<InitData> init{};
    std::optional<DoneData> done{};
    std::vector<FrameData> frames{};
};

struct SimResult {
    double mean_energy;
    double standard_error;
    double acceptance_rate;
};