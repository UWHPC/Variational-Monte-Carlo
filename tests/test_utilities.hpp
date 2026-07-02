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

#ifdef FP_64
constexpr real_t DEFAULT_TOLERANCE{1e-12_r};
#else
constexpr real_t DEFAULT_TOLERANCE{1e-5_r};
#endif

inline void require_near(real_t actual, real_t expected, real_t tolerance = DEFAULT_TOLERANCE) {
  REQUIRE(std::abs(actual - expected) <= tolerance);
}

inline void copy_positions(const Particles& source, Particles& dest) {
  const std::size_t N{source.size()};
  std::copy_n(source.pos().x_, N, dest.pos().x_);
  std::copy_n(source.pos().y_, N, dest.pos().y_);
  std::copy_n(source.pos().z_, N, dest.pos().z_);
}

inline void copy_derivatives(const Particles& source, Particles& dest) {
  const std::size_t N{source.size()};
  std::copy_n(source.grad_log_psi().x_, N, dest.grad_log_psi().x_);
  std::copy_n(source.grad_log_psi().y_, N, dest.grad_log_psi().y_);
  std::copy_n(source.grad_log_psi().z_, N, dest.grad_log_psi().z_);
  std::copy_n(source.lap_log_psi(), N, dest.lap_log_psi());
}

inline Particles copy_particle_positions(const Particles& source) {
  Particles copy{source.size()};
  copy_positions(source, copy);
  return copy;
}

inline std::size_t matrix_index(std::size_t row, std::size_t col, std::size_t n) {
  return row * n + col;
}

inline real_t determinant_3x3(const real_t* matrix, std::size_t stride) {
  return
    matrix[0 * stride + 0] * (
      matrix[1 * stride + 1] * matrix[2 * stride + 2] -
      matrix[1 * stride + 2] * matrix[2 * stride + 1]
    ) -
    matrix[0 * stride + 1] * (
      matrix[1 * stride + 0] * matrix[2 * stride + 2] -
      matrix[1 * stride + 2] * matrix[2 * stride + 0]
    ) +
    matrix[0 * stride + 2] * (
      matrix[1 * stride + 0] * matrix[2 * stride + 1] -
      matrix[1 * stride + 1] * matrix[2 * stride + 0]
    );
}

inline real_t slater_identity_residual(const SlaterPlaneWave& slater) {
  const std::size_t N{slater.num_orbitals()};
  const std::size_t S{slater.matrix_row_stride()};
  real_t max_residual{};

  for (std::size_t row = 0; row < N; ++row) {
    for (std::size_t col = 0; col < N; ++col) {
      real_t value{};

      for (std::size_t k = 0; k < N; ++k) {
        value += slater.determinant()[row * S + k] * slater.inv_determinant()[col * S + k];
      }

      const real_t expected{row == col ? 1.0_r : 0.0_r};
      const real_t residual{std::abs(value - expected)};
      max_residual = std::max(max_residual, residual);
    }
  }

  return max_residual;
}

inline real_t wrap_coordinate(real_t value, real_t box_length) {
  return value - box_length * std::floor(value / box_length);
}

inline real_t minimum_image(real_t dx, real_t box_length) {
  const real_t HALF_LENGTH{0.5_r * box_length};

  if (dx <= -HALF_LENGTH) {
    dx += box_length;
  } else if (dx > HALF_LENGTH) {
    dx -= box_length;
  }

  return dx;
}

inline real_t exact_kinetic_energy(const SlaterPlaneWave& slater) {
  const std::size_t N{slater.num_orbitals()};
  const real_t* k_x{slater.k_vector().x_};
  const real_t* k_y{slater.k_vector().y_};
  const real_t* k_z{slater.k_vector().z_};
  const auto& K_INDEX{slater.orbital_k_index()};

  real_t T_exact{};
  for (std::size_t j = 0; j < N; ++j) {
    const std::size_t IDX{K_INDEX[j]};
    T_exact += 0.5_r * (k_x[IDX] * k_x[IDX] + k_y[IDX] * k_y[IDX] + k_z[IDX] * k_z[IDX]);
  }
  return T_exact;
}

inline real_t local_kinetic_energy(const Particles& particles) {
  const std::size_t N{particles.size()};
  real_t T_local{};
  for (std::size_t i = 0; i < N; ++i) {
    const real_t GX{particles.grad_log_psi().x_[i]};
    const real_t GY{particles.grad_log_psi().y_[i]};
    const real_t GZ{particles.grad_log_psi().z_[i]};
    const real_t LAP{particles.lap_log_psi()[i]};
    T_local += -0.5_r * (LAP + GX * GX + GY * GY + GZ * GZ);
  }
  return T_local;
}

inline real_t box_length_from_rs(real_t r_s, std::size_t N) {
  return std::cbrt(4.0_r * std::numbers::pi_v<real_t> * static_cast<real_t>(N) / 3.0_r) * r_s;
}

inline void set_stable_closed_shell_positions(Particles& particles) {
  const std::size_t N{particles.size()};

  for (std::size_t i = 0; i < N; ++i) {
    particles.pos().x_[i] = 1.0_r + static_cast<real_t>(i) * 1.1_r;
    particles.pos().y_[i] = 0.5_r + static_cast<real_t>(i) * 0.7_r;
    particles.pos().z_[i] = 0.3_r + static_cast<real_t>(i) * 1.3_r;
  }
}

template <typename Fn> std::string capture_stdout(Fn&& fn) {
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
inline Config make_config(
  std::size_t num_particles,
  real_t box_length,
  std::size_t warmup_steps,
  std::size_t measure_steps,
  real_t step_size,
  uint64_t master_seed,
  std::size_t block_size,
  std::size_t num_threads = 1U,
  bool is_master_thread = false
) {
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
  real_t mean_energy;
  real_t standard_error;
  real_t acceptance_rate;
};