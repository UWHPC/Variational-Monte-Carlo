#pragma once

#include "particles/particles.hpp"
#include "slater_plane_wave/slater_plane_wave.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace validation {

template <std::size_t particle_count>
void set_positions(
  Particles& particles,
  const std::array<std::array<fp_t, idx(Axis::NUM)>, particle_count>& values
) {
  auto positions{particles.pos()};
  std::array<fp_t, particle_count> axis_values{};

  for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
    for (auto particle{0uz}; particle < particle_count; ++particle) {
      axis_values[particle] = values[particle][axis];
    }
    xpu::copy_n(positions[axis], axis_values.data(), axis_values.size());
  }
}

inline void zero_derivatives(Particles& particles) {
  auto derivatives{particles.derivatives()};
  for (auto component{idx(Derivatives::GRAD_X)};
       component < idx(Derivatives::NUM);
       ++component) {
    xpu::zero_n(derivatives[component], particles.count());
  }
}

[[nodiscard]] inline fp_t local_kinetic_energy(const Particles& particles) {
  std::array<std::vector<fp_t>, idx(Derivatives::NUM)> derivatives{};
  for (auto component{idx(Derivatives::GRAD_X)};
       component < idx(Derivatives::NUM);
       ++component) {
    derivatives[component].resize(particles.count());
    xpu::copy_n(
      derivatives[component].data(),
      particles.derivatives()[component],
      particles.count()
    );
  }

  fp_t energy{};
  for (auto particle{0uz}; particle < particles.count(); ++particle) {
    const auto gradient_squared{
      derivatives[idx(Derivatives::GRAD_X)][particle] *
        derivatives[idx(Derivatives::GRAD_X)][particle] +
      derivatives[idx(Derivatives::GRAD_Y)][particle] *
        derivatives[idx(Derivatives::GRAD_Y)][particle] +
      derivatives[idx(Derivatives::GRAD_Z)][particle] *
        derivatives[idx(Derivatives::GRAD_Z)][particle]
    };
    energy -= 0.5_fp * (
      derivatives[idx(Derivatives::LAP)][particle] + gradient_squared
    );
  }
  return energy;
}

[[nodiscard]] inline fp_t exact_free_gas_kinetic_energy(
  const SlaterPlaneWave& slater
) {
  const auto orbital_count{slater.num_orbitals()};
  const auto wave_vector_count{slater.num_unique_k()};
  std::vector<std::size_t> orbital_wave_vector(orbital_count);
  std::array<std::vector<fp_t>, idx(Axis::NUM)> wave_vectors{};

  xpu::copy_n(
    orbital_wave_vector.data(),
    slater.orbital_k_index(),
    orbital_count
  );
  for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
    wave_vectors[axis].resize(wave_vector_count);
    xpu::copy_n(
      wave_vectors[axis].data(),
      slater.k_vector()[axis],
      wave_vector_count
    );
  }

  fp_t energy{};
  for (auto orbital{0uz}; orbital < orbital_count; ++orbital) {
    const auto wave_vector{orbital_wave_vector[orbital]};
    for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
      energy += 0.5_fp * wave_vectors[axis][wave_vector] *
        wave_vectors[axis][wave_vector];
    }
  }
  return energy;
}

[[nodiscard]] inline fp_t wrap_coordinate(fp_t value, fp_t box_length) {
  return value - box_length * std::floor(value / box_length);
}

[[nodiscard]] inline fp_t determinant_3x3(
  const std::vector<fp_t>& matrix,
  std::size_t stride
) {
  return
    matrix[0uz * stride + 0uz] * (
      matrix[1uz * stride + 1uz] * matrix[2uz * stride + 2uz] -
      matrix[1uz * stride + 2uz] * matrix[2uz * stride + 1uz]
    ) -
    matrix[0uz * stride + 1uz] * (
      matrix[1uz * stride + 0uz] * matrix[2uz * stride + 2uz] -
      matrix[1uz * stride + 2uz] * matrix[2uz * stride + 0uz]
    ) +
    matrix[0uz * stride + 2uz] * (
      matrix[1uz * stride + 0uz] * matrix[2uz * stride + 1uz] -
      matrix[1uz * stride + 1uz] * matrix[2uz * stride + 0uz]
    );
}

[[nodiscard]] inline std::vector<fp_t> determinant_matrix(
  const SlaterPlaneWave& slater
) {
  std::vector<fp_t> matrix(slater.matrix_size());
  xpu::copy_n(matrix.data(), slater.determinant(), matrix.size());
  return matrix;
}

} // namespace validation
