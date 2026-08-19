#pragma once

#include <xpu/xpu.hpp>
#include <xpu/random.hpp>
#include "../energy_tracking/energy_tracking.hpp"
#include "../energy_tracking/energy_tracking_kernels.hpp"
#include "../jastrow_pade/jastrow_pade_kernels.hpp"
#include "../slater_plane_wave/slater_plane_wave.hpp"
#include "../slater_plane_wave/slater_plane_wave_kernels.hpp"
#include "simulation_types.hpp"

#if defined(XPU_CUDA)
  #include <cuda/std/limits>
  #include <cuda/std/numbers>
#endif

#include <cstddef>
#include <cstdint>

namespace stencil {
namespace simulation {

struct RandomProposal {
  std::size_t particle;
  xpu::array<real_t, idx(Axis::NUM)> displacement;
  real_t acceptance;
};

[[nodiscard]] DEVICE_ONLY
inline RandomProposal generate_random_proposal(
  xpu::random::generator& generator,
  std::size_t num_particles,
  real_t step_size
) noexcept {
  const auto particle{generator.uniform_index(num_particles)};

  xpu::array<real_t, idx(Axis::NUM)> displacement{};
  if (step_size > 0.0_r) {
    for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
      displacement[axis] = generator.uniform(-step_size, step_size);
    }
  }

  return {
    particle,
    displacement,
    generator.uniform<real_t>()
  };
}

} // namespace stencil::simulation
} // namespace stencil

namespace kernel {
namespace simulation {

#if defined(XPU_CUDA)

namespace {

__global__
void cudaSeedGenerator(
  xpu::random::generator* generator,
  std::uint64_t master_seed,
  std::uint64_t walker_id
) {
  generator->seed(master_seed, walker_id);
}

__global__
void cudaInitializePositions(
  xpu::random::generator* generator,
  real_t box_length,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos
) {
  if (threadIdx.x != 0u) { return; }

  auto local_generator{*generator};

  for (auto particle{0uz}; particle < pos.count(); ++particle) {
    pos[idx(Axis::X)][particle] = local_generator.uniform<real_t>() * box_length;
    pos[idx(Axis::Y)][particle] = local_generator.uniform<real_t>() * box_length;
    pos[idx(Axis::Z)][particle] = local_generator.uniform<real_t>() * box_length;
  }

  *generator = local_generator;
}

__global__
void cudaMetropolisStep(
  xpu::random::generator* generator,
  real_t step_size,
  real_t L,
  real_t jastrow_a,
  real_t jastrow_b,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  SlaterPlaneWave::View slater,
  EnergyTracker::View energy,
  ::simulation::StepResult* result
) {
  const auto current_thread{static_cast<std::size_t>(threadIdx.x)};
  const auto thread_stride{static_cast<std::size_t>(blockDim.x)};
  const auto num_orbitals{slater.num_orbitals};
  const xpu::soa_view<const real_t, idx(Axis::NUM)> const_pos{
    pos[idx(Axis::X)], pos.count()
  };

  __shared__ stencil::simulation::RandomProposal random;
  __shared__ ::simulation::StepResult step_result;
  __shared__ real_t ratio;
  __shared__ real_t delta_jastrow;
  __shared__ real_t real_energy_delta;
  __shared__ real_t reciprocal_sum;

  if (current_thread == 0uz) {
    auto local_generator{*generator};
    random = stencil::simulation::generate_random_proposal(
      local_generator,
      num_orbitals,
      step_size
    );
    *generator = local_generator;

    const auto particle{random.particle};
    const auto inv_L{1.0_r / L};

    step_result = {
      false,
      particle,
      {
        pos[idx(Axis::X)][particle],
        pos[idx(Axis::Y)][particle],
        pos[idx(Axis::Z)][particle]
      },
      {},
      0.0_r,
      0.0_r,
      0.0_r
    };

    for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
      auto new_position{step_result.old_pos[axis] + random.displacement[axis]};
      new_position -= L * xpu::floor(new_position * inv_L);

      pos[axis][particle] = new_position;
      step_result.new_pos[axis] = new_position;
    }

    ratio = 0.0_r;
    delta_jastrow = 0.0_r;
  }
  __syncthreads();

  const auto particle{random.particle};
  const auto trig_row_offset{particle * slater.trig_row_stride};

  for (auto i{current_thread}; i < slater.num_unique_k; i += thread_stride) {
    const auto trig_cache_index{trig_row_offset + i};
    slater.sin_saved[i] = slater.sin_cache[trig_cache_index];
    slater.cos_saved[i] = slater.cos_cache[trig_cache_index];

    stencil::slater::update_trig_cache(
      i, trig_row_offset, particle,
      pos, slater.k_vector,
      slater.sin_cache, slater.cos_cache
    );
  }
  __syncthreads();

  for (auto i{current_thread}; i < num_orbitals; i += thread_stride) {
    stencil::slater::build_row(
      i, particle, slater.trig_row_stride,
      slater.sin_cache, slater.cos_cache,
      slater.orbital_k_index, slater.orbital_type,
      slater.new_row
    );
  }
  __syncthreads();

  const auto half_L{0.5_r * L};
  for (auto i{current_thread}; i < num_orbitals; i += thread_stride) {
    stencil::slater::determinant_ratio(
      i, particle, slater.matrix_row_stride,
      slater.new_row, slater.inv_determinant,
      &ratio
    );
    stencil::jpade::delta_value(
      i, particle,
      L, half_L,
      jastrow_a, jastrow_b,
      step_result.old_pos, step_result.new_pos, pos,
      &delta_jastrow
    );
  }
  __syncthreads();

  if (current_thread == 0uz) {
    const auto log_psi_delta{xpu::log(xpu::abs(ratio)) + delta_jastrow};
    const auto log_ratio_sq{2.0_r * log_psi_delta};
    const auto uniform{xpu::max(random.acceptance, xstd::numeric_limits<real_t>::min())};
    const auto log_uniform{xpu::log(uniform)};
    const auto min_term{xpu::min(0.0_r, log_ratio_sq)};

    step_result.accepted = log_uniform < min_term;
    step_result.log_psi_delta = log_psi_delta;
  }
  __syncthreads();

  if (!step_result.accepted) {
    for (auto i{current_thread}; i < slater.num_unique_k; i += thread_stride) {
      const auto trig_cache_index{trig_row_offset + i};
      slater.sin_cache[trig_cache_index] = slater.sin_saved[i];
      slater.cos_cache[trig_cache_index] = slater.cos_saved[i];
    }

    if (current_thread == 0uz) {
      for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
        pos[axis][particle] = step_result.old_pos[axis];
      }
      *result = step_result;
    }
    return;
  }

  const auto determinant_row_offset{particle * slater.matrix_row_stride};
  for (auto i{current_thread}; i < num_orbitals; i += thread_stride) {
    const auto determinant_index{determinant_row_offset + i};
    slater.inv_d_col[i] = slater.inv_determinant[determinant_index];
    slater.solution[i] = 0.0_r;
  }
  __syncthreads();

  for (auto j{current_thread}; j < num_orbitals; j += thread_stride) {
    if (j == particle) { continue; }

    auto solution{0.0_r};
    const auto inverse_row_offset{j * slater.matrix_row_stride};
    for (auto i{0uz}; i < num_orbitals; ++i) {
      const auto inverse_index{inverse_row_offset + i};
      solution += slater.new_row[i] * slater.inv_determinant[inverse_index];
    }
    slater.solution[j] = solution;
  }
  __syncthreads();

  const auto inv_ratio{1.0_r / ratio};
  const auto matrix_elements{num_orbitals * num_orbitals};
  for (auto element{current_thread}; element < matrix_elements; element += thread_stride) {
    const auto j{element / num_orbitals};
    const auto i{element - j * num_orbitals};

    stencil::slater::k_update_inverse(
      i, j,
      particle, slater.matrix_row_stride, inv_ratio,
      slater.inv_d_col, slater.solution,
      slater.inv_determinant
    );
  }

  for (auto i{current_thread}; i < num_orbitals; i += thread_stride) {
    const auto determinant_index{determinant_row_offset + i};
    slater.determinant[determinant_index] = slater.new_row[i];
  }
  __syncthreads();

  for (auto i{current_thread}; i < energy.num_g_vectors; i += thread_stride) {
    stencil::energy::update_structure_factors(
      i,
      step_result.old_pos, step_result.new_pos,
      energy.g_vector,
      energy.sum_real, energy.sum_imag
    );
  }
  __syncthreads();

  if (current_thread == 0uz) {
    real_energy_delta = 0.0_r;
    reciprocal_sum = 0.0_r;
  }
  __syncthreads();

  for (auto i{current_thread}; i < energy.num_g_vectors; i += thread_stride) {
    stencil::energy::initialize_reciprocal_energy(
      i,
      energy.g_weights,
      energy.sum_real, energy.sum_imag,
      &reciprocal_sum
    );
  }

  for (auto i{current_thread}; i < num_orbitals; i += thread_stride) {
    stencil::energy::update_real_energy(
      i, particle,
      L, half_L, energy.ewald_alpha,
      step_result.old_pos, step_result.new_pos, const_pos,
      &real_energy_delta
    );
  }
  __syncthreads();

  if (current_thread == 0uz) {
    const auto reciprocal_prefactor{1.0_r / (2.0_r * xstd::numbers::pi_v<real_t> * L * L * L)};
    step_result.real_energy_delta = real_energy_delta;
    step_result.reciprocal_energy = reciprocal_prefactor * reciprocal_sum;
    *result = step_result;
  }
}

} // namespace

inline void seed_generator(
  xpu::random::generator* generator,
  std::uint64_t master_seed,
  std::uint64_t walker_id
) {
  dim3 seedGeneratorThreads{1u};
  dim3 seedGeneratorBlocks{1u};

  cudaSeedGenerator<<<
    seedGeneratorBlocks, seedGeneratorThreads
  >>>(
    generator,
    master_seed,
    walker_id
  );
  xpu::cu_check(cudaGetLastError());
}

inline void initialize_positions(
  xpu::random::generator* generator,
  real_t box_length,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos
) {
  dim3 initializePositionsThreads{1u};
  dim3 initializePositionsBlocks{1u};

  cudaInitializePositions<<<
    initializePositionsBlocks, initializePositionsThreads
  >>>(
    generator,
    box_length,
    pos
  );
  xpu::cu_check(cudaGetLastError());
}

inline ::simulation::StepResult metropolis_step(
  xpu::random::generator* generator,
  real_t step_size,
  real_t L,
  real_t jastrow_a,
  real_t jastrow_b,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  SlaterPlaneWave::View slater,
  EnergyTracker::View energy,
  ::simulation::StepResult* result_scratch
) {
  dim3 metropolisStepThreads{256u};
  dim3 metropolisStepBlocks{1u};

  cudaMetropolisStep<<<
    metropolisStepBlocks, metropolisStepThreads
  >>>(
    generator,
    step_size,
    L,
    jastrow_a,
    jastrow_b,
    pos,
    slater,
    energy,
    result_scratch
  );
  xpu::cu_check(cudaGetLastError());

  ::simulation::StepResult result{};
  xpu::copy_n(&result, result_scratch, 1uz);
  return result;
}

#endif

} // namespace kernel::simulation
} // namespace kernel
