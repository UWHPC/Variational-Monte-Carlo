#pragma once

#include <xpu/xpu.hpp>
#include <xpu/random.hpp>
#include "../energy_tracking/energy_tracking.hpp"
#include "../energy_tracking/energy_tracking_kernels.hpp"
#include "../jastrow_pade/jastrow_pade_kernels.hpp"
#include "../slater_plane_wave/slater_plane_wave.hpp"
#include "../slater_plane_wave/slater_plane_wave_kernels.hpp"
#include "../wavefunction/wavefunction_kernels.hpp"
#include "../utilities/execution.hpp"
#include "simulation.hpp"

#if defined(XPU_CUDA)
  #include <cuda/std/limits>
  #include <cuda/std/numbers>
#else
  #include <limits>
  #include <numbers>
#endif

#include <cstddef>
#include <cstdint>

namespace stencil {
namespace simulation {

CUDA_CALLABLE
inline void initialize_walker_positions(
  const Simulation::View simulation,
  fp_t box_length
) noexcept {
  auto& generator{*simulation.generator};
  auto positions{simulation.particles.pos};

  for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
    for (auto particle{0uz}; particle < positions.count(); ++particle) {
      positions[axis][particle] = generator.uniform<fp_t>() * box_length;
    }
  }
}

[[nodiscard]] CUDA_CALLABLE
inline Simulation::RandomProposal generate_random_proposal(
  xpu::random::generator& generator,
  std::size_t num_particles,
  fp_t step_size
) noexcept {
  const auto particle{generator.uniform_index(num_particles)};

  xpu::array<fp_t, idx(Axis::NUM)> displacement{};
  if (step_size > 0.0_fp) {
    for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
      displacement[axis] = generator.uniform(-step_size, step_size);
    }
  }

  return {
    particle,
    displacement,
    generator.uniform<fp_t>()
  };
}

CUDA_CALLABLE
inline void propose_move(
  Simulation::View simulation,
  fp_t step_size,
  Simulation::MetropolisScratch& scratch
) noexcept {
  if (execution::thread() != 0uz) { return; }

  auto& generator{*simulation.generator};
  const auto L{simulation.wave_function.jastrow.box_length};
  auto pos{simulation.particles.pos};
  auto local_generator{generator};
  scratch.proposal = generate_random_proposal(
    local_generator, pos.count(), step_size
  );
  generator = local_generator;

  const auto particle{scratch.proposal.particle};
  const auto inv_L{1.0_fp / L};

  scratch.result = {
    false,
    particle,
    {
      pos[idx(Axis::X)][particle],
      pos[idx(Axis::Y)][particle],
      pos[idx(Axis::Z)][particle]
    },
    {},
    0.0_fp,
    0.0_fp,
    0.0_fp
  };

  for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
    auto new_position{
      scratch.result.old_pos[axis] + scratch.proposal.displacement[axis]
    };
    new_position -= L * xpu::floor(new_position * inv_L);

    pos[axis][particle] = new_position;
    scratch.result.new_pos[axis] = new_position;
  }

  scratch.slater_ratio = 0.0_fp;
  scratch.jastrow_delta = 0.0_fp;
}

CUDA_CALLABLE
inline void update_trig_cache(
  Simulation::View simulation,
  const Simulation::MetropolisScratch& scratch
) noexcept {
  const auto slater{simulation.wave_function.slater};
  const auto particle{scratch.proposal.particle};
  const auto trig_row_offset{particle * slater.trig_row_stride};

  for (auto i{execution::thread()}; i < slater.num_unique_k; i += execution::stride()) {
    const auto trig_cache_index{trig_row_offset + i};
    slater.sin_saved[i] = slater.sin_cache[trig_cache_index];
    slater.cos_saved[i] = slater.cos_cache[trig_cache_index];

    stencil::slater::update_trig_cache(
      i,
      particle,
      slater,
      simulation.particles
    );
  }
}

CUDA_CALLABLE
inline void build_slater_row(
  Simulation::View simulation,
  const Simulation::MetropolisScratch& scratch
) noexcept {
  const auto slater{simulation.wave_function.slater};
  const auto particle{scratch.proposal.particle};

  for (auto i{execution::thread()}; i < slater.num_orbitals; i += execution::stride()) {
    stencil::slater::build_row(
      i,
      particle,
      slater
    );
  }
}

CUDA_CALLABLE
inline void calculate_probability_ratio(
  Simulation::View simulation,
  Simulation::MetropolisScratch& scratch
) noexcept {
  const auto jastrow{simulation.wave_function.jastrow};
  const auto slater{simulation.wave_function.slater};
  const auto particle{scratch.proposal.particle};

  for (auto i{execution::thread()}; i < slater.num_orbitals; i += execution::stride()) {
    stencil::slater::determinant_ratio(
      i,
      particle,
      slater,
      &scratch.slater_ratio
    );
    stencil::jpade::delta_value(
      i,
      particle,
      jastrow,
      simulation.particles,
      scratch.result.old_pos,
      scratch.result.new_pos,
      &scratch.jastrow_delta
    );
  }
}

CUDA_CALLABLE
inline void decide_move(
  Simulation::MetropolisScratch& scratch
) noexcept {
  if (execution::thread() != 0uz) { return; }

  const auto log_psi_delta{
    xpu::log(xpu::abs(scratch.slater_ratio)) + scratch.jastrow_delta
  };
  const auto log_ratio_sq{2.0_fp * log_psi_delta};
  const auto uniform{
    xpu::max(
      scratch.proposal.acceptance,
      xstd::numeric_limits<fp_t>::min()
    )
  };

  scratch.result.accepted = xpu::log(uniform) < xpu::min(0.0_fp, log_ratio_sq);
  scratch.result.log_psi_delta = log_psi_delta;
}

CUDA_CALLABLE
inline void reject_move(
  Simulation::View simulation,
  const Simulation::MetropolisScratch& scratch
) noexcept {
  auto pos{simulation.particles.pos};
  const auto slater{simulation.wave_function.slater};
  const auto particle{scratch.proposal.particle};
  const auto trig_row_offset{particle * slater.trig_row_stride};

  for (auto i{execution::thread()}; i < slater.num_unique_k; i += execution::stride()) {
    const auto trig_cache_index{trig_row_offset + i};
    slater.sin_cache[trig_cache_index] = slater.sin_saved[i];
    slater.cos_cache[trig_cache_index] = slater.cos_saved[i];
  }

  if (execution::thread() == 0uz) {
    for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
      pos[axis][particle] = scratch.result.old_pos[axis];
    }
  }
}

CUDA_CALLABLE
inline void prepare_inverse_update(
  Simulation::View simulation,
  const Simulation::MetropolisScratch& scratch
) noexcept {
  const auto slater{simulation.wave_function.slater};
  const auto particle{scratch.proposal.particle};
  const auto determinant_row_offset{particle * slater.matrix_row_stride};

  for (auto i{execution::thread()}; i < slater.num_orbitals; i += execution::stride()) {
    slater.inv_d_col[i] = slater.inv_determinant[determinant_row_offset + i];
    slater.solution[i] = 0.0_fp;
  }
}

CUDA_CALLABLE
inline void calculate_inverse_solution(
  Simulation::View simulation,
  const Simulation::MetropolisScratch& scratch
) noexcept {
  const auto slater{simulation.wave_function.slater};
  const auto particle{scratch.proposal.particle};

  for (auto j{execution::thread()}; j < slater.num_orbitals; j += execution::stride()) {
    if (j == particle) { continue; }

    auto solution{0.0_fp};
    const auto inverse_row_offset{j * slater.matrix_row_stride};
    for (auto i{0uz}; i < slater.num_orbitals; ++i) {
      solution += slater.new_row[i] * slater.inv_determinant[inverse_row_offset + i];
    }
    slater.solution[j] = solution;
  }
}

CUDA_CALLABLE
inline void commit_slater_move(
  Simulation::View simulation,
  const Simulation::MetropolisScratch& scratch
) noexcept {
  const auto slater{simulation.wave_function.slater};
  const auto particle{scratch.proposal.particle};
  const auto inv_ratio{1.0_fp / scratch.slater_ratio};
  const auto matrix_elements{slater.num_orbitals * slater.num_orbitals};

  for (auto element{execution::thread()}; element < matrix_elements; element += execution::stride()) {
    const auto j{element / slater.num_orbitals};
    const auto i{element - j * slater.num_orbitals};

    stencil::slater::k_update_inverse(
      i,
      j,
      particle,
      inv_ratio,
      slater
    );
  }

  const auto determinant_row_offset{particle * slater.matrix_row_stride};
  for (auto i{execution::thread()}; i < slater.num_orbitals; i += execution::stride()) {
    slater.determinant[determinant_row_offset + i] = slater.new_row[i];
  }
}

CUDA_CALLABLE
inline void update_structure_factors(
  Simulation::View simulation,
  const Simulation::MetropolisScratch& scratch
) noexcept {
  const auto energy{simulation.energy_tracker};
  for (auto i{execution::thread()}; i < energy.num_g_vectors; i += execution::stride()) {
    stencil::energy::update_structure_factors(
      i,
      energy,
      scratch.result.old_pos,
      scratch.result.new_pos
    );
  }
}

CUDA_CALLABLE
inline void reset_energy_reductions(
  Simulation::MetropolisScratch& scratch
) noexcept {
  if (execution::thread() != 0uz) { return; }
  scratch.real_energy_delta = 0.0_fp;
  scratch.reciprocal_sum = 0.0_fp;
}

CUDA_CALLABLE
inline void calculate_energy_deltas(
  Simulation::View simulation,
  Simulation::MetropolisScratch& scratch
) noexcept {
  const auto slater{simulation.wave_function.slater};
  const auto energy{simulation.energy_tracker};
  const auto particle{scratch.proposal.particle};

  for (auto i{execution::thread()}; i < energy.num_g_vectors; i += execution::stride()) {
    stencil::energy::initialize_reciprocal_energy(
      i,
      energy,
      &scratch.reciprocal_sum
    );
  }

  for (auto i{execution::thread()}; i < slater.num_orbitals; i += execution::stride()) {
    stencil::energy::update_real_energy(
      i,
      particle,
      energy,
      simulation.particles,
      scratch.result.old_pos,
      scratch.result.new_pos,
      &scratch.real_energy_delta
    );
  }
}

CUDA_CALLABLE
inline void finalize_energy(
  Simulation::View simulation,
  Simulation::MetropolisScratch& scratch
) noexcept {
  if (execution::thread() != 0uz) { return; }

  const auto reciprocal_prefactor{
    1.0_fp / (
      2.0_fp * xstd::numbers::pi_v<fp_t>
      * simulation.energy_tracker.box_length
      * simulation.energy_tracker.box_length
      * simulation.energy_tracker.box_length
    )
  };
  scratch.result.real_energy_delta = scratch.real_energy_delta;
  scratch.result.reciprocal_energy = reciprocal_prefactor * scratch.reciprocal_sum;
}

CUDA_CALLABLE
inline void commit_energy(
  Simulation::View simulation,
  const Simulation::MetropolisScratch& scratch
) noexcept {
  if (execution::thread() != 0uz) { return; }

  *simulation.energy_tracker.real_energy += scratch.result.real_energy_delta;
  *simulation.energy_tracker.reciprocal_energy = scratch.result.reciprocal_energy;
}

CUDA_CALLABLE
inline void metropolis_step(
  Simulation::View simulation,
  fp_t step_size,
  Simulation::MetropolisScratch& scratch
) noexcept {
  propose_move(simulation, step_size, scratch);
  execution::sync();

  update_trig_cache(simulation, scratch);
  execution::sync();

  build_slater_row(simulation, scratch);
  execution::sync();

  calculate_probability_ratio(simulation, scratch);
  execution::sync();

  decide_move(scratch);
  execution::sync();

  if (!scratch.result.accepted) {
    reject_move(simulation, scratch);
    execution::sync();
    return;
  }

  prepare_inverse_update(simulation, scratch);
  execution::sync();

  calculate_inverse_solution(simulation, scratch);
  execution::sync();

  commit_slater_move(simulation, scratch);
  execution::sync();

  update_structure_factors(simulation, scratch);
  execution::sync();

  reset_energy_reductions(scratch);
  execution::sync();

  calculate_energy_deltas(simulation, scratch);
  execution::sync();

  finalize_energy(simulation, scratch);
  execution::sync();

  commit_energy(simulation, scratch);
  execution::sync();
}

CUDA_CALLABLE
inline void measure_walker(Simulation::View simulation) noexcept {
  stencil::wavefunction::evaluate_derivatives(
    simulation.wave_function,
    simulation.particles
  );
  execution::sync();

  stencil::energy::evaluate_local_energy(
    simulation.energy_tracker,
    simulation.particles,
    simulation.local_energy
  );
  execution::sync();
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
void cudaInitializeWalkerPositions(
  const Simulation::View simulation,
  fp_t box_length
) {
  stencil::simulation::initialize_walker_positions(
    simulation,
    box_length
  );
}

__global__
void cudaInitializePositions(
  const Simulation::View* simulations,
  std::size_t walker_count,
  fp_t box_length
) {
  const auto [walker]{xpu::global_index<1>()};
  if (walker >= walker_count) { return; }

  stencil::simulation::initialize_walker_positions(
    simulations[walker],
    box_length
  );
}

__global__
void cudaMetropolisStep(
  Simulation::View simulation,
  fp_t step_size
) {
  __shared__ Simulation::MetropolisScratch scratch;

  stencil::simulation::metropolis_step(
    simulation,
    step_size,
    scratch
  );

  if (execution::thread() == 0uz) {
    *simulation.step_result = scratch.result;
  }
}

__global__
void cudaMetropolisSweep(
  Simulation::View* simulations,
  std::size_t walker_count,
  std::size_t proposals_per_walker,
  fp_t step_size,
  Simulation::SweepResult* result
) {
  const auto walker{scast<std::size_t>(blockIdx.x)};
  if (walker >= walker_count) { return; }

  __shared__ Simulation::MetropolisScratch scratch;
  __shared__ Simulation::SweepResult local_result;

  if (execution::thread() == 0uz) {
    local_result = {};
  }
  execution::sync();

  for (auto proposal{0uz}; proposal < proposals_per_walker; ++proposal) {
    stencil::simulation::metropolis_step(
      simulations[walker],
      step_size,
      scratch
    );

    if (execution::thread() == 0uz) {
      ++local_result.proposed;
      local_result.accepted += scast<std::size_t>(scratch.result.accepted);
    }
    execution::sync();
  }

  if (execution::thread() == 0uz) {
    static_assert(sizeof(std::size_t) == sizeof(unsigned long long));
    atomicAdd(
      reinterpret_cast<unsigned long long*>(&result->proposed),
      scast<unsigned long long>(local_result.proposed)
    );
    atomicAdd(
      reinterpret_cast<unsigned long long*>(&result->accepted),
      scast<unsigned long long>(local_result.accepted)
    );
  }
}

__global__
void cudaMeasureWalkers(
  Simulation::View* simulations,
  std::size_t walker_count
) {
  const auto walker{scast<std::size_t>(blockIdx.x)};
  if (walker >= walker_count) { return; }

  stencil::simulation::measure_walker(simulations[walker]);
}

} // namespace

#endif

inline void seed_generator(
  xpu::random::generator* generator,
  std::uint64_t master_seed,
  std::uint64_t walker_id
) {
#if defined(XPU_CUDA)
  constexpr dim3 seedGeneratorThreads{1u};
  constexpr dim3 seedGeneratorBlocks{1u};

  cudaSeedGenerator<<<
    seedGeneratorBlocks, seedGeneratorThreads
  >>>(
    generator,
    master_seed,
    walker_id
  );
  xpu::cu_check(cudaGetLastError());
#else
  generator->seed(master_seed, walker_id);
#endif
}

inline void initialize_positions(
  const Simulation::View simulation,
  fp_t box_length
) {
#if defined(XPU_CUDA)
  constexpr dim3 initializePositionsThreads{1u};
  constexpr dim3 initializePositionsBlocks{1u};

  cudaInitializeWalkerPositions<<<
    initializePositionsBlocks, initializePositionsThreads
  >>>(
    simulation,
    box_length
  );
  xpu::cu_check(cudaGetLastError());
#else
  stencil::simulation::initialize_walker_positions(
    simulation,
    box_length
  );
#endif
}

inline void initialize_positions(
  const Simulation::View* simulations,
  std::size_t walker_count,
  fp_t box_length
) {
#if defined(XPU_CUDA)
  constexpr dim3 initializePositionsThreads{256u};
  const dim3 initializePositionsBlocks{
    xpu::block_per_dim(walker_count, initializePositionsThreads.x)
  };

  cudaInitializePositions<<<
    initializePositionsBlocks, initializePositionsThreads
  >>>(
    simulations,
    walker_count,
    box_length
  );
  xpu::cu_check(cudaGetLastError());
#else
  for (auto walker{0uz}; walker < walker_count; ++walker) {
    stencil::simulation::initialize_walker_positions(
      simulations[walker],
      box_length
    );
  }
#endif
}

inline Simulation::StepResult metropolis_step(
  Simulation::View simulation,
  fp_t step_size
) {
#if defined(XPU_CUDA)
  dim3 metropolisStepThreads{256u};
  dim3 metropolisStepBlocks{1u};
  cudaMetropolisStep<<<
    metropolisStepBlocks, metropolisStepThreads
  >>>(
    simulation,
    step_size
  );
  xpu::cu_check(cudaGetLastError());

  Simulation::StepResult result{};
  xpu::copy_n(&result, simulation.step_result, 1uz);
  return result;
#else
  Simulation::MetropolisScratch scratch{};
  stencil::simulation::metropolis_step(
    simulation,
    step_size,
    scratch
  );
  *simulation.step_result = scratch.result;
  return scratch.result;
#endif
}

inline void metropolis_sweep(
  Simulation::View* simulations,
  std::size_t walker_count,
  std::size_t proposals_per_walker,
  fp_t step_size,
  Simulation::SweepResult* result
) {
  xpu::memset(result, 0, sizeof(Simulation::SweepResult));

#if defined(XPU_CUDA)
  dim3 metropolisSweepThreads{256u};
  dim3 metropolisSweepBlocks{scast<unsigned int>(walker_count)};
  cudaMetropolisSweep<<<
    metropolisSweepBlocks, metropolisSweepThreads
  >>>(
    simulations,
    walker_count,
    proposals_per_walker,
    step_size,
    result
  );
  xpu::cu_check(cudaGetLastError());
#else
  for (auto walker{0uz}; walker < walker_count; ++walker) {
    Simulation::SweepResult local_result{};

    for (auto proposal{0uz}; proposal < proposals_per_walker; ++proposal) {
      Simulation::MetropolisScratch scratch{};
      stencil::simulation::metropolis_step(
        simulations[walker],
        step_size,
        scratch
      );

      ++local_result.proposed;
      local_result.accepted += scast<std::size_t>(scratch.result.accepted);
    }

    result->proposed += local_result.proposed;
    result->accepted += local_result.accepted;
  }
#endif
}

inline void measure_walkers(
  Simulation::View* simulations,
  std::size_t walker_count
) {
#if defined(XPU_CUDA)
  dim3 measureWalkersThreads{256u};
  dim3 measureWalkersBlocks{scast<unsigned int>(walker_count)};
  cudaMeasureWalkers<<<
    measureWalkersBlocks, measureWalkersThreads
  >>>(
    simulations,
    walker_count
  );
  xpu::cu_check(cudaGetLastError());
#else
  for (auto walker{0uz}; walker < walker_count; ++walker) {
    stencil::simulation::measure_walker(simulations[walker]);
  }
#endif
}

} // namespace kernel::simulation
} // namespace kernel
