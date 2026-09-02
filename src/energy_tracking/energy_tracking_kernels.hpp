#pragma once

#include <xpu/xpu.hpp>
#include "energy_tracking.hpp"
#include "../utilities/components.hpp"
#include "../utilities/macros.hpp"

namespace stencil {
namespace energy {

namespace {

CUDA_CALLABLE
inline fp_t evaluate_pair_energy(
  std::size_t other,
  fp_t L, fp_t half_L, fp_t alpha,
  const xpu::array<fp_t, idx(Axis::NUM)>& particle_pos,
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos
) noexcept {
  xpu::array<fp_t, idx(Axis::NUM)> displacement{};

  for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
    auto delta{particle_pos[axis] - pos[axis][other]};

    delta += L * (delta <= -half_L) - L * (delta > half_L);

    displacement[axis] = delta;
  }

  const auto distance{
    xpu::sqrt(
      displacement[idx(Axis::X)] * displacement[idx(Axis::X)] +
      displacement[idx(Axis::Y)] * displacement[idx(Axis::Y)] +
      displacement[idx(Axis::Z)] * displacement[idx(Axis::Z)]
    )
  };
  const auto inverse_distance{
    (distance < 1e-12_fp) ? 1.0_fp : 1.0_fp / distance
  };

  return xpu::erfc(alpha * distance) * inverse_distance;
}

} // namespace

CUDA_CALLABLE
inline void initialize_reciprocal_energy(
  std::size_t i,
  const fp_t* RESTRICT g_weights,
  const fp_t* RESTRICT sum_real,
  const fp_t* RESTRICT sum_imag,
  fp_t* RESTRICT reciprocal_sum
) noexcept {
  const auto contribution{
    g_weights[i] * (
      sum_real[i] * sum_real[i] +
      sum_imag[i] * sum_imag[i]
    )
  };

#if defined(__CUDA_ARCH__)
  atomicAdd(reciprocal_sum, contribution);
#else
  *reciprocal_sum += contribution;
#endif
}

CUDA_CALLABLE
inline void initialize_real_energy(
  std::size_t other,
  fp_t L, fp_t half_L, fp_t alpha,
  const xpu::array<fp_t, idx(Axis::NUM)>& particle_pos,
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos,
  fp_t* RESTRICT real_sum
) noexcept {
  const auto pair_energy{
    evaluate_pair_energy(
      other,
      L, half_L, alpha,
      particle_pos, pos
    )
  };

#if defined(__CUDA_ARCH__)
  atomicAdd(real_sum, pair_energy);
#else
  *real_sum += pair_energy;
#endif
}

CUDA_CALLABLE
inline void update_real_energy(
  std::size_t other, std::size_t moved,
  fp_t L, fp_t half_L, fp_t alpha,
  const xpu::array<fp_t, idx(Axis::NUM)>& old_pos,
  const xpu::array<fp_t, idx(Axis::NUM)>& new_pos,
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos,
  fp_t* RESTRICT delta
) noexcept {
  if (other == moved) { return; }

  const auto old_pair_energy{
    evaluate_pair_energy(
      other,
      L, half_L, alpha,
      old_pos, pos
    )
  };
  const auto new_pair_energy{
    evaluate_pair_energy(
      other,
      L, half_L, alpha,
      new_pos, pos
    )
  };

#if defined(__CUDA_ARCH__)
  atomicAdd(delta, new_pair_energy - old_pair_energy);
#else
  *delta += new_pair_energy - old_pair_energy;
#endif
}

CUDA_CALLABLE
inline void kinetic_energy(
  std::size_t i,
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> derivatives,
  fp_t* RESTRICT kinetic_sum
) noexcept {
  const auto gradient_x{derivatives[idx(Derivatives::GRAD_X)][i]};
  const auto gradient_y{derivatives[idx(Derivatives::GRAD_Y)][i]};
  const auto gradient_z{derivatives[idx(Derivatives::GRAD_Z)][i]};
  const auto gradient_squared{
    gradient_x * gradient_x +
    gradient_y * gradient_y +
    gradient_z * gradient_z
  };
  const auto contribution{
    derivatives[idx(Derivatives::LAP)][i] + gradient_squared
  };

#if defined(__CUDA_ARCH__)
  atomicAdd(kinetic_sum, contribution);
#else
  *kinetic_sum += contribution;
#endif
}

CUDA_CALLABLE
inline void initialize_structure_factors(
  std::size_t g, std::size_t particle,
  xpu::soa_view<fp_t, idx(Axis::NUM)> g_vector,
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos,
  fp_t* RESTRICT sum_real,
  fp_t* RESTRICT sum_imag
) noexcept {
  const auto g_dot_r{
    g_vector[idx(Axis::X)][g] * pos[idx(Axis::X)][particle] +
    g_vector[idx(Axis::Y)][g] * pos[idx(Axis::Y)][particle] +
    g_vector[idx(Axis::Z)][g] * pos[idx(Axis::Z)][particle]
  };

  auto sin_term{0.0_fp};
  auto cos_term{0.0_fp};
  xpu::sincos(g_dot_r, &sin_term, &cos_term);

#if defined(__CUDA_ARCH__)
  atomicAdd(sum_real, cos_term);
  atomicAdd(sum_imag, sin_term);
#else
  *sum_real += cos_term;
  *sum_imag += sin_term;
#endif
}

CUDA_CALLABLE
inline void update_structure_factors(
  std::size_t i,
  const xpu::array<fp_t, idx(Axis::NUM)>& old_pos,
  const xpu::array<fp_t, idx(Axis::NUM)>& new_pos,
  xpu::soa_view<fp_t, idx(Axis::NUM)> g_vector,
  fp_t* RESTRICT sum_real,
  fp_t* RESTRICT sum_imag
) noexcept {
  const auto old_dot{
    g_vector[idx(Axis::X)][i] * old_pos[idx(Axis::X)] +
    g_vector[idx(Axis::Y)][i] * old_pos[idx(Axis::Y)] +
    g_vector[idx(Axis::Z)][i] * old_pos[idx(Axis::Z)]
  };
  const auto new_dot{
    g_vector[idx(Axis::X)][i] * new_pos[idx(Axis::X)] +
    g_vector[idx(Axis::Y)][i] * new_pos[idx(Axis::Y)] +
    g_vector[idx(Axis::Z)][i] * new_pos[idx(Axis::Z)]
  };

  auto old_sin{0.0_fp};
  auto old_cos{0.0_fp};
  auto new_sin{0.0_fp};
  auto new_cos{0.0_fp};

  xpu::sincos(old_dot, &old_sin, &old_cos);
  xpu::sincos(new_dot, &new_sin, &new_cos);

  sum_real[i] += new_cos - old_cos;
  sum_imag[i] += new_sin - old_sin;
}

} // namespace stencil::energy
} // namespace stencil

namespace kernel {
namespace energy {

namespace {

#if defined(XPU_CUDA)
__global__
void cudaInitializeReciprocalEnergy(
  std::size_t num_g_vectors,
  const fp_t* RESTRICT g_weights,
  const fp_t* RESTRICT sum_real,
  const fp_t* RESTRICT sum_imag,
  fp_t* RESTRICT reciprocal_sum
) {
  const auto [i]{xpu::global_index<1>()};
  if (i >= num_g_vectors) { return; }

  stencil::energy::initialize_reciprocal_energy(
    i,
    g_weights, sum_real, sum_imag,
    reciprocal_sum
  );
}

__global__
void cudaInitializeRealEnergy(
  fp_t L, fp_t half_L, fp_t alpha,
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos,
  fp_t* RESTRICT real_sum
) {
  const auto [i, j]{xpu::global_index<2>()};
  if (i >= pos.count() || j >= pos.count()) { return; }
  if (i >= j) { return; }

  xpu::array<fp_t, idx(Axis::NUM)> particle_pos{
    pos[idx(Axis::X)][i],
    pos[idx(Axis::Y)][i],
    pos[idx(Axis::Z)][i]
  };

  stencil::energy::initialize_real_energy(
    j,
    L, half_L, alpha,
    particle_pos, pos,
    real_sum
  );
}

__global__
void cudaUpdateRealEnergy(
  std::size_t moved,
  fp_t L, fp_t half_L, fp_t alpha,
  xpu::array<fp_t, idx(Axis::NUM)> old_pos,
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos,
  fp_t* RESTRICT delta
) {
  const auto [i]{xpu::global_index<1>()};
  if (i >= pos.count()) { return; }

  xpu::array<fp_t, idx(Axis::NUM)> new_pos{
    pos[idx(Axis::X)][moved],
    pos[idx(Axis::Y)][moved],
    pos[idx(Axis::Z)][moved]
  };

  stencil::energy::update_real_energy(
    i, moved,
    L, half_L, alpha,
    old_pos, new_pos, pos,
    delta
  );
}

__global__
void cudaKineticEnergy(
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> derivatives,
  fp_t* RESTRICT kinetic_sum
) {
  const auto [i]{xpu::global_index<1>()};
  if (i >= derivatives.count()) { return; }

  stencil::energy::kinetic_energy(
    i,
    derivatives,
    kinetic_sum
  );
}

__global__
void cudaInitializeStructureFactors(
  xpu::soa_view<fp_t, idx(Axis::NUM)> g_vector,
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos,
  fp_t* RESTRICT sum_real,
  fp_t* RESTRICT sum_imag
) {
  const auto [i, j]{xpu::global_index<2>()};
  if (i >= g_vector.count() || j >= pos.count()) { return; }

  stencil::energy::initialize_structure_factors(
    i, j,
    g_vector, pos,
    &sum_real[i], &sum_imag[i]
  );
}

__global__
void cudaUpdateStructureFactors(
  xpu::array<fp_t, idx(Axis::NUM)> old_pos,
  xpu::array<fp_t, idx(Axis::NUM)> new_pos,
  xpu::soa_view<fp_t, idx(Axis::NUM)> g_vector,
  fp_t* RESTRICT sum_real,
  fp_t* RESTRICT sum_imag
) {
  const auto [i]{xpu::global_index<1>()};
  if (i >= g_vector.count()) { return; }

  stencil::energy::update_structure_factors(
    i,
    old_pos, new_pos,
    g_vector,
    sum_real, sum_imag
  );
}
#endif

} // namespace

inline fp_t initialize_reciprocal_energy(
  std::size_t num_g_vectors,
  const fp_t* RESTRICT g_weights,
  const fp_t* RESTRICT sum_real,
  const fp_t* RESTRICT sum_imag,
  fp_t* RESTRICT reciprocal_sum_scratch
) noexcept {
#if defined(XPU_CUDA)
  xpu::zero_n(reciprocal_sum_scratch, 1uz);

  dim3 initializeReciprocalEnergyThreads{256u};
  dim3 initializeReciprocalEnergyBlocks{
    xpu::block_per_dim(
      num_g_vectors, initializeReciprocalEnergyThreads.x
    )
  };

  cudaInitializeReciprocalEnergy<<<
    initializeReciprocalEnergyBlocks, initializeReciprocalEnergyThreads
  >>>(
    num_g_vectors,
    g_weights, sum_real, sum_imag,
    reciprocal_sum_scratch
  );
  xpu::cu_check(cudaGetLastError());

  auto reciprocal_sum{0.0_fp};
  xpu::copy_n(&reciprocal_sum, reciprocal_sum_scratch, 1uz);
  return reciprocal_sum;
#else
  static_cast<void>(reciprocal_sum_scratch);

  auto reciprocal_sum{0.0_fp};

  #pragma omp simd reduction(+ : reciprocal_sum)
  for (auto i = 0uz; i < num_g_vectors; ++i) {
    stencil::energy::initialize_reciprocal_energy(
      i,
      g_weights, sum_real, sum_imag,
      &reciprocal_sum
    );
  }

  return reciprocal_sum;
#endif
}

inline fp_t initialize_real_energy(
  fp_t L, fp_t half_L, fp_t alpha,
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos,
  fp_t* RESTRICT real_sum_scratch
) noexcept {
#if defined(XPU_CUDA)
  xpu::zero_n(real_sum_scratch, 1uz);

  dim3 initializeRealEnergyThreads{16u, 16u};
  dim3 initializeRealEnergyBlocks{
    xpu::block_per_dim(pos.count(), initializeRealEnergyThreads.x),
    xpu::block_per_dim(pos.count(), initializeRealEnergyThreads.y)
  };

  cudaInitializeRealEnergy<<<
    initializeRealEnergyBlocks, initializeRealEnergyThreads
  >>>(
    L, half_L, alpha,
    pos,
    real_sum_scratch
  );
  xpu::cu_check(cudaGetLastError());

  auto real_sum{0.0_fp};
  xpu::copy_n(&real_sum, real_sum_scratch, 1uz);
  return real_sum;
#else
  static_cast<void>(real_sum_scratch);

  auto real_sum{0.0_fp};

  for (auto i = 0uz; i < pos.count(); ++i) {
    xpu::array<fp_t, idx(Axis::NUM)> particle_pos{
      pos[idx(Axis::X)][i],
      pos[idx(Axis::Y)][i],
      pos[idx(Axis::Z)][i]
    };

    auto local_sum{0.0_fp};

    #pragma omp simd reduction(+ : local_sum)
    for (auto j = i + 1uz; j < pos.count(); ++j) {
      stencil::energy::initialize_real_energy(
        j,
        L, half_L, alpha,
        particle_pos, pos,
        &local_sum
      );
    }

    real_sum += local_sum;
  }

  return real_sum;
#endif
}

inline fp_t update_real_energy(
  std::size_t moved,
  fp_t L, fp_t half_L, fp_t alpha,
  xpu::array<fp_t, idx(Axis::NUM)> old_pos,
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos,
  fp_t* RESTRICT delta_scratch
) noexcept {
#if defined(XPU_CUDA)
  xpu::zero_n(delta_scratch, 1uz);

  dim3 updateRealEnergyThreads{256u};
  dim3 updateRealEnergyBlocks{
    xpu::block_per_dim(pos.count(), updateRealEnergyThreads.x)
  };

  cudaUpdateRealEnergy<<<
    updateRealEnergyBlocks, updateRealEnergyThreads
  >>>(
    moved,
    L, half_L, alpha,
    old_pos, pos,
    delta_scratch
  );
  xpu::cu_check(cudaGetLastError());

  auto delta{0.0_fp};
  xpu::copy_n(&delta, delta_scratch, 1uz);
  return delta;
#else
  static_cast<void>(delta_scratch);

  xpu::array<fp_t, idx(Axis::NUM)> new_pos{
    pos[idx(Axis::X)][moved],
    pos[idx(Axis::Y)][moved],
    pos[idx(Axis::Z)][moved]
  };

  auto delta{0.0_fp};

  #pragma omp simd reduction(+ : delta)
  for (auto i = 0uz; i < pos.count(); ++i) {
    stencil::energy::update_real_energy(
      i, moved,
      L, half_L, alpha,
      old_pos, new_pos, pos,
      &delta
    );
  }

  return delta;
#endif
}

inline fp_t kinetic_energy(
  xpu::soa_view<fp_t, idx(Derivatives::NUM)> derivatives,
  fp_t* RESTRICT kinetic_sum_scratch
) noexcept {
#if defined(XPU_CUDA)
  xpu::zero_n(kinetic_sum_scratch, 1uz);

  dim3 kineticEnergyThreads{256u};
  dim3 kineticEnergyBlocks{
    xpu::block_per_dim(derivatives.count(), kineticEnergyThreads.x)
  };

  cudaKineticEnergy<<<
    kineticEnergyBlocks, kineticEnergyThreads
  >>>(
    derivatives,
    kinetic_sum_scratch
  );
  xpu::cu_check(cudaGetLastError());

  auto kinetic_sum{0.0_fp};
  xpu::copy_n(&kinetic_sum, kinetic_sum_scratch, 1uz);
  return -0.5_fp * kinetic_sum;
#else
  static_cast<void>(kinetic_sum_scratch);

  auto kinetic_sum{0.0_fp};

  #pragma omp simd reduction(+ : kinetic_sum)
  for (auto i = 0uz; i < derivatives.count(); ++i) {
    stencil::energy::kinetic_energy(
      i,
      derivatives,
      &kinetic_sum
    );
  }

  return -0.5_fp * kinetic_sum;
#endif
}

inline void initialize_structure_factors(
  xpu::soa_view<fp_t, idx(Axis::NUM)> g_vector,
  xpu::soa_view<fp_t, idx(Axis::NUM)> pos,
  fp_t* RESTRICT sum_real,
  fp_t* RESTRICT sum_imag
) noexcept {
  xpu::zero_n(sum_real, g_vector.count());
  xpu::zero_n(sum_imag, g_vector.count());

#if defined(XPU_CUDA)
  dim3 initializeStructureFactorsThreads{16u, 16u};
  dim3 initializeStructureFactorsBlocks{
    xpu::block_per_dim(
      g_vector.count(), initializeStructureFactorsThreads.x
    ),
    xpu::block_per_dim(
      pos.count(), initializeStructureFactorsThreads.y
    )
  };

  cudaInitializeStructureFactors<<<
    initializeStructureFactorsBlocks, initializeStructureFactorsThreads
  >>>(
    g_vector, pos,
    sum_real, sum_imag
  );
  xpu::cu_check(cudaGetLastError());
#else
  for (auto i = 0uz; i < g_vector.count(); ++i) {
    auto& real_sum{sum_real[i]};
    auto& imag_sum{sum_imag[i]};

    #pragma omp simd reduction(+ : real_sum, imag_sum)
    for (auto j = 0uz; j < pos.count(); ++j) {
      stencil::energy::initialize_structure_factors(
        i, j,
        g_vector, pos,
        &real_sum, &imag_sum
      );
    }
  }
#endif
}

inline void update_structure_factors(
  xpu::array<fp_t, idx(Axis::NUM)> old_pos,
  xpu::array<fp_t, idx(Axis::NUM)> new_pos,
  xpu::soa_view<fp_t, idx(Axis::NUM)> g_vector,
  fp_t* RESTRICT sum_real,
  fp_t* RESTRICT sum_imag
) noexcept {
#if defined(XPU_CUDA)
  dim3 updateStructureFactorsThreads{256u};
  dim3 updateStructureFactorsBlocks{
    xpu::block_per_dim(
      g_vector.count(), updateStructureFactorsThreads.x
    )
  };

  cudaUpdateStructureFactors<<<
    updateStructureFactorsBlocks, updateStructureFactorsThreads
  >>>(
    old_pos, new_pos,
    g_vector,
    sum_real, sum_imag
  );
  xpu::cu_check(cudaGetLastError());
#else
  #pragma omp simd
  for (auto i = 0uz; i < g_vector.count(); ++i) {
    stencil::energy::update_structure_factors(
      i,
      old_pos, new_pos,
      g_vector,
      sum_real, sum_imag
    );
  }
#endif
}

inline fp_t initialize_reciprocal_energy(
  EnergyTracker::View energy
) noexcept {
  return initialize_reciprocal_energy(
    energy.num_g_vectors,
    energy.g_weights,
    energy.sum_real,
    energy.sum_imag,
    energy.reduction_scratch
  );
}

inline fp_t initialize_real_energy(
  EnergyTracker::View energy,
  Particles::View particles
) noexcept {
  return initialize_real_energy(
    energy.box_length,
    0.5_fp * energy.box_length,
    energy.ewald_alpha,
    particles.pos,
    energy.reduction_scratch
  );
}

inline fp_t update_real_energy(
  EnergyTracker::View energy,
  std::size_t moved,
  xpu::array<fp_t, idx(Axis::NUM)> old_pos,
  Particles::View particles
) noexcept {
  return update_real_energy(
    moved,
    energy.box_length,
    0.5_fp * energy.box_length,
    energy.ewald_alpha,
    old_pos,
    particles.pos,
    energy.reduction_scratch
  );
}

inline fp_t kinetic_energy(
  EnergyTracker::View energy,
  Particles::View particles
) noexcept {
  return kinetic_energy(
    particles.derivatives,
    energy.reduction_scratch
  );
}

inline void initialize_structure_factors(
  EnergyTracker::View energy,
  Particles::View particles
) noexcept {
  initialize_structure_factors(
    energy.g_vector,
    particles.pos,
    energy.sum_real,
    energy.sum_imag
  );
}

inline void update_structure_factors(
  EnergyTracker::View energy,
  xpu::array<fp_t, idx(Axis::NUM)> old_pos,
  xpu::array<fp_t, idx(Axis::NUM)> new_pos
) noexcept {
  update_structure_factors(
    old_pos,
    new_pos,
    energy.g_vector,
    energy.sum_real,
    energy.sum_imag
  );
}

} // namespace kernel::energy
} // namespace kernel
