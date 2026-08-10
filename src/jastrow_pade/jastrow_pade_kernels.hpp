#pragma once

#include <xpu/xpu.hpp>
#include "../utilities/components.hpp"
#include "../utilities/macros.hpp"

namespace stencil {
namespace jpade {

namespace {

CUDA_CALLABLE
inline real_t evaluate_pair_value(
  const xpu::array<real_t, idx(Axis::NUM)>& particle_pos,
  std::size_t other,
  real_t L, real_t half_L,
  real_t a, real_t b,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos
) noexcept {
  xpu::array<real_t, idx(Axis::NUM)> displacement{};

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

  return a * distance / (1.0_r + b * distance);
}

CUDA_CALLABLE
inline void evaluate_derivative_factors(
  const real_t distance,
  const real_t a,
  const real_t b,
  const real_t neg_two_ab,
  real_t* gradient_factor,
  real_t* laplacian_factor
) noexcept {
  if (distance < EPSILON) {
    *gradient_factor = 0.0_r;
    *laplacian_factor = 0.0_r;
    return;
  }

  const auto inverse_distance{1.0_r / distance};
  const auto inverse_denominator{1.0_r / (1.0_r + b * distance)};
  const auto inverse_denominator_squared{
    inverse_denominator * inverse_denominator
  };

  const auto first_deriv{a * inverse_denominator_squared};
  const auto second_deriv{
    neg_two_ab * inverse_denominator_squared * inverse_denominator
  };

  *gradient_factor = first_deriv * inverse_distance;
  *laplacian_factor = second_deriv + 2.0_r * first_deriv * inverse_distance;
}

CUDA_CALLABLE
inline void evaluate_pair_derivatives(
  const xpu::array<real_t, idx(Axis::NUM)>& particle_pos,
  std::size_t other,
  real_t L, real_t half_L,
  real_t a, real_t b, real_t neg_two_ab,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  xpu::array<real_t, idx(Axis::NUM)>* gradient,
  real_t* laplacian
) noexcept {
  xpu::array<real_t, idx(Axis::NUM)> displacement{};

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

  auto gradient_factor{0.0_r};

  evaluate_derivative_factors(
    distance, a, b, neg_two_ab,
    &gradient_factor, laplacian
  );

  for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
    (*gradient)[axis] = gradient_factor * displacement[axis];
  }
}

} // namespace

CUDA_CALLABLE
inline void value(
  std::size_t other,
  real_t L, real_t half_L,
  real_t a, real_t b,
  const xpu::array<real_t, idx(Axis::NUM)>& particle_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  real_t* jastrow_value
) noexcept {
  const auto pair_value{
    evaluate_pair_value(
      particle_pos, other,
      L, half_L, a, b,
      pos
    )
  };

#if defined(__CUDA_ARCH__)
  atomicAdd(jastrow_value, pair_value);
#else
  *jastrow_value += pair_value;
#endif
}

CUDA_CALLABLE
inline void delta_value(
  std::size_t other,
  std::size_t moved,
  real_t L, real_t half_L,
  real_t a, real_t b,
  const xpu::array<real_t, idx(Axis::NUM)>& old_pos,
  const xpu::array<real_t, idx(Axis::NUM)>& new_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  real_t* delta
) noexcept {
  if (other == moved) { return; }

  const auto value_old{
    evaluate_pair_value(
      old_pos, other,
      L, half_L, a, b,
      pos
    )
  };
  const auto value_new{
    evaluate_pair_value(
      new_pos, other,
      L, half_L, a, b,
      pos
    )
  };

#if defined(__CUDA_ARCH__)
  atomicAdd(delta, value_new - value_old);
#else
  *delta += value_new - value_old;
#endif
}

CUDA_CALLABLE
inline void compute_derivatives(
  std::size_t i,
  std::size_t moved,
  real_t L, real_t half_L,
  real_t a, real_t b, real_t neg_two_ab,
  const xpu::array<real_t, idx(Axis::NUM)>& old_pos,
  const xpu::array<real_t, idx(Axis::NUM)>& new_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  xpu::soa_view<real_t, idx(Derivatives::NUM)> derivatives,
  real_t* moved_gx, real_t* moved_gy, real_t* moved_gz,real_t* moved_lap
) noexcept {
  if (i == moved) { return; }

  xpu::array<real_t, idx(Axis::NUM)> gradient_old{};
  xpu::array<real_t, idx(Axis::NUM)> gradient_new{};

  auto laplacian_old{0.0_r};
  auto laplacian_new{0.0_r};

  evaluate_pair_derivatives(
    old_pos, i,
    L, half_L, a, b, neg_two_ab,
    pos, &gradient_old, &laplacian_old
  );
  evaluate_pair_derivatives(
    new_pos, i,
    L, half_L, a, b, neg_two_ab,
    pos, &gradient_new, &laplacian_new
  );

  for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
    derivatives[axis][i] += gradient_old[axis] - gradient_new[axis];
  }
  derivatives[idx(Derivatives::LAP)][i] += laplacian_new - laplacian_old;

#if defined(__CUDA_ARCH__)
  atomicAdd(moved_gx, gradient_new[idx(Axis::X)] - gradient_old[idx(Axis::X)]);
  atomicAdd(moved_gy, gradient_new[idx(Axis::Y)] - gradient_old[idx(Axis::Y)]);
  atomicAdd(moved_gz, gradient_new[idx(Axis::Z)] - gradient_old[idx(Axis::Z)]);
  atomicAdd(moved_lap, laplacian_new - laplacian_old);
#else
  *moved_gx += gradient_new[idx(Axis::X)] - gradient_old[idx(Axis::X)];
  *moved_gy += gradient_new[idx(Axis::Y)] - gradient_old[idx(Axis::Y)];
  *moved_gz += gradient_new[idx(Axis::Z)] - gradient_old[idx(Axis::Z)];
  *moved_lap += laplacian_new - laplacian_old;
#endif
}

CUDA_CALLABLE
inline void add_derivatives(
  std::size_t other,
  real_t L, real_t half_L,
  real_t a, real_t b, real_t neg_two_ab,
  const xpu::array<real_t, idx(Axis::NUM)>& particle_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  real_t* gradient_x, real_t* gradient_y, real_t* gradient_z,
  real_t* laplacian
) noexcept {
  xpu::array<real_t, idx(Axis::NUM)> pair_gradient{};
  auto pair_laplacian{0.0_r};

  evaluate_pair_derivatives(
    particle_pos, other,
    L, half_L, a, b, neg_two_ab,
    pos, &pair_gradient, &pair_laplacian
  );

#if defined(__CUDA_ARCH__)
  atomicAdd(gradient_x, pair_gradient[idx(Axis::X)]);
  atomicAdd(gradient_y, pair_gradient[idx(Axis::Y)]);
  atomicAdd(gradient_z, pair_gradient[idx(Axis::Z)]);
  atomicAdd(laplacian, pair_laplacian);
#else
  *gradient_x += pair_gradient[idx(Axis::X)];
  *gradient_y += pair_gradient[idx(Axis::Y)];
  *gradient_z += pair_gradient[idx(Axis::Z)];
  *laplacian += pair_laplacian;
#endif
}

} // namespace stencil::jpade
} // namespace stencil

namespace kernel {
namespace jpade {

namespace {

#if defined(XPU_CUDA)
__global__
void cudaValue(
  real_t L, real_t half_L,
  real_t a, real_t b,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  real_t* jastrow_value
) {
  const auto [i, j]{xpu::global_index<2>()};
  if (i >= pos.count() || j >= pos.count()) { return; }
  if (i >= j) { return; }

  xpu::array<real_t, idx(Axis::NUM)> particle_pos{
    pos[idx(Axis::X)][i],
    pos[idx(Axis::Y)][i],
    pos[idx(Axis::Z)][i]
  };

  stencil::jpade::value(
    j,
    L, half_L, a, b,
    particle_pos, pos,
    jastrow_value
  );
}

__global__
void cudaDeltaValue(
  std::size_t moved,
  real_t L, real_t half_L,
  real_t a, real_t b,
  xpu::array<real_t, idx(Axis::NUM)> old_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  real_t* delta
) {
  const auto [i]{xpu::global_index<1>()};
  if (i >= pos.count()) { return; }

  xpu::array<real_t, idx(Axis::NUM)> new_pos{
    pos[idx(Axis::X)][moved],
    pos[idx(Axis::Y)][moved],
    pos[idx(Axis::Z)][moved]
  };

  stencil::jpade::delta_value(
    i, moved,
    L, half_L, a, b,
    old_pos, new_pos, pos,
    delta
  );
}

__global__
void cudaComputeDerivatives(
  std::size_t moved,
  real_t L, real_t half_L,
  real_t a, real_t b, real_t neg_two_ab,
  xpu::array<real_t, idx(Axis::NUM)> old_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  xpu::soa_view<real_t, idx(Derivatives::NUM)> derivatives
) {
  const auto [i]{xpu::global_index<1>()};
  if (i >= pos.count()) { return; }

  xpu::array<real_t, idx(Axis::NUM)> new_pos{
    pos[idx(Axis::X)][moved],
    pos[idx(Axis::Y)][moved],
    pos[idx(Axis::Z)][moved]
  };

  stencil::jpade::compute_derivatives(
    i, moved,
    L, half_L, a, b, neg_two_ab,
    old_pos, new_pos, pos, derivatives,
    &derivatives[idx(Derivatives::GRAD_X)][moved],
    &derivatives[idx(Derivatives::GRAD_Y)][moved],
    &derivatives[idx(Derivatives::GRAD_Z)][moved],
    &derivatives[idx(Derivatives::LAP)][moved]
  );
}

__global__
void cudaAddDerivatives(
  real_t L, real_t half_L,
  real_t a, real_t b, real_t neg_two_ab,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  xpu::soa_view<real_t, idx(Derivatives::NUM)> derivatives
) {
  const auto [i, j]{xpu::global_index<2>()};
  if (i >= pos.count() || j >= pos.count()) { return; }
  if (i == j) { return; }

  xpu::array<real_t, idx(Axis::NUM)> particle_pos{
    pos[idx(Axis::X)][i],
    pos[idx(Axis::Y)][i],
    pos[idx(Axis::Z)][i]
  };

  stencil::jpade::add_derivatives(
    j,
    L, half_L, a, b, neg_two_ab,
    particle_pos, pos,
    &derivatives[idx(Derivatives::GRAD_X)][i],
    &derivatives[idx(Derivatives::GRAD_Y)][i],
    &derivatives[idx(Derivatives::GRAD_Z)][i],
    &derivatives[idx(Derivatives::LAP)][i]
  );
}
#endif

} // namespace

inline real_t value(
  real_t L, real_t half_L,
  real_t a, real_t b,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos
) {
#if defined(XPU_CUDA)
  xpu::buffer<real_t> jastrow_value{1};

  dim3 valueThreads{16, 16};
  dim3 valueBlocks{
    xpu::block_per_dim(pos.count(), valueThreads.x),
    xpu::block_per_dim(pos.count(), valueThreads.y)
  };

  cudaValue<<<
    valueBlocks, valueThreads
  >>>(
    L, half_L, a, b,
    pos, jastrow_value.data()
  );
  xpu::cu_check(cudaGetLastError());

  real_t value_host{};
  xpu::cu_check(cudaMemcpy(
    &value_host, jastrow_value.data(),
    sizeof(real_t), cudaMemcpyDeviceToHost
  ));
  return value_host;
#else
  auto jastrow_value{0.0_r};

  for (auto i = 0uz; i < pos.count(); ++i) {
    xpu::array<real_t, idx(Axis::NUM)> particle_pos{
      pos[idx(Axis::X)][i],
      pos[idx(Axis::Y)][i],
      pos[idx(Axis::Z)][i]
    };

    auto local_value{0.0_r};

    #pragma omp simd reduction(+ : local_value)
    for (auto j = i + 1uz; j < pos.count(); ++j) {
      stencil::jpade::value(
        j,
        L, half_L, a, b,
        particle_pos, pos,
        &local_value
      );
    }

    jastrow_value += local_value;
  }

  return jastrow_value;
#endif
}

inline real_t delta_value(
  std::size_t moved,
  real_t L, real_t half_L,
  real_t a, real_t b,
  xpu::array<real_t, idx(Axis::NUM)> old_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos
) {
#if defined(XPU_CUDA)
  xpu::buffer<real_t> delta{1};

  dim3 deltaValueThreads{256};
  dim3 deltaValueBlocks{
    xpu::block_per_dim(pos.count(), deltaValueThreads.x)
  };

  cudaDeltaValue<<<
    deltaValueBlocks, deltaValueThreads
  >>>(
    moved,
    L, half_L, a, b,
    old_pos, pos, delta.data()
  );
  xpu::cu_check(cudaGetLastError());

  real_t delta_host{};
  xpu::cu_check(cudaMemcpy(
    &delta_host, delta.data(),
    sizeof(real_t), cudaMemcpyDeviceToHost
  ));
  return delta_host;
#else
  xpu::array<real_t, idx(Axis::NUM)> new_pos{
    pos[idx(Axis::X)][moved],
    pos[idx(Axis::Y)][moved],
    pos[idx(Axis::Z)][moved]
  };

  auto delta{0.0_r};

  #pragma omp simd reduction(+ : delta)
  for (auto i = 0uz; i < pos.count(); ++i) {
    stencil::jpade::delta_value(
      i, moved,
      L, half_L, a, b,
      old_pos, new_pos, pos,
      &delta
    );
  }

  return delta;
#endif
}

inline void compute_derivatives(
  std::size_t moved,
  real_t L, real_t half_L,
  real_t a, real_t b, real_t neg_two_ab,
  xpu::array<real_t, idx(Axis::NUM)> old_pos,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  xpu::soa_view<real_t, idx(Derivatives::NUM)> derivatives
) {
#if defined(XPU_CUDA)
  dim3 computeDerivativesThreads{256};
  dim3 computeDerivativesBlocks{
    xpu::block_per_dim(pos.count(), computeDerivativesThreads.x)
  };

  cudaComputeDerivatives<<<
    computeDerivativesBlocks, computeDerivativesThreads
  >>>(
    moved, 
    L, half_L, a, b, neg_two_ab,
    old_pos, pos, derivatives
  );
  xpu::cu_check(cudaGetLastError());
#else
  xpu::array<real_t, idx(Axis::NUM)> new_pos{
    pos[idx(Axis::X)][moved],
    pos[idx(Axis::Y)][moved],
    pos[idx(Axis::Z)][moved]
  };

  auto& moved_gx{derivatives[idx(Derivatives::GRAD_X)][moved]};
  auto& moved_gy{derivatives[idx(Derivatives::GRAD_Y)][moved]};
  auto& moved_gz{derivatives[idx(Derivatives::GRAD_Z)][moved]};
  auto& moved_lap{derivatives[idx(Derivatives::LAP)][moved]};

  #pragma omp simd reduction(+ : moved_gx, moved_gy, moved_gz, moved_lap)
  for (auto i = 0uz; i < pos.count(); ++i) {
    stencil::jpade::compute_derivatives(
      i, moved,
      L, half_L, a, b, neg_two_ab,
      old_pos, new_pos, pos, derivatives,
      &moved_gx, &moved_gy, &moved_gz, &moved_lap
    );
  }
#endif
}

inline void add_derivatives(
  real_t L, real_t half_L,
  real_t a, real_t b, real_t neg_two_ab,
  xpu::soa_view<real_t, idx(Axis::NUM)> pos,
  xpu::soa_view<real_t, idx(Derivatives::NUM)> derivatives
) {
#if defined(XPU_CUDA)
  dim3 addDerivativesThreads{16, 16};
  dim3 addDerivativesBlocks{
    xpu::block_per_dim(pos.count(), addDerivativesThreads.x),
    xpu::block_per_dim(pos.count(), addDerivativesThreads.y)
  };

  cudaAddDerivatives<<<
    addDerivativesBlocks, addDerivativesThreads
  >>>(
    L, half_L, a, b, neg_two_ab,
    pos, derivatives
  );
  xpu::cu_check(cudaGetLastError());
#else
  for (auto i = 0uz; i < pos.count(); ++i) {
    xpu::array<real_t, idx(Axis::NUM)> particle_pos{
      pos[idx(Axis::X)][i],
      pos[idx(Axis::Y)][i],
      pos[idx(Axis::Z)][i]
    };

    auto& gradient_x{derivatives[idx(Derivatives::GRAD_X)][i]};
    auto& gradient_y{derivatives[idx(Derivatives::GRAD_Y)][i]};
    auto& gradient_z{derivatives[idx(Derivatives::GRAD_Z)][i]};
    auto& laplacian{derivatives[idx(Derivatives::LAP)][i]};

    #pragma omp simd reduction(+ : gradient_x, gradient_y, gradient_z, laplacian)
    for (auto j = 0uz; j < pos.count(); ++j) {
      stencil::jpade::add_derivatives(
        j,
        L, half_L, a, b, neg_two_ab,
        particle_pos, pos,
        &gradient_x, &gradient_y, &gradient_z, &laplacian
      );
    }
  }
#endif
}

} // namespace kernel::jpade
} // namespace kernel
