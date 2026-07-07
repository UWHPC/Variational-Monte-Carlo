#include "slater_plane_wave.cuh"
#include "../utilities/matrix.hpp"
#include "particles/particles.cuh"
#include "utilities/aligned_soa.cuh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <vector>

// Real plane-wave basis from complex exponentials
// The complex basis e^{ik·r} = cos(k·r) + i·sin(k·r) gives two
// linearly independent real orbitals per k-vector:
//   φ_cos(r) = cos(k·r)
//   φ_sin(r) = sin(k·r)
//
// Since cos(k·r) = cos(-k·r) and sin(k·r) = -sin(-k·r),
// the ±k pair maps to the same {cos, sin} pair. We deduplicate
// by keeping only the canonical +n representative (first nonzero
// component positive).
//
// k = 0 is special: sin(0) = 0 for all r, so it contributes
// only one orbital (cos(0) = 1).
//
// Orbital count: N = 1 + 2 × (number of nonzero canonical k-vectors)
// Closed shells: N = 1, 7, 19, 27, 33, 57, ...
SlaterPlaneWave::SlaterPlaneWave(const Particles& particles, real_t box_lengthL)
: num_orbitals_{particles.size()}
, matrix_row_stride_{AlignedSoA<real_t>::round_up(particles.size())}
, matrix_size_{matrix_row_stride_ * particles.size()}
, box_length_{box_lengthL}
, orbital_k_index_(particles.size(), NUM_ORB_K)
, orbital_type_(particles.size(), NUM_ORB_TYPE)
, int_vec_{particles.size(), NUM_INT_VECTORS}
, fp_vec_{particles.size(), NUM_DOUBLE_VECTORS}
, trig_cache_{}
, matrices_{matrix_row_stride_ * particles.size(), NUM_MATRIX} {

  const std::size_t N{num_orbitals()};
  const std::size_t num_particles{particles.size()};

  // Generate deduplicated canonical n-vectors
  // Only keep canonical representatives: first nonzero component positive
  // k=0 produces 1 orbital (cos = 1)
  // Each nonzero k produces 2 orbitals (cos, sin)

  struct nVectorCandidate {
    int n_cand_x;
    int n_cand_y;
    int n_cand_z;
    int n_mag_sq;
  };

  // Increase vector size to be safe
  const int N_MAX{static_cast<int>(vmc::ceil(vmc::cbrt(static_cast<real_t>(N)))) + 2};
  const std::size_t side{static_cast<std::size_t>((2 * N_MAX + 1))};

  // Vector for all possible n-vector candidates
  std::vector<nVectorCandidate> n_candidates{};
  n_candidates.reserve(side * side * side);

  // Fill candidate vector with all possible states from [-N, N]
  for (int new_x = -N_MAX; new_x <= N_MAX; ++new_x) {
    for (int new_y = -N_MAX; new_y <= N_MAX; ++new_y) {
      for (int new_z = -N_MAX; new_z <= N_MAX; ++new_z) {
        if (!is_canonical(new_x, new_y, new_z))
          continue;

        const int new_mag_sq{
          new_x * new_x +
          new_y * new_y +
          new_z * new_z
        };
        n_candidates.emplace_back(new_x, new_y, new_z, new_mag_sq);
      }
    }
  }

  // Sort n-vector states to go from smallest magnitude to largest
  std::sort(
    n_candidates.begin(),
    n_candidates.end(),
    [](const nVectorCandidate& a, const nVectorCandidate& b) {
      return
        std::tie(a.n_mag_sq, a.n_cand_x, a.n_cand_y, a.n_cand_z) <
        std::tie(b.n_mag_sq, b.n_cand_x, b.n_cand_y, b.n_cand_z);
    }
  );

  // Assign orbitals
  // Orbital 0: k=0 -> cos(0 dot r) = cos(0) = 1
  // For each nonzero canonical k: orbital 2m-1 -> cos(k dot r), orbital 2m -> sin(k dot r)
  const auto nv{n_vector().align()};

  auto* orb_k_idx{orbital_k_index()};
  auto* orb_type{orbital_type()};

  std::size_t orb_idx{};
  std::size_t k_idx{};

  for (std::size_t c = 0; c < n_candidates.size() && orb_idx < N; ++c) {
    const auto& cand{n_candidates[c]};
    const bool mag_not_zero{cand.n_mag_sq != 0};

    nv.x_[k_idx] = cand.n_cand_x;
    nv.y_[k_idx] = cand.n_cand_y;
    nv.z_[k_idx] = cand.n_cand_z;

    // Cos orbital (every k-vector gets one)
    orb_k_idx[orb_idx] = k_idx;
    orb_type[orb_idx] = 0;
    ++orb_idx;

    // Sin orbital (only for non-zero k - sin(0) = 0 is singular)
    if (mag_not_zero && orb_idx < N) {
      orb_k_idx[orb_idx] = k_idx;
      orb_type[orb_idx] = 1;
      ++orb_idx;
    }

    ++k_idx;
  }

  num_unique_k_ = k_idx;
  trig_row_stride_ = AlignedSoA<real_t>::round_up(k_idx);

  const auto kv{k_vector().align()};

  const real_t inv_L{1.0_r / box_length()};

  // Follows the calculation: K = (2pi/L) * n;
  #pragma omp simd
  for (std::size_t i = 0; i < k_idx; ++i) {
    kv.x_[i] = 2.0_r * std::numbers::pi_v<real_t> * inv_L * static_cast<real_t>(nv.x_[i]);
    kv.y_[i] = 2.0_r * std::numbers::pi_v<real_t> * inv_L * static_cast<real_t>(nv.y_[i]);
    kv.z_[i] = 2.0_r * std::numbers::pi_v<real_t> * inv_L * static_cast<real_t>(nv.z_[i]);
  }

  trig_cache_ = AlignedSoA<real_t>(num_particles * trig_row_stride(), NUM_TRIG_ARRAYS);
  trig_scratch_ = AlignedSoA<real_t>(num_unique_k(), NUM_SCRATCH_TRIG);

  std::size_t trig_size{trig_cache_.num_elements()};
  std::fill_n(sin_cache(), trig_size, 0.0_r);
  std::fill_n(cos_cache(), trig_size, 0.0_r);
};


real_t SlaterPlaneWave::determinant_ratio(
  std::size_t particle,
  const real_t* new_row
) const noexcept {
  const std::size_t N{num_orbitals()};
  const std::size_t S{matrix_row_stride()};
  const real_t* RESTRICT inv_det{inv_determinant()};
  ASSUME_ALIGNED(inv_det, SIMD_BYTES);

  real_t ratio{};
  #pragma omp simd reduction(+ : ratio)
  for (std::size_t j = 0; j < N; ++j) {
    ratio += new_row[j] * inv_det[particle * S + j];
  }

  return ratio;
}

void SlaterPlaneWave::accept_move(
  std::size_t particle,
  const real_t* new_row,
  real_t ratio
) noexcept {
  const std::size_t N{num_orbitals()};
  const std::size_t S{matrix_row_stride()};

  real_t* RESTRICT inv_det{inv_determinant()};
  real_t* RESTRICT det_matrix{determinant()};
  real_t* RESTRICT inv_d_col{this->inv_d_col()};

  ASSUME_ALIGNED(inv_det, SIMD_BYTES);
  ASSUME_ALIGNED(det_matrix, SIMD_BYTES);
  ASSUME_ALIGNED(inv_d_col, SIMD_BYTES);

  const real_t inv_ratio{1.0_r / ratio};

  if (!std::isfinite(inv_ratio)) {
    return;
  }

  const std::size_t p_offset{particle * S}; // Pre-calculate particle row offset

  // Cache particle row column j for inv_D before changing
  std::memcpy(inv_d_col, &inv_det[p_offset], N * sizeof(real_t));

  // Follows Sherman-Morrison update (branchless)
  for (std::size_t k = 0; k < N; ++k) {
    // Skip the special row for now
    if (k == particle)
      continue;

    const std::size_t k_offset{k * S}; // Pre-calculate row offset
    real_t s_k{};

    // Compiler can now easily auto-vectorize this
    #pragma omp simd reduction(+ : s_k)
    for (std::size_t m = 0; m < N; ++m) {
      s_k += new_row[m] * inv_det[k_offset + m];
    }

    const real_t factor{s_k * inv_ratio};

    // Not vectorized: loop-carried data dependency
    #pragma omp simd
    for (std::size_t j = 0; j < N; ++j) {
      inv_det[k_offset + j] -= inv_d_col[j] * factor;
    }
  }

  // Special case: handle the particle row outside the loop
  #pragma omp simd
  for (std::size_t j = 0; j < N; ++j) {
    inv_det[p_offset + j] = inv_d_col[j] * inv_ratio;
  }

  // Patch row `particle` of D to match the new positions:
  std::memcpy(&det_matrix[p_offset], new_row, N * sizeof(real_t));
}

void SlaterPlaneWave::add_derivatives(
  real_t* RESTRICT grad_x,
  real_t* RESTRICT grad_y,
  real_t* RESTRICT grad_z,
  real_t* RESTRICT laplacian
) const noexcept {
  const std::size_t N{num_orbitals()};
  const std::size_t S{matrix_row_stride()};

  const auto kv{k_vector().align()};

  const auto& k_index{orbital_k_index()};
  const auto& o_type{orbital_type()};

  const real_t* RESTRICT inv_det{inv_determinant()};
  const real_t* RESTRICT cos_cache{this->cos_cache()};
  const real_t* RESTRICT sin_cache{this->sin_cache()};

  ASSUME_ALIGNED(inv_det, SIMD_BYTES);
  ASSUME_ALIGNED(cos_cache, SIMD_BYTES);
  ASSUME_ALIGNED(sin_cache, SIMD_BYTES);

  ASSUME_ALIGNED(grad_x, SIMD_BYTES);
  ASSUME_ALIGNED(grad_y, SIMD_BYTES);
  ASSUME_ALIGNED(grad_z, SIMD_BYTES);
  ASSUME_ALIGNED(laplacian, SIMD_BYTES);

  for (std::size_t particle = 0; particle < N; ++particle) {
    real_t d_log_det_dx{}, d_log_det_dy{}, d_log_det_dz{};
    real_t laplace_det_term{};

    const std::size_t offset{particle * trig_row_stride()};
    const std::size_t p_offset{particle * S};

    // Added reduction clauses so the compiler can safely vectorize the accumulators
    // Not vectorized: loop-carried data dependency
    #pragma omp simd reduction(+ : d_log_det_dx, d_log_det_dy, d_log_det_dz, laplace_det_term)
    for (std::size_t orbital = 0; orbital < N; ++orbital) {
      const std::size_t k_idx{k_index[orbital]};

      const real_t k_x_orbital{kv.x_[k_idx]};
      const real_t k_y_orbital{kv.y_[k_idx]};
      const real_t k_z_orbital{kv.z_[k_idx]};

      const real_t k_sq{
        k_x_orbital * k_x_orbital +
        k_y_orbital * k_y_orbital +
        k_z_orbital * k_z_orbital
      };

      const real_t type{static_cast<real_t>(o_type[orbital])};

      const real_t cos_term{cos_cache[offset + k_idx]};
      const real_t sin_term{sin_cache[offset + k_idx]};

      const real_t grad_factor{-sin_term + type * (sin_term + cos_term)};
      const real_t lap_factor{-cos_term + type * (cos_term - sin_term)};

      const real_t weight{inv_det[p_offset + orbital]};

      d_log_det_dx += weight * k_x_orbital * grad_factor;
      d_log_det_dy += weight * k_y_orbital * grad_factor;
      d_log_det_dz += weight * k_z_orbital * grad_factor;
      laplace_det_term += weight * k_sq * lap_factor;
    }

    const real_t grad_sq{
      d_log_det_dx * d_log_det_dx +
      d_log_det_dy * d_log_det_dy +
      d_log_det_dz * d_log_det_dz
    };

    grad_x[particle] += d_log_det_dx;
    grad_y[particle] += d_log_det_dy;
    grad_z[particle] += d_log_det_dz;

    laplacian[particle] += (laplace_det_term - grad_sq);
  }
}
