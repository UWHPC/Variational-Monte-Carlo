#include "slater_plane_wave.cuh"
#include "../utilities/matrix.hpp"

#ifdef VMC_CUDA_BACKEND
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <numbers>
#include <tuple>

namespace {

__global__
void cudaFindNVecCandidates(
  int n_max,
  int* RESTRICT nv_x,
  int* RESTRICT nv_y,
  int* RESTRICT nv_z,
  int* RESTRICT nv_mag_sq,
  std::size_t* RESTRICT counter
) {
  const auto [raw_i, raw_j, raw_k]{vmc::cudaThreadIdx<3>()};
  const auto side{static_cast<std::size_t>(2 * n_max + 1)};
  if (raw_i >= side || raw_j >= side || raw_k >= side) { return; }

  const auto i{static_cast<int>(raw_i) - n_max};
  const auto j{static_cast<int>(raw_j) - n_max};
  const auto k{static_cast<int>(raw_k) - n_max};

  if (!is_canonical(i, j, k)) { return; }

  const auto c{atomicAdd(counter, std::size_t{1})};
  nv_x[c] = i;
  nv_y[c] = j;
  nv_z[c] = k;
  nv_mag_sq[c] = (
    i * i +
    j * j +
    k * k
  );
}

__global__
void cudaAssignOrbitals(
  std::size_t N, std::size_t num_unique_k, real_t two_pi_inv_L,
  const std::size_t* RESTRICT order,
  const int* RESTRICT cand_x, const int* RESTRICT cand_y, const int* RESTRICT cand_z,
  int* RESTRICT nv_x, int* RESTRICT nv_y, int* RESTRICT nv_z,
  real_t* RESTRICT k_x, real_t* RESTRICT k_y, real_t* RESTRICT k_z,
  std::size_t* RESTRICT orb_k_idx, std::uint8_t* RESTRICT orb_type
) {
  const auto [i]{vmc::cudaThreadIdx<1>()};
  if (i >= num_unique_k) { return; }

  const std::size_t src{order[i]};
  const int nx{cand_x[src]};
  const int ny{cand_y[src]};
  const int nz{cand_z[src]};

  nv_x[i] = nx;
  nv_y[i] = ny;
  nv_z[i] = nz;

  k_x[i] = two_pi_inv_L * static_cast<real_t>(nx);
  k_y[i] = two_pi_inv_L * static_cast<real_t>(ny);
  k_z[i] = two_pi_inv_L * static_cast<real_t>(nz);

  const std::size_t cos_orb{(i == 0) ? 0 : 2 * i - 1};
  orb_k_idx[cos_orb] = i;
  orb_type[cos_orb] = 0;

  const std::size_t sin_orb{2 * i};
  if (i != 0 && sin_orb < N) {
    orb_k_idx[sin_orb] = i;
    orb_type[sin_orb] = 1;
  }
}

} // namespace
#else
#include "../particles/particles.cuh"
#include "../utilities/aligned_soa.cuh"

#include <algorithm>
#include <cstddef>
#include <numbers>
#include <tuple>
#include <vector>
#endif

void SlaterPlaneWave::initialize(const Particles& particles) {
#ifdef VMC_CUDA_BACKEND
  const std::size_t N{this->num_orbitals()};
  const std::size_t num_particles{particles.size()};

  enum nVectorCandidate : std::size_t {
    CAND_X,
    CAND_Y,
    CAND_Z,
    CAND_MAG_SQ,
    N_VEC_ARR
  };

  const int n_max{static_cast<int>(vmc::ceil(vmc::cbrt(static_cast<real_t>(N)))) + 2};
  const std::size_t side{static_cast<std::size_t>(2 * n_max + 1)};
  const std::size_t max_candidates{side * side * side};

  AlignedSoA<int> n_vec_tmp{max_candidates, N_VEC_ARR};
  AlignedSoA<std::size_t> counter{1, 1};

  dim3 findNVecCandidatesThreads(8, 8, 8);
  dim3 findNVecCandidatesBlocks(
    vmc::cudaNumBlocks(side, findNVecCandidatesThreads.x),
    vmc::cudaNumBlocks(side, findNVecCandidatesThreads.y),
    vmc::cudaNumBlocks(side, findNVecCandidatesThreads.z)
  );
  cudaFindNVecCandidates<<<findNVecCandidatesBlocks, findNVecCandidatesThreads>>>(
    n_max,
    n_vec_tmp[CAND_X], n_vec_tmp[CAND_Y], n_vec_tmp[CAND_Z], n_vec_tmp[CAND_MAG_SQ],
    counter[0]
  );
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  auto* RESTRICT cand_x{n_vec_tmp[CAND_X]};
  auto* RESTRICT cand_y{n_vec_tmp[CAND_Y]};
  auto* RESTRICT cand_z{n_vec_tmp[CAND_Z]};
  auto* RESTRICT cand_mag_sq{n_vec_tmp[CAND_MAG_SQ]};

  AlignedSoA<std::size_t> order{*counter[0], 1};
  std::iota(order[0], order[0] + *counter[0], std::size_t{0});

  std::sort(
    order[0],
    order[0] + *counter[0],
    [&](std::size_t a, std::size_t b) {
      return
        std::tie(cand_mag_sq[a], cand_x[a], cand_y[a], cand_z[a]) <
        std::tie(cand_mag_sq[b], cand_x[b], cand_y[b], cand_z[b]);
    }
  );

  num_unique_k_ = N / 2 + 1;
  trig_row_stride_ = this->num_unique_k();

  const real_t two_pi_inv_L{
    2.0_r * std::numbers::pi_v<real_t> / this->box_length()
  };

  const auto nv{this->n_vector()};
  const auto kv{this->k_vector()};

  dim3 assignOrbitalsThreads(256);
  dim3 assignOrbitalsBlocks(
    vmc::cudaNumBlocks(this->num_unique_k(), assignOrbitalsThreads.x)
  );
  cudaAssignOrbitals<<<assignOrbitalsBlocks, assignOrbitalsThreads>>>(
    N, this->num_unique_k(), two_pi_inv_L,
    order[0],
    cand_x, cand_y, cand_z,
    nv.x_, nv.y_, nv.z_,
    kv.x_, kv.y_, kv.z_,
    this->orbital_k_index(), this->orbital_type()
  );
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  trig_cache_ = AlignedSoA<real_t>(num_particles * this->trig_row_stride(), NUM_TRIG_ARRAYS);
  trig_scratch_ = AlignedSoA<real_t>(this->num_unique_k(), NUM_SCRATCH_TRIG);
#else
  const std::size_t N{num_orbitals()};
  const std::size_t num_particles{particles.size()};

  struct nVectorCandidate {
    int n_cand_x;
    int n_cand_y;
    int n_cand_z;
    int n_mag_sq;
  };

  const int N_MAX{static_cast<int>(vmc::ceil(vmc::cbrt(static_cast<real_t>(N)))) + 2};
  const std::size_t side{static_cast<std::size_t>((2 * N_MAX + 1))};

  std::vector<nVectorCandidate> n_candidates{};
  n_candidates.reserve(side * side * side);

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

  std::sort(
    n_candidates.begin(),
    n_candidates.end(),
    [](const nVectorCandidate& a, const nVectorCandidate& b) {
      return
        std::tie(a.n_mag_sq, a.n_cand_x, a.n_cand_y, a.n_cand_z) <
        std::tie(b.n_mag_sq, b.n_cand_x, b.n_cand_y, b.n_cand_z);
    }
  );

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

    orb_k_idx[orb_idx] = k_idx;
    orb_type[orb_idx] = 0;
    ++orb_idx;

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
#endif
}
