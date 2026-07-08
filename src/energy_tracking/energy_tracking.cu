#include "energy_tracking.cuh"

#include <numbers>

EnergyTracker::EnergyTracker(real_t box_length, real_t num_particles)
: box_length_{box_length}
, ewald_alpha_{6.0_r / box_length}
, ewald_correction_{-6.0_r * num_particles / (vmc::sqrt(std::numbers::pi_v<real_t>) * box_length)}
, ewald_background_{-std::numbers::pi_v<real_t> * num_particles * num_particles / (72.0_r * box_length)}
, V_recip_{}
, V_real_{}
, num_g_vectors_{}
, data_{} {
  const real_t two_pi_over_L{2.0_r * std::numbers::pi_v<real_t> / box_length};
  const real_t four_alpha_sq{4.0_r * ewald_alpha_ * ewald_alpha_};
  const real_t cutoff_factor{-vmc::log(EWALD_RECIPROCAL_TOLERANCE)};

  const real_t g_max_mag_sq{four_alpha_sq * cutoff_factor};
  const int m_max{static_cast<int>(vmc::ceil(vmc::sqrt(g_max_mag_sq) / two_pi_over_L)) + 1};

  std::vector<real_t> tmp_x, tmp_y, tmp_z, tmp_w;

  auto& g_x{tmp_x};
  auto& g_y{tmp_y};
  auto& g_z{tmp_z};
  auto& weights{tmp_w};

  // G = 2*pi / L
  for (int m_x = -m_max; m_x <= m_max; ++m_x) {
    for (int m_y = -m_max; m_y <= m_max; ++m_y) {
      for (int m_z = -m_max; m_z <= m_max; ++m_z) {
        // Only keeps "positive half-sphere"
        // Make up for this by 2x the weights
        if (m_x < 0)
          continue;
        if (m_x == 0 && m_y < 0)
          continue;
        if (m_x == 0 && m_y == 0 && m_z <= 0)
          continue;

        const real_t g_cand_x{two_pi_over_L * static_cast<real_t>(m_x)};
        const real_t g_cand_y{two_pi_over_L * static_cast<real_t>(m_y)};
        const real_t g_cand_z{two_pi_over_L * static_cast<real_t>(m_z)};
        const real_t g_cand_mag_sq{
          g_cand_x * g_cand_x +
          g_cand_y * g_cand_y +
          g_cand_z * g_cand_z
        };

        if (g_cand_mag_sq > g_max_mag_sq)
          continue;

        g_x.emplace_back(g_cand_x);
        g_y.emplace_back(g_cand_y);
        g_z.emplace_back(g_cand_z);
        weights.emplace_back(
          8.0_r * std::numbers::pi_v<real_t> * std::numbers::pi_v<real_t> / g_cand_mag_sq *
          vmc::exp(-g_cand_mag_sq / four_alpha_sq)
        );
      }
    }
  }

  num_g_vectors_ = g_x.size();

  data_ = AlignedSoA<real_t>(num_g_vectors_, NUM_ARRAYS);

  const auto g_dst{g_vector()};
  std::copy_n(tmp_x.data(), num_g_vectors_, g_dst.x_);
  std::copy_n(tmp_y.data(), num_g_vectors_, g_dst.y_);
  std::copy_n(tmp_z.data(), num_g_vectors_, g_dst.z_);
  std::copy_n(tmp_w.data(), num_g_vectors_, g_weights());
}

void EnergyTracker::initialize_reciprocal_energy() noexcept {
  const real_t L{box_length_};
  const real_t prefactor{1.0_r / (2.0_r * std::numbers::pi_v<real_t> * L * L * L)};

  const std::size_t num_G{num_g_vectors_};
  const real_t* RESTRICT g_weights{this->g_weights()};
  const real_t* RESTRICT sum_real{this->sum_real()};
  const real_t* RESTRICT sum_imag{this->sum_imag()};

  ASSUME_ALIGNED(g_weights, SIMD_BYTES);
  ASSUME_ALIGNED(sum_real, SIMD_BYTES);
  ASSUME_ALIGNED(sum_imag, SIMD_BYTES);

  real_t sum{};
  #pragma omp simd reduction(+ : sum)
  for (std::size_t g = 0; g < num_G; ++g) {
    sum += g_weights[g] * (sum_real[g] * sum_real[g] + sum_imag[g] * sum_imag[g]);
  }
  V_recip_ = prefactor * sum;
}

void EnergyTracker::initialize_real_energy(const Particles& particles) noexcept {
  const std::size_t N{particles.size()};
  const real_t L{box_length_};
  const real_t neg_L{-1.0_r * L};
  const real_t half_L{0.5_r * L};
  const real_t neg_half_L{-1.0_r * half_L};
  const real_t alpha{ewald_alpha_};

  const auto pos{particles.pos().align()};

  real_t sum{};
  for (std::size_t i = 0; i < N; ++i) {
    real_t local_sum{};

    // Not vectorzied: loop contains a mathematical function
    #pragma omp simd reduction(+ : local_sum)
    for (std::size_t j = i + 1; j < N; ++j) {
      real_t dx{pos.x_[i] - pos.x_[j]};
      real_t dy{pos.y_[i] - pos.y_[j]};
      real_t dz{pos.z_[i] - pos.z_[j]};

      dx += L * (dx <= neg_half_L) + neg_L * (dx > half_L);
      dy += L * (dy <= neg_half_L) + neg_L * (dy > half_L);
      dz += L * (dz <= neg_half_L) + neg_L * (dz > half_L);

      const real_t r{
        vmc::sqrt(
          dx * dx +
          dy * dy +
          dz * dz
        )
      };
      const real_t inv_r{(r < 1e-12_r) ? 1.0_r : 1.0_r / r};
      local_sum += vmc::erfc(alpha * r) * inv_r;
    }
    sum += local_sum;
  }
  V_real_ = sum;
}

void EnergyTracker::initialize_structure_factors(const Particles& particles) noexcept {
  const std::size_t N{particles.size()};
  const std::size_t num_G{num_g_vectors_};

  const auto pos{particles.pos().align()};
  const auto gv{g_vector().align()};

  real_t* RESTRICT sum_real{this->sum_real()};
  real_t* RESTRICT sum_imag{this->sum_imag()};

  ASSUME_ALIGNED(sum_real, SIMD_BYTES);
  ASSUME_ALIGNED(sum_imag, SIMD_BYTES);

  for (std::size_t g = 0; g < num_G; ++g) {
    real_t cos_sum{};
    real_t sin_sum{};

    // Not vectorzied: loop contains a mathematical function
    #pragma omp simd reduction(+ : cos_sum, sin_sum)
    for (std::size_t j = 0; j < N; ++j) {
      const real_t G_dot_r{
        gv.x_[g] * pos.x_[j] +
        gv.y_[g] * pos.y_[j] +
        gv.z_[g] * pos.z_[j]
      };
      real_t cos_temp{};
      real_t sin_temp{};

      vmc::sincos(G_dot_r, &sin_temp, &cos_temp);

      cos_sum += cos_temp;
      sin_sum += sin_temp;
    }

    sum_real[g] = cos_sum;
    sum_imag[g] = sin_sum;
  }
}

namespace {

__inline__ __device__ double warp_reduce_sum(double val) {
    for (int offset = 16; offset > 0; offset /= 2) {
        val += __shfl_down_sync(0xffffffff, val, offset);
    }
    return val;
}

__inline__ __device__ double block_reduce_sum(double val) {
    __shared__ double shared[32];

    int lane = threadIdx.x % 32;
    int warp_id = threadIdx.x / 32;

    val = warp_reduce_sum(val);

    if (lane == 0) {
        shared[warp_id] = val;
    }
    __syncthreads();

    int num_warps = (blockDim.x + 31) / 32;
    val = (threadIdx.x < num_warps) ? shared[lane] : 0.0;

    if (warp_id == 0) {
        val = warp_reduce_sum(val);
    }
    return val;
}

constexpr int kBlockSize{256};

[[nodiscard]] int grid_size_for(std::size_t n, int block_size) {
    return static_cast<int>((n + static_cast<std::size_t>(block_size) - 1) /
                            static_cast<std::size_t>(block_size));
}
// Kernel: incremental reciprocal-space structure-factor + energy update.
// Mirrors EnergyTracker::update_structure_factors on the CPU.

__global__ void update_structure_factors_kernel(
    const double* g_x, const double* g_y, const double* g_z, const double* g_weights,
    double* sum_real, double* sum_imag,
    std::size_t num_g,
    double old_x, double old_y, double old_z,
    double new_x, double new_y, double new_z,
    double prefactor,
    double* d_v_recip)
{
    std::size_t g = blockIdx.x * blockDim.x + threadIdx.x;

    double local_delta = 0.0;

    if (g < num_g) {
        double old_dot = g_x[g] * old_x + g_y[g] * old_y + g_z[g] * old_z;
        double new_dot = g_x[g] * new_x + g_y[g] * new_y + g_z[g] * new_z;

        double old_sin, old_cos, new_sin, new_cos;
        sincos(old_dot, &old_sin, &old_cos);
        sincos(new_dot, &new_sin, &new_cos);

        double dr = new_cos - old_cos;
        double di = new_sin - old_sin;

        double sr = sum_real[g];
        double si = sum_imag[g];

        local_delta = g_weights[g] * (2.0 * (sr * dr + si * di) + dr * dr + di * di);

        sum_real[g] = sr + dr;
        sum_imag[g] = si + di;
    }

    double block_sum = block_reduce_sum(local_delta);

    if (threadIdx.x == 0) {
        atomicAdd(d_v_recip, prefactor * block_sum);
    }
}

} // namespace

// Device memory lifecycle - called from the constructor / destructor

void EnergyTracker::init_device_gpu() {
    const std::size_t num_g{num_g_vectors_get()};
    const std::size_t bytes{num_g * sizeof(double)};

    cudaMalloc(&d_g_x_, bytes);
    cudaMalloc(&d_g_y_, bytes);
    cudaMalloc(&d_g_z_, bytes);
    cudaMalloc(&d_g_weights_, bytes);
    cudaMalloc(&d_sum_real_, bytes);
    cudaMalloc(&d_sum_imag_, bytes);
    cudaMalloc(&d_v_recip_, sizeof(double));

    // G-vectors are static after construction - copy once.
    cudaMemcpy(d_g_x_, G_vector_x_get(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_g_y_, G_vector_y_get(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_g_z_, G_vector_z_get(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_g_weights_, G_vector_weights_get(), bytes, cudaMemcpyHostToDevice);

    // sum_real/sum_imag start at zero, same as the host AlignedSoA does.
    cudaMemset(d_sum_real_, 0, bytes);
    cudaMemset(d_sum_imag_, 0, bytes);

    const double zero{0.0};
    cudaMemcpy(d_v_recip_, &zero, sizeof(double), cudaMemcpyHostToDevice);
}

void EnergyTracker::free_device_gpu() noexcept {
    cudaFree(d_g_x_);
    cudaFree(d_g_y_);
    cudaFree(d_g_z_);
    cudaFree(d_g_weights_);
    cudaFree(d_sum_real_);
    cudaFree(d_sum_imag_);
    cudaFree(d_v_recip_);

    d_g_x_ = nullptr;
    d_g_y_ = nullptr;
    d_g_z_ = nullptr;
    d_g_weights_ = nullptr;
    d_sum_real_ = nullptr;
    d_sum_imag_ = nullptr;
    d_v_recip_ = nullptr;
}

// Public GPU method - launches the kernel above
void EnergyTracker::update_structure_factors_gpu(double old_x, double old_y, double old_z,
                                                 double new_x, double new_y, double new_z,
                                                 cudaStream_t stream) noexcept {
    const std::size_t num_g{num_g_vectors_get()};
    if (num_g == 0U) {
        return;
    }

    const double L{box_length_get()};
    const double prefactor{1.0 / (2.0 * std::numbers::pi * L * L * L)};

    const int grid{grid_size_for(num_g, kBlockSize)};

    update_structure_factors_kernel<<<grid, kBlockSize, 0, stream>>>(
        d_g_x_, d_g_y_, d_g_z_, d_g_weights_, d_sum_real_, d_sum_imag_, num_g,
        old_x, old_y, old_z, new_x, new_y, new_z, prefactor, d_v_recip_);
}

real_t EnergyTracker::kinetic_energy(const Particles& particles) const noexcept {
  const auto grad{particles.grad_log_psi().align()};
  const real_t* RESTRICT lap{particles.lap_log_psi()};

  ASSUME_ALIGNED(lap, SIMD_BYTES);

  // Kinetic
  real_t T_sum{};
  const std::size_t N{particles.size()};

  #pragma omp simd reduction(+ : T_sum)
  for (std::size_t i = 0; i < N; ++i) {
    // Computes ||Grad(logPsi)||^2
    const real_t grad_sq{
      grad.x_[i] * grad.x_[i] +
      grad.y_[i] * grad.y_[i] +
      grad.z_[i] * grad.z_[i]
    };

    // Accumulate Lapl(LogPsi) + ||Grad(LogPsi)||^2
    T_sum += (lap[i] + grad_sq);
  }

  return -0.5_r * T_sum;
}

__global__ void update_real_energy_kernel(
    const double* p_x, const double* p_y, const double* p_z,
    std::size_t N, std::size_t moved_idx,
    double old_x, double old_y, double old_z,
    double L, double alpha,
    double* d_v_real)
{
    std::size_t j = blockIdx.x * blockDim.x + threadIdx.x;

    const double neg_L{-L};
    const double half_L{0.5 * L};
    const double neg_half_L{-half_L};

    // Must match the CPU constants exactly.
    const double p{0.3275911};
    const double a1{0.254829592};
    const double a2{-0.284496736};
    const double a3{1.421413741};
    const double a4{-1.453152027};
    const double a5{1.061405429};

    double local_delta = 0.0;

    if (j < N) {
        const double valid_mask{(j == moved_idx) ? 0.0 : 1.0};

        // New position of the moved particle - read from the (already
        // updated) position array, same convention as the CPU function.
        const double new_x{p_x[moved_idx]};
        const double new_y{p_y[moved_idx]};
        const double new_z{p_z[moved_idx]};

        // Old pair
        double dx_old{old_x - p_x[j]};
        double dy_old{old_y - p_y[j]};
        double dz_old{old_z - p_z[j]};

        dx_old += L * (dx_old <= neg_half_L) + neg_L * (dx_old > half_L);
        dy_old += L * (dy_old <= neg_half_L) + neg_L * (dy_old > half_L);
        dz_old += L * (dz_old <= neg_half_L) + neg_L * (dz_old > half_L);

        const double r_old{sqrt(dx_old * dx_old + dy_old * dy_old + dz_old * dz_old)};
        const double inv_r_old{(r_old < 1e-12) ? 1.0 : 1.0 / r_old};

        const double erfc_arg_old{alpha * r_old};
        const double t_old{1.0 / (1.0 + p * erfc_arg_old)};
        const double tau_old{t_old * (a1 + t_old * (a2 + t_old * (a3 + t_old * (a4 + t_old * a5))))};
        const double erfc_old{tau_old * exp(-erfc_arg_old * erfc_arg_old) * inv_r_old};

        // New pair
        double dx_new{new_x - p_x[j]};
        double dy_new{new_y - p_y[j]};
        double dz_new{new_z - p_z[j]};

        dx_new += L * (dx_new <= neg_half_L) + neg_L * (dx_new > half_L);
        dy_new += L * (dy_new <= neg_half_L) + neg_L * (dy_new > half_L);
        dz_new += L * (dz_new <= neg_half_L) + neg_L * (dz_new > half_L);

        const double r_new{sqrt(dx_new * dx_new + dy_new * dy_new + dz_new * dz_new)};
        const double inv_r_new{(r_new < 1e-12) ? 1.0 : 1.0 / r_new};

        const double erfc_arg_new{alpha * r_new};
        const double t_new{1.0 / (1.0 + p * erfc_arg_new)};
        const double tau_new{t_new * (a1 + t_new * (a2 + t_new * (a3 + t_new * (a4 + t_new * a5))))};
        const double erfc_new{tau_new * exp(-erfc_arg_new * erfc_arg_new) * inv_r_new};

        local_delta = valid_mask * (erfc_new - erfc_old);
    }

    double block_sum = block_reduce_sum(local_delta);

    if (threadIdx.x == 0) {
        atomicAdd(d_v_real, block_sum);
    }
}

real_t EnergyTracker::potential_energy() const noexcept {
  // Ewald constants:
  const real_t ewald_self_correction_term{ewald_correction_};
  const real_t ewald_background{ewald_background_};

  // Potentials calcualted in cache:
  const real_t V_recip{V_recip_};
  const real_t V_real{V_real_};

  // Self + background:
  return V_real + V_recip + ewald_self_correction_term + ewald_background;
}

real_t EnergyTracker::eval_total_energy(const Particles& particles) const noexcept {
  return kinetic_energy(particles) + potential_energy();
}
