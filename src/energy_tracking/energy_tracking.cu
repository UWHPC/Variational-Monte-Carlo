#include <xpu/xpu.hpp>
#include "energy_tracking.hpp"
#include "energy_tracking_kernels.hpp"

#include <numbers>
#include <vector>

EnergyTracker::EnergyTracker(fp_t box_length, std::size_t num_particles)
: box_length_{box_length}
, ewald_alpha_{6.0_fp / box_length}
, ewald_correction_{-6.0_fp * static_cast<fp_t>(num_particles)/(xpu::sqrt(std::numbers::pi_v<fp_t>) * box_length)}
, ewald_background_{-std::numbers::pi_v<fp_t>*static_cast<fp_t>(num_particles)*static_cast<fp_t>(num_particles)/(72.0_fp * box_length)}
, V_recip_{}
, V_real_{}
, num_g_vectors_{}
, data_{0uz}
, reduction_scratch_{1uz}
{
  const fp_t two_pi_over_L{2.0_fp * std::numbers::pi_v<fp_t> / box_length};
  const fp_t four_alpha_sq{4.0_fp * ewald_alpha_ * ewald_alpha_};
  const fp_t cutoff_factor{-xpu::log(EWALD_RECIPROCAL_TOLERANCE)};

  const fp_t g_max_mag_sq{four_alpha_sq * cutoff_factor};
  const int m_max{static_cast<int>(xpu::ceil(xpu::sqrt(g_max_mag_sq) / two_pi_over_L)) + 1};

  std::vector<fp_t> g_x{};
  std::vector<fp_t> g_y{};
  std::vector<fp_t> g_z{};
  std::vector<fp_t> weights{};

  for (int m_x = -m_max; m_x <= m_max; ++m_x) {
    for (int m_y = -m_max; m_y <= m_max; ++m_y) {
      for (int m_z = -m_max; m_z <= m_max; ++m_z) {
        if (m_x < 0) { continue; }
        if (m_x == 0 && m_y < 0) { continue; }
        if (m_x == 0 && m_y == 0 && m_z <= 0) { continue; }

        const fp_t g_cand_x{two_pi_over_L * static_cast<fp_t>(m_x)};
        const fp_t g_cand_y{two_pi_over_L * static_cast<fp_t>(m_y)};
        const fp_t g_cand_z{two_pi_over_L * static_cast<fp_t>(m_z)};
        const fp_t g_cand_mag_sq{
          g_cand_x * g_cand_x +
          g_cand_y * g_cand_y +
          g_cand_z * g_cand_z
        };

        if (g_cand_mag_sq > g_max_mag_sq) { continue; }

        g_x.emplace_back(g_cand_x);
        g_y.emplace_back(g_cand_y);
        g_z.emplace_back(g_cand_z);
        weights.emplace_back(
          8.0_fp *
          std::numbers::pi_v<fp_t> * std::numbers::pi_v<fp_t> /
          g_cand_mag_sq *
          xpu::exp(-g_cand_mag_sq / four_alpha_sq)
        );
      }
    }
  }

  num_g_vectors_ = g_x.size();

  data_ = xpu::soa<fp_t, idx(ArrayIndex::NUM_ARRAYS)>{
    this->num_g_vectors()
  };

  auto g_vector{this->g_vector()};

  xpu::copy_n(
    g_vector[idx(Axis::X)],
    g_x.data(),
    this->num_g_vectors()
  );
  xpu::copy_n(
    g_vector[idx(Axis::Y)],
    g_y.data(),
    this->num_g_vectors()
  );
  xpu::copy_n(
    g_vector[idx(Axis::Z)],
    g_z.data(),
    this->num_g_vectors()
  );
  xpu::copy_n(
    this->g_weights(),
    weights.data(),
    this->num_g_vectors()
  );
}

void EnergyTracker::initialize_reciprocal_energy() noexcept {
  const fp_t L{box_length_};
  const fp_t prefactor{
    1.0_fp / (2.0_fp * std::numbers::pi_v<fp_t> * L * L * L)
  };
  const auto reciprocal_sum{
    kernel::energy::initialize_reciprocal_energy(
      this->num_g_vectors(),
      this->g_weights(),
      this->sum_real(), this->sum_imag(),
      this->reduction_scratch()
    )
  };

  V_recip_ = prefactor * reciprocal_sum;
}

void EnergyTracker::initialize_real_energy(const Particles& particles) noexcept {
  const fp_t L{box_length_};
  const fp_t half_L{0.5_fp * L};

  V_real_ = kernel::energy::initialize_real_energy(
    L, half_L, ewald_alpha_,
    particles.pos(),
    this->reduction_scratch()
  );
}

void EnergyTracker::initialize_structure_factors(const Particles& particles) noexcept {
  kernel::energy::initialize_structure_factors(
    this->g_vector(), particles.pos(),
    this->sum_real(), this->sum_imag()
  );
}

void EnergyTracker::update_structure_factors(
  xpu::array<fp_t, idx(Axis::NUM)> old_pos,
  xpu::array<fp_t, idx(Axis::NUM)> new_pos
) noexcept {
  kernel::energy::update_structure_factors(
    old_pos, new_pos,
    this->g_vector(),
    this->sum_real(), this->sum_imag()
  );

  this->initialize_reciprocal_energy();
}
fp_t EnergyTracker::kinetic_energy(const Particles& particles) const noexcept {
  return kernel::energy::kinetic_energy(
    particles.derivatives(),
    this->reduction_scratch()
  );
}

void EnergyTracker::update_real_energy(
  std::size_t moved,
  xpu::array<fp_t, idx(Axis::NUM)> old_pos,
  const Particles& particles
) noexcept {
  const fp_t L{box_length_};
  const fp_t half_L{0.5_fp * L};

  V_real_ += kernel::energy::update_real_energy(
    moved,
    L, half_L, ewald_alpha_,
    old_pos, particles.pos(),
    this->reduction_scratch()
  );
}
