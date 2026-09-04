#include <xpu/xpu.hpp>
#include "energy_tracking.hpp"
#include "energy_tracking_kernels.hpp"

#include <numbers>
#include <vector>

EnergyTracker::EnergyTracker(
  fp_t box_length,
  std::size_t num_particles,
  std::size_t num_walkers
)
: box_length_{box_length}
, ewald_alpha_{6.0_fp / box_length}
, ewald_correction_{-6.0_fp * scast<fp_t>(num_particles)/(xpu::sqrt(std::numbers::pi_v<fp_t>) * box_length)}
, ewald_background_{-std::numbers::pi_v<fp_t>*scast<fp_t>(num_particles)*scast<fp_t>(num_particles)/(72.0_fp * box_length)}
, num_g_vectors_{}
, num_walkers_{num_walkers}
, shared_data_{0uz}
, walker_data_{num_walkers, 0uz}
, walker_scalars_{num_walkers}
, reduction_scratch_{num_walkers}
{
  const fp_t two_pi_over_L{2.0_fp * std::numbers::pi_v<fp_t> / box_length};
  const fp_t four_alpha_sq{4.0_fp * ewald_alpha_ * ewald_alpha_};
  const fp_t cutoff_factor{-xpu::log(EWALD_RECIPROCAL_TOLERANCE)};

  const fp_t g_max_mag_sq{four_alpha_sq * cutoff_factor};
  const int m_max{scast<int>(xpu::ceil(xpu::sqrt(g_max_mag_sq) / two_pi_over_L)) + 1};

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

        const fp_t g_cand_x{two_pi_over_L * scast<fp_t>(m_x)};
        const fp_t g_cand_y{two_pi_over_L * scast<fp_t>(m_y)};
        const fp_t g_cand_z{two_pi_over_L * scast<fp_t>(m_z)};
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

  shared_data_ = xpu::soa<fp_t, idx(SharedArray::NUM_ARRAYS)>{
    this->num_g_vectors()
  };
  walker_data_ = xpu::soa_batch<fp_t, idx(WalkerArray::NUM_ARRAYS)>{
    this->walker_count(),
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

void EnergyTracker::initialize_reciprocal_energy(std::size_t walker) noexcept {
  auto energy{this->view(walker)};
  const fp_t prefactor{
    1.0_fp / (2.0_fp * std::numbers::pi_v<fp_t>
      * energy.box_length * energy.box_length * energy.box_length)
  };
  const auto reciprocal_sum{
    kernel::energy::initialize_reciprocal_energy(energy)
  };

  *energy.reciprocal_energy = prefactor * reciprocal_sum;
}

void EnergyTracker::initialize_real_energy(
  Particles::View particles,
  std::size_t walker
) noexcept {
  auto energy{this->view(walker)};
  *energy.real_energy = kernel::energy::initialize_real_energy(energy, particles);
}

void EnergyTracker::initialize_structure_factors(
  Particles::View particles,
  std::size_t walker
) noexcept {
  kernel::energy::initialize_structure_factors(this->view(walker), particles);
}

void EnergyTracker::update_structure_factors(
  xpu::array<fp_t, idx(Axis::NUM)> old_pos,
  xpu::array<fp_t, idx(Axis::NUM)> new_pos,
  std::size_t walker
) noexcept {
  kernel::energy::update_structure_factors(this->view(walker), old_pos, new_pos);

  this->initialize_reciprocal_energy(walker);
}
fp_t EnergyTracker::kinetic_energy(
  Particles::View particles,
  std::size_t walker
) noexcept {
  return kernel::energy::kinetic_energy(
    this->view(walker),
    particles
  );
}

void EnergyTracker::update_real_energy(
  std::size_t moved,
  xpu::array<fp_t, idx(Axis::NUM)> old_pos,
  Particles::View particles,
  std::size_t walker
) noexcept {
  auto energy{this->view(walker)};
  *energy.real_energy += kernel::energy::update_real_energy(
    energy, moved, old_pos, particles
  );
}
