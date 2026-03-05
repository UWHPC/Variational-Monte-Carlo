#include "local_energy.hpp"

double EnergyTracker::kinetic_energy(const Particles& particles) const noexcept {
    const double* RESTRICT gx{particles.grad_log_psi_x_get()};
    const double* RESTRICT gy{particles.grad_log_psi_y_get()};
    const double* RESTRICT gz{particles.grad_log_psi_z_get()};
    const double* RESTRICT lap{particles.lap_log_psi_get()};

    double t_sum{};
    const std::size_t N{particles.num_particles_get()};

    for (std::size_t i{}; i<N; i++) {
      // Computes ||Grad(logPsi)||^2
      double grad_sq{gx[i]*gx[i] + gy[i]*gy[i] + gz[i]*gz[i]};
      
      // Accumulate Lapl(LogPsi) + ||Grad(LogPsi)||^2
      t_sum += (lap[i] + grad_sq);
    }

    return -0.5 * t_sum;
}

double EnergyTracker::potential_energy(const Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept {
    const double* RESTRICT px{particles.pos_x_get()};
    const double* RESTRICT py{particles.pos_y_get()};
    const double* RESTRICT pz{particles.pos_z_get()};

    double v_sum{};
    const std::size_t N{particles.num_particles_get()};

    for (std::size_t i{}; i < N ; i++) {
      for (std::size_t j{i + 1}; j < N; j++) {
        double r_ij{pbc.distance(px[i], py[i], pz[i], px[j], py[j], pz[j])};
        v_sum += 1.0 / r_ij;
      }
    }

    return v_sum;
}
double EnergyTracker::eval_total_energy(const Particles& p, const PeriodicBoundaryCondition& pbc) const noexcept {
  return kinetic_energy(p) + potential_energy(p, pbc);
}