#include "local_energy.hpp"

double EnergyTracker::kinetic_energy(const Particles& particles) const noexcept {
    const double* RESTRICT grad_x{particles.grad_log_psi_x_get()};
    const double* RESTRICT grad_y{particles.grad_log_psi_y_get()};
    const double* RESTRICT grad_z{particles.grad_log_psi_z_get()};
    const double* RESTRICT lap{particles.lap_log_psi_get()};

    // Kinetic
    double T_sum{};
    const std::size_t N{particles.num_particles_get()};

    for (std::size_t i{}; i < N; i++) {
        // Computes ||Grad(logPsi)||^2
        const double grad_sq{grad_x[i] * grad_x[i] + grad_y[i] * grad_y[i] + grad_z[i] * grad_z[i]};

        // Accumulate Lapl(LogPsi) + ||Grad(LogPsi)||^2
        T_sum += (lap[i] + grad_sq);
    }

    return -0.5 * T_sum;
}

double EnergyTracker::potential_energy(const Particles& particles,
                                       const PeriodicBoundaryCondition& pbc) const noexcept {
    const double* RESTRICT px{particles.pos_x_get()};
    const double* RESTRICT py{particles.pos_y_get()};
    const double* RESTRICT pz{particles.pos_z_get()};

    // Potential
    double V_sum{};
    const std::size_t N{particles.num_particles_get()};

    for (std::size_t i{}; i < N; i++) {
        for (std::size_t j{i + 1}; j < N; j++) {
            const double r_ij{pbc.distance(px[i], py[i], pz[i], px[j], py[j], pz[j])};
            V_sum += 1.0 / r_ij;
        }
    }

    return V_sum;
}
double EnergyTracker::eval_total_energy(const Particles& particles,
                                        const PeriodicBoundaryCondition& pbc) const noexcept {
    return kinetic_energy(particles) + potential_energy(particles, pbc);
}