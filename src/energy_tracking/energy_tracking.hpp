#pragma once

#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"
#include "../utilities/aligned_soa.hpp"

#include <cstddef>
#include <vector>

class EnergyTracker {
private:
    // Constants for Ewald Energy:
    double ewald_alpha_;
    double ewald_correction_;
    double ewald_background_;

    // G-vectors for Ewald energy:
    std::vector<double> G_vector_x_;
    std::vector<double> G_vector_y_;
    std::vector<double> G_vector_z_;
    std::vector<double> G_vector_weights_;

public:
    explicit EnergyTracker(double box_length, double num_particles);

    double eval_total_energy(const Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept;

private:
    [[nodiscard]] double ewald_alpha_get() const noexcept { return ewald_alpha_; }
    [[nodiscard]] double ewald_correction_get() const noexcept { return ewald_correction_; }
    [[nodiscard]] double ewald_background_get() const noexcept { return ewald_background_; }

    [[nodiscard]] std::vector<double>& G_vector_x_get() noexcept { return G_vector_x_; }
    [[nodiscard]] std::vector<double>& G_vector_y_get() noexcept { return G_vector_y_; }
    [[nodiscard]] std::vector<double>& G_vector_z_get() noexcept { return G_vector_z_; }
    [[nodiscard]] std::vector<double>& G_vector_weights_get() noexcept { return G_vector_weights_; }

    [[nodiscard]] const std::vector<double>& G_vector_x_get() const noexcept { return G_vector_x_; }
    [[nodiscard]] const std::vector<double>& G_vector_y_get() const noexcept { return G_vector_y_; }
    [[nodiscard]] const std::vector<double>& G_vector_z_get() const noexcept { return G_vector_z_; }
    [[nodiscard]] const std::vector<double>& G_vector_weights_get() const noexcept { return G_vector_weights_; }

    [[nodiscard]] std::size_t num_g_vectors_get() const noexcept { return G_vector_x_.size(); }

    double kinetic_energy(const Particles& particles) const noexcept;
    double potential_energy(const Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept;
};