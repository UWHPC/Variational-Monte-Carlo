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

    // Cached structure factors: Sum(G) = sum_j exp(i G dot r_j)
    // sum_real_[g] = sum_j cos(G_g dot r_j)
    // sum_imag_[g] = sum_j sin(G_g dot r_j)
    std::vector<double> sum_real_;
    std::vector<double> sum_imag_;

public:
    explicit EnergyTracker(double box_length, double num_particles);

    // Full O(num_G * N) calculation - only call once after positions are initialized
    void initialize_structure_factors(const Particles& particles) noexcept;

    // O(num_G) calculation - called for after positions are first initialized
    void update_structure_factors(double old_x, double old_y, double old_z, double new_x, double new_y,
                                  double new_z) noexcept;

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

    [[nodiscard]] std::vector<double>& sum_real_get() noexcept { return sum_real_; }
    [[nodiscard]] std::vector<double>& sum_imag_get() noexcept { return sum_imag_; }

    [[nodiscard]] const std::vector<double>& sum_real_get() const noexcept { return sum_real_; }
    [[nodiscard]] const std::vector<double>& sum_imag_get() const noexcept { return sum_imag_; }

    [[nodiscard]] std::size_t num_g_vectors_get() const noexcept { return G_vector_x_.size(); }

    double kinetic_energy(const Particles& particles) const noexcept;
    double potential_energy(const Particles& particles, const PeriodicBoundaryCondition& pbc) const noexcept;
};