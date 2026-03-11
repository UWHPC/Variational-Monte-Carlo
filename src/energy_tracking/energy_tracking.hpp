#pragma once

#include "../particles/particles.hpp"
#include "../utilities/aligned_soa.hpp"

#include <cstddef>
#include <vector>

class EnergyTracker {
private:
    // Box length:
    double box_length_;

    // Constants for Ewald Energy:
    static constexpr double EWALD_RECIPROCAL_TOLERANCE{1.0e-6};
    double ewald_alpha_;
    double ewald_correction_;
    double ewald_background_;

    // Energy caching:
    double V_recip_;
    double V_real_;

    // Number of G-vectors:
    std::size_t num_g_vectors_;

    // G-vector and structure factors:
    // Cached structure factors: Sum(G) = sum_j exp(i G dot r_j)
    // sum_real_[g] = sum_j cos(G_g dot r_j)
    // sum_imag_[g] = sum_j sin(G_g dot r_j)
    enum ArrayIndex : std::size_t { G_X_, G_Y_, G_Z_, G_WEIGHTS_, S_REAL_, S_IMAG_, NUM_ARRAYS_ };
    AlignedSoA<double> data_;

public:
    explicit EnergyTracker(double box_length, double num_particles);

    // Initialize energies - allows caching them
    void initialize_reciprocal_energy() noexcept;
    void initialize_real_energy(const Particles& particles) noexcept;

    // Full O(num_G * N) calculation - only call once after positions are initialized
    void initialize_structure_factors(const Particles& particles) noexcept;

    // O(num_G) calculation - called for after positions are first initialized
    void update_structure_factors(double old_x, double old_y, double old_z, double new_x, double new_y,
                                  double new_z) noexcept;

    // Updates the real energy term
    void update_real_energy(std::size_t moved_idx, double old_x, double old_y, double old_z,
                            const Particles& particles) noexcept;

    double eval_total_energy(const Particles& particles) const noexcept;

private:
    [[nodiscard]] double ewald_alpha_get() const noexcept { return ewald_alpha_; }
    [[nodiscard]] double ewald_correction_get() const noexcept { return ewald_correction_; }
    [[nodiscard]] double ewald_background_get() const noexcept { return ewald_background_; }

    [[nodiscard]] double* G_vector_x_get() noexcept { return data_[G_X_]; }
    [[nodiscard]] double* G_vector_y_get() noexcept { return data_[G_Y_]; }
    [[nodiscard]] double* G_vector_z_get() noexcept { return data_[G_Z_]; }
    [[nodiscard]] double* G_vector_weights_get() noexcept { return data_[G_WEIGHTS_]; }

    [[nodiscard]] const double* G_vector_x_get() const noexcept { return data_[G_X_]; }
    [[nodiscard]] const double* G_vector_y_get() const noexcept { return data_[G_Y_]; }
    [[nodiscard]] const double* G_vector_z_get() const noexcept { return data_[G_Z_]; }
    [[nodiscard]] const double* G_vector_weights_get() const noexcept { return data_[G_WEIGHTS_]; }

    [[nodiscard]] double* sum_real_get() noexcept { return data_[S_REAL_]; }
    [[nodiscard]] double* sum_imag_get() noexcept { return data_[S_IMAG_]; }

    [[nodiscard]] const double* sum_real_get() const noexcept { return data_[S_REAL_]; }
    [[nodiscard]] const double* sum_imag_get() const noexcept { return data_[S_IMAG_]; }

    [[nodiscard]] std::size_t num_g_vectors_get() const noexcept { return num_g_vectors_; }
    [[nodiscard]] std::size_t& num_g_vectors_set() noexcept { return num_g_vectors_; }

    [[nodiscard]] double& V_recip_set() { return V_recip_; }
    [[nodiscard]] double V_recip_get() const { return V_recip_; }

    [[nodiscard]] double& V_real_set() { return V_real_; }
    [[nodiscard]] double V_real_get() const { return V_real_; }

    double kinetic_energy(const Particles& particles) const noexcept;
    double potential_energy() const noexcept;
};