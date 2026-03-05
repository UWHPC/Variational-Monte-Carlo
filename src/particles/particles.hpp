#pragma once

#include "../utilities/aligned_soa.hpp"
#include "../utilities/memory.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>

class Particles {
private:
    // All sub-arrays owned by Particles:
    enum ArrayIndex : std::size_t {
        POS_X_,
        POS_Y_,
        POS_Z_,
        GRAD_X_,
        GRAD_Y_,
        GRAD_Z_,
        LOG_PSI_,
        LAP_LOG_PSI_,
        NUM_SUB_ARRAYS_
    };
    // Number of particles:
    std::size_t num_particles_;

    // Aligned memory block:
    AlignedSoA<double> particle_data_;

public:
    explicit Particles(std::size_t num_particles)
        : num_particles_{num_particles}, particle_data_{num_particles, NUM_SUB_ARRAYS_} {}

    // Physical number of particles
    [[nodiscard]] std::size_t num_particles_get() const { return num_particles_; }

    // Length of the padded stride
    [[nodiscard]] std::size_t padding_stride_get() const { return particle_data_.stride(); }

    // Raw Pointers:
    // X position of particle
    [[nodiscard]] double* pos_x_get() noexcept { return particle_data_[POS_X_]; }
    [[nodiscard]] double const* pos_x_get() const noexcept { return particle_data_[POS_X_]; }

    // Y position of particle
    [[nodiscard]] double* pos_y_get() noexcept { return particle_data_[POS_Y_]; }
    [[nodiscard]] double const* pos_y_get() const noexcept { return particle_data_[POS_Y_]; }

    // Z position of particle
    [[nodiscard]] double* pos_z_get() noexcept { return particle_data_[POS_Z_]; }
    [[nodiscard]] double const* pos_z_get() const noexcept { return particle_data_[POS_Z_]; }

    // X component of gradient( log|PSI| )
    [[nodiscard]] double* grad_log_psi_x_get() noexcept { return particle_data_[GRAD_X_]; }
    [[nodiscard]] double const* grad_log_psi_x_get() const noexcept { return particle_data_[GRAD_X_]; }

    // Y component of gradient( log|PSI| )
    [[nodiscard]] double* grad_log_psi_y_get() noexcept { return particle_data_[GRAD_Y_]; }
    [[nodiscard]] double const* grad_log_psi_y_get() const noexcept { return particle_data_[GRAD_Y_]; }

    // Z component of gradient( log|PSI| )
    [[nodiscard]] double* grad_log_psi_z_get() noexcept { return particle_data_[GRAD_Z_]; }
    [[nodiscard]] double const* grad_log_psi_z_get() const noexcept { return particle_data_[GRAD_Z_]; }

    // Log|PSI|
    [[nodiscard]] double* log_psi_get() noexcept { return particle_data_[LOG_PSI_]; }
    [[nodiscard]] double const* log_psi_get() const noexcept { return particle_data_[LOG_PSI_]; }

    // Laplacian of Log|PSI|
    [[nodiscard]] double* lap_log_psi_get() noexcept { return particle_data_[LAP_LOG_PSI_]; }
    [[nodiscard]] double const* lap_log_psi_get() const noexcept { return particle_data_[LAP_LOG_PSI_]; }
};
