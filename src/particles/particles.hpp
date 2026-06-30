#pragma once

#include "../utilities/aligned_soa.hpp"
#include "../utilities/macros.hpp"
#include "../utilities/ptr3d.hpp"

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

    [[nodiscard]]
    Ptr3D<const double> pos() const noexcept {
        return {
            particle_data_[POS_X_],
            particle_data_[POS_Y_],
            particle_data_[POS_Z_]
        };
    }

    [[nodiscard]]
    Ptr3D<double> pos() noexcept {
        return {
            particle_data_[POS_X_],
            particle_data_[POS_Y_],
            particle_data_[POS_Z_]
        };
    }

    // Gradient( log|PSI| )
    [[nodiscard]]
    Ptr3D<const double> grad_log_psi() const noexcept {
        return {
            particle_data_[GRAD_X_],
            particle_data_[GRAD_Y_],
            particle_data_[GRAD_Z_]
        };
    }

    [[nodiscard]]
    Ptr3D<double> grad_log_psi() noexcept {
        return {
            particle_data_[GRAD_X_],
            particle_data_[GRAD_Y_],
            particle_data_[GRAD_Z_]
        };
    }

    // Laplacian of Log|PSI|
    [[nodiscard]] double* lap_log_psi_get() noexcept { return particle_data_[LAP_LOG_PSI_]; }
    [[nodiscard]] double const* lap_log_psi_get() const noexcept {
        return particle_data_[LAP_LOG_PSI_];
    }
};
