#include <catch2/catch_test_macros.hpp>

#include "energy_tracking/energy_tracking.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace {

void require_near_energy(double actual, double expected, double tolerance = 1e-10) {
    REQUIRE(std::abs(actual - expected) <= tolerance);
}

void copy_positions(const Particles& source, Particles& dest) {
    const std::size_t n{source.num_particles_get()};
    std::copy_n(source.pos_x_get(), n, dest.pos_x_get());
    std::copy_n(source.pos_y_get(), n, dest.pos_y_get());
    std::copy_n(source.pos_z_get(), n, dest.pos_z_get());
}

void copy_derivatives(const Particles& source, Particles& dest) {
    const std::size_t n{source.num_particles_get()};
    std::copy_n(source.grad_log_psi_x_get(), n, dest.grad_log_psi_x_get());
    std::copy_n(source.grad_log_psi_y_get(), n, dest.grad_log_psi_y_get());
    std::copy_n(source.grad_log_psi_z_get(), n, dest.grad_log_psi_z_get());
    std::copy_n(source.lap_log_psi_get(), n, dest.lap_log_psi_get());
}

} // namespace

TEST_CASE("EnergyTracker total energy changes by expected kinetic contribution", "[energy]") {
    constexpr std::size_t n{3U};
    const PeriodicBoundaryCondition pbc{8.0};
    const EnergyTracker tracker{pbc.L_get(), static_cast<double>(n)};

    Particles reference{n};
    reference.pos_x_get()[0] = 0.3;
    reference.pos_y_get()[0] = 1.5;
    reference.pos_z_get()[0] = 2.4;
    reference.pos_x_get()[1] = 3.7;
    reference.pos_y_get()[1] = 0.2;
    reference.pos_z_get()[1] = 5.1;
    reference.pos_x_get()[2] = 6.9;
    reference.pos_y_get()[2] = 7.2;
    reference.pos_z_get()[2] = 1.1;

    Particles with_derivatives{n};
    copy_positions(reference, with_derivatives);

    with_derivatives.grad_log_psi_x_get()[0] = 0.2;
    with_derivatives.grad_log_psi_y_get()[0] = -0.1;
    with_derivatives.grad_log_psi_z_get()[0] = 0.3;
    with_derivatives.lap_log_psi_get()[0] = 0.5;

    with_derivatives.grad_log_psi_x_get()[1] = -0.4;
    with_derivatives.grad_log_psi_y_get()[1] = 0.0;
    with_derivatives.grad_log_psi_z_get()[1] = 0.1;
    with_derivatives.lap_log_psi_get()[1] = -0.2;

    with_derivatives.grad_log_psi_x_get()[2] = 0.3;
    with_derivatives.grad_log_psi_y_get()[2] = -0.2;
    with_derivatives.grad_log_psi_z_get()[2] = -0.5;
    with_derivatives.lap_log_psi_get()[2] = 0.7;

    const double energy_without_derivatives{tracker.eval_total_energy(reference, pbc)};
    const double energy_with_derivatives{tracker.eval_total_energy(with_derivatives, pbc)};

    double expected_kinetic{};
    for (std::size_t i = 0; i < n; ++i) {
        const double gx{with_derivatives.grad_log_psi_x_get()[i]};
        const double gy{with_derivatives.grad_log_psi_y_get()[i]};
        const double gz{with_derivatives.grad_log_psi_z_get()[i]};
        const double lap{with_derivatives.lap_log_psi_get()[i]};
        expected_kinetic += -0.5 * (lap + gx * gx + gy * gy + gz * gz);
    }

    require_near_energy(energy_with_derivatives - energy_without_derivatives, expected_kinetic);
}

TEST_CASE("EnergyTracker is invariant under box-periodic particle translations", "[energy]") {
    constexpr std::size_t n{3U};
    const PeriodicBoundaryCondition pbc{7.5};
    const EnergyTracker tracker{pbc.L_get(), static_cast<double>(n)};

    Particles particles{n};
    particles.pos_x_get()[0] = 1.1;
    particles.pos_y_get()[0] = 2.3;
    particles.pos_z_get()[0] = 3.7;
    particles.pos_x_get()[1] = 6.4;
    particles.pos_y_get()[1] = 5.6;
    particles.pos_z_get()[1] = 0.4;
    particles.pos_x_get()[2] = 2.9;
    particles.pos_y_get()[2] = 1.2;
    particles.pos_z_get()[2] = 7.0;

    particles.grad_log_psi_x_get()[0] = 0.12;
    particles.grad_log_psi_y_get()[0] = -0.08;
    particles.grad_log_psi_z_get()[0] = 0.05;
    particles.lap_log_psi_get()[0] = 0.2;

    particles.grad_log_psi_x_get()[1] = -0.15;
    particles.grad_log_psi_y_get()[1] = 0.11;
    particles.grad_log_psi_z_get()[1] = -0.03;
    particles.lap_log_psi_get()[1] = -0.1;

    particles.grad_log_psi_x_get()[2] = 0.07;
    particles.grad_log_psi_y_get()[2] = 0.09;
    particles.grad_log_psi_z_get()[2] = -0.04;
    particles.lap_log_psi_get()[2] = 0.05;

    Particles translated{n};
    copy_positions(particles, translated);
    copy_derivatives(particles, translated);

    const double L{pbc.L_get()};
    translated.pos_x_get()[0] += L;
    translated.pos_y_get()[0] -= 2.0 * L;
    translated.pos_z_get()[1] += 3.0 * L;
    translated.pos_x_get()[2] -= L;

    const double baseline{tracker.eval_total_energy(particles, pbc)};
    const double shifted{tracker.eval_total_energy(translated, pbc)};
    require_near_energy(shifted, baseline, 1e-9);
}

TEST_CASE("EnergyTracker handles degenerate positions and permutation symmetry", "[energy]") {
    constexpr std::size_t n{2U};
    const PeriodicBoundaryCondition pbc{6.0};
    const EnergyTracker tracker{pbc.L_get(), static_cast<double>(n)};

    Particles particles{n};
    particles.pos_x_get()[0] = 2.0;
    particles.pos_y_get()[0] = 2.0;
    particles.pos_z_get()[0] = 2.0;
    particles.pos_x_get()[1] = 2.0;
    particles.pos_y_get()[1] = 2.0;
    particles.pos_z_get()[1] = 2.0;

    particles.grad_log_psi_x_get()[0] = 0.1;
    particles.grad_log_psi_y_get()[0] = 0.2;
    particles.grad_log_psi_z_get()[0] = 0.3;
    particles.lap_log_psi_get()[0] = 0.4;
    particles.grad_log_psi_x_get()[1] = -0.3;
    particles.grad_log_psi_y_get()[1] = -0.1;
    particles.grad_log_psi_z_get()[1] = 0.2;
    particles.lap_log_psi_get()[1] = -0.5;

    const double degenerate_energy{tracker.eval_total_energy(particles, pbc)};
    REQUIRE(std::isfinite(degenerate_energy));

    Particles permuted{n};
    permuted.pos_x_get()[0] = particles.pos_x_get()[1];
    permuted.pos_y_get()[0] = particles.pos_y_get()[1];
    permuted.pos_z_get()[0] = particles.pos_z_get()[1];
    permuted.pos_x_get()[1] = particles.pos_x_get()[0];
    permuted.pos_y_get()[1] = particles.pos_y_get()[0];
    permuted.pos_z_get()[1] = particles.pos_z_get()[0];

    permuted.grad_log_psi_x_get()[0] = particles.grad_log_psi_x_get()[1];
    permuted.grad_log_psi_y_get()[0] = particles.grad_log_psi_y_get()[1];
    permuted.grad_log_psi_z_get()[0] = particles.grad_log_psi_z_get()[1];
    permuted.lap_log_psi_get()[0] = particles.lap_log_psi_get()[1];
    permuted.grad_log_psi_x_get()[1] = particles.grad_log_psi_x_get()[0];
    permuted.grad_log_psi_y_get()[1] = particles.grad_log_psi_y_get()[0];
    permuted.grad_log_psi_z_get()[1] = particles.grad_log_psi_z_get()[0];
    permuted.lap_log_psi_get()[1] = particles.lap_log_psi_get()[0];

    const double permuted_energy{tracker.eval_total_energy(permuted, pbc)};
    require_near_energy(permuted_energy, degenerate_energy, 1e-10);
}
