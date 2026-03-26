#include "test_utilities.hpp"

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include "energy_tracking/energy_tracking.hpp"
#include "jastrow_pade/jastrow_pade.hpp"
#include "particles/particles.hpp"
#include "slater_plane_wave/slater_plane_wave.hpp"
#include "wavefunction/wavefunction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <random>
#include <vector>

namespace {

constexpr double EWALD_RECIPROCAL_TOLERANCE{1.0e-6};

double phy_determinant3x3(const double* matrix) {
    return matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7]) -
           matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6]) +
           matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6]);
}

void phy_copy_positions(const Particles& source, Particles& dest) {
    const std::size_t n{source.num_particles_get()};

    for (std::size_t i = 0; i < n; ++i) {
        dest.pos_x_get()[i] = source.pos_x_get()[i];
        dest.pos_y_get()[i] = source.pos_y_get()[i];
        dest.pos_z_get()[i] = source.pos_z_get()[i];
    }
}

double exact_real_potential(const Particles& particles, double box_length) {
    const std::size_t n{particles.num_particles_get()};
    const double alpha{6.0 / box_length};

    double total{};

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double dx{
                minimum_image(particles.pos_x_get()[i] - particles.pos_x_get()[j], box_length)};
            const double dy{
                minimum_image(particles.pos_y_get()[i] - particles.pos_y_get()[j], box_length)};
            const double dz{
                minimum_image(particles.pos_z_get()[i] - particles.pos_z_get()[j], box_length)};

            const double r{std::sqrt(dx * dx + dy * dy + dz * dz)};
            total += std::erfc(alpha * r) / r;
        }
    }

    return total;
}

double exact_reciprocal_potential(const Particles& particles, double box_length) {
    const std::size_t n{particles.num_particles_get()};
    const double alpha{6.0 / box_length};
    const double two_pi_over_l{2.0 * std::numbers::pi / box_length};
    const double four_alpha_sq{4.0 * alpha * alpha};
    const double cutoff_factor{-std::log(EWALD_RECIPROCAL_TOLERANCE)};
    const double g_max_mag_sq{four_alpha_sq * cutoff_factor};
    const int m_max{
        static_cast<int>(std::ceil(std::sqrt(g_max_mag_sq) / two_pi_over_l)) + 1};

    double weighted_sum{};

    for (int mx = -m_max; mx <= m_max; ++mx) {
        for (int my = -m_max; my <= m_max; ++my) {
            for (int mz = -m_max; mz <= m_max; ++mz) {
                if (mx < 0) {
                    continue;
                }
                if (mx == 0 && my < 0) {
                    continue;
                }
                if (mx == 0 && my == 0 && mz <= 0) {
                    continue;
                }

                const double gx{two_pi_over_l * static_cast<double>(mx)};
                const double gy{two_pi_over_l * static_cast<double>(my)};
                const double gz{two_pi_over_l * static_cast<double>(mz)};
                const double g_sq{gx * gx + gy * gy + gz * gz};

                if (g_sq > g_max_mag_sq) {
                    continue;
                }

                const double weight{
                    8.0 * std::numbers::pi * std::numbers::pi * std::exp(-g_sq / four_alpha_sq) /
                    g_sq};

                double s_real{};
                double s_imag{};

                for (std::size_t i = 0; i < n; ++i) {
                    const double phase{gx * particles.pos_x_get()[i] +
                                       gy * particles.pos_y_get()[i] +
                                       gz * particles.pos_z_get()[i]};
                    s_real += std::cos(phase);
                    s_imag += std::sin(phase);
                }

                weighted_sum += weight * (s_real * s_real + s_imag * s_imag);
            }
        }
    }

    return weighted_sum /
           (2.0 * std::numbers::pi * box_length * box_length * box_length);
}

double exact_total_potential(const Particles& particles, double box_length) {
    const double n{static_cast<double>(particles.num_particles_get())};
    const double self_correction{
        -6.0 * n / (std::sqrt(std::numbers::pi) * box_length)};
    const double background{
        -std::numbers::pi * n * n / (72.0 * box_length)};

    return exact_real_potential(particles, box_length) +
           exact_reciprocal_potential(particles, box_length) +
           self_correction + background;
}

} // namespace

TEST_CASE("Fully polarized Jastrow cusp matches the same-spin analytical form near contact",
          "[rigorous][jastrow]") {
    constexpr double box_length{50.0};
    const JastrowPade jastrow{box_length, 0.25, 1.0};
    Particles particles{2U};

    particles.pos_x_get()[0] = 0.0;
    particles.pos_y_get()[0] = 0.0;
    particles.pos_z_get()[0] = 0.0;

    for (const double r : {1.0e-5, 2.0e-5, 5.0e-5, 1.0e-4}) {
        particles.pos_x_get()[1] = r;
        particles.pos_y_get()[1] = 0.0;
        particles.pos_z_get()[1] = 0.0;

        const double expected_value{0.25 * r / (1.0 + r)};
        const double actual_value{jastrow.value(particles)};

        INFO("Checking Jastrow value against the exact same-spin Padé form near coalescence.");
        CAPTURE(r, actual_value, expected_value);
        REQUIRE(std::abs(actual_value - expected_value) <= 1e-14);

        std::vector<double> grad_x(particles.padding_stride_get(), 0.0);
        std::vector<double> grad_y(particles.padding_stride_get(), 0.0);
        std::vector<double> grad_z(particles.padding_stride_get(), 0.0);
        std::vector<double> lap(particles.padding_stride_get(), 0.0);

        jastrow.add_derivatives(particles, grad_x.data(), grad_y.data(), grad_z.data(), lap.data());

        const double first_derivative{0.25 / ((1.0 + r) * (1.0 + r))};
        const double second_derivative{-0.5 / ((1.0 + r) * (1.0 + r) * (1.0 + r))};
        const double expected_grad_particle_0{-first_derivative};
        const double expected_grad_particle_1{first_derivative};
        const double expected_laplacian{
            second_derivative + 2.0 * first_derivative / r};

        INFO("Checking same-spin Jastrow gradient and Laplacian near contact.");
        CAPTURE(r, first_derivative, second_derivative, expected_laplacian);
        REQUIRE(std::abs(grad_x[0] - expected_grad_particle_0) <= 1e-10);
        REQUIRE(std::abs(grad_x[1] - expected_grad_particle_1) <= 1e-10);
        REQUIRE(std::abs(grad_y[0]) <= 1e-14);
        REQUIRE(std::abs(grad_y[1]) <= 1e-14);
        REQUIRE(std::abs(grad_z[0]) <= 1e-14);
        REQUIRE(std::abs(grad_z[1]) <= 1e-14);
        REQUIRE(std::abs(lap[0] - expected_laplacian) <= 1e-7);
        REQUIRE(std::abs(lap[1] - expected_laplacian) <= 1e-7);
    }
}

TEST_CASE("Slater determinant changes sign under particle exchange while log_abs_det stays invariant",
          "[rigorous][slater]") {
    constexpr std::size_t n{3U};
    constexpr double box_length{9.0};

    Particles original{n};
    original.pos_x_get()[0] = 0.7;
    original.pos_y_get()[0] = 1.4;
    original.pos_z_get()[0] = 2.1;

    original.pos_x_get()[1] = 2.8;
    original.pos_y_get()[1] = 0.9;
    original.pos_z_get()[1] = 1.6;

    original.pos_x_get()[2] = 4.3;
    original.pos_y_get()[2] = 3.2;
    original.pos_z_get()[2] = 0.5;

    SlaterPlaneWave slater_original{original, box_length};
    const double log_original{slater_original.log_abs_det(original)};
    REQUIRE(std::isfinite(log_original));

    const double det_original{phy_determinant3x3(slater_original.determinant_get())};

    Particles swapped{n};
    phy_copy_positions(original, swapped);
    std::swap(swapped.pos_x_get()[0], swapped.pos_x_get()[1]);
    std::swap(swapped.pos_y_get()[0], swapped.pos_y_get()[1]);
    std::swap(swapped.pos_z_get()[0], swapped.pos_z_get()[1]);

    SlaterPlaneWave slater_swapped{swapped, box_length};
    const double log_swapped{slater_swapped.log_abs_det(swapped)};
    REQUIRE(std::isfinite(log_swapped));

    const double det_swapped{phy_determinant3x3(slater_swapped.determinant_get())};

    INFO("Swapping two particles must negate the Slater determinant but preserve log|det|.");
    CAPTURE(det_original, det_swapped, log_original, log_swapped);
    REQUIRE(std::abs(det_swapped + det_original) <= 1e-10);
    REQUIRE(std::abs(log_swapped - log_original) <= 1e-12);
}

TEST_CASE("SlaterPlaneWave is periodic under box translations of a single particle",
          "[rigorous][slater]") {
    constexpr std::size_t n{3U};
    constexpr double box_length{10.0};

    Particles particles{n};
    particles.pos_x_get()[0] = 1.1;
    particles.pos_y_get()[0] = 2.2;
    particles.pos_z_get()[0] = 3.3;

    particles.pos_x_get()[1] = 4.4;
    particles.pos_y_get()[1] = 5.5;
    particles.pos_z_get()[1] = 6.6;

    particles.pos_x_get()[2] = 7.7;
    particles.pos_y_get()[2] = 1.8;
    particles.pos_z_get()[2] = 2.9;

    SlaterPlaneWave baseline{particles, box_length};
    const double baseline_log_det{baseline.log_abs_det(particles)};
    REQUIRE(std::isfinite(baseline_log_det));

    Particles shifted{n};
    phy_copy_positions(particles, shifted);
    shifted.pos_x_get()[1] += box_length;
    shifted.pos_y_get()[1] -= box_length;
    shifted.pos_z_get()[1] += 2.0 * box_length;

    SlaterPlaneWave translated{shifted, box_length};
    const double translated_log_det{translated.log_abs_det(shifted)};
    REQUIRE(std::isfinite(translated_log_det));

    INFO("Plane-wave Slater determinant must be invariant under integer box translations.");
    CAPTURE(baseline_log_det, translated_log_det);
    REQUIRE(std::abs(translated_log_det - baseline_log_det) <= 1e-12);

    for (std::size_t idx = 0; idx < n * n; ++idx) {
        CAPTURE(idx);
        REQUIRE(std::abs(baseline.determinant_get()[idx] - translated.determinant_get()[idx]) <=
                1e-12);
    }
}

TEST_CASE("Randomized determinant_ratio matches exact determinant ratio over many configurations",
          "[rigorous][slater]") {
    constexpr std::size_t n{3U};
    constexpr double box_length{11.0};
    constexpr std::size_t samples{40U};

    std::mt19937_64 rng{123456789ULL};
    std::uniform_real_distribution<double> position_dist{0.0, box_length};
    std::uniform_real_distribution<double> move_dist{-0.20, 0.20};

    for (std::size_t sample = 0; sample < samples; ++sample) {
        Particles particles{n};

        for (std::size_t i = 0; i < n; ++i) {
            particles.pos_x_get()[i] = position_dist(rng);
            particles.pos_y_get()[i] = position_dist(rng);
            particles.pos_z_get()[i] = position_dist(rng);
        }

        SlaterPlaneWave slater{particles, box_length};
        const double baseline_log_det{slater.log_abs_det(particles)};

        INFO("Baseline Slater rebuild must be finite for randomized determinant-ratio testing.");
        CAPTURE(sample, baseline_log_det);
        REQUIRE(std::isfinite(baseline_log_det));

        const double det_old{phy_determinant3x3(slater.determinant_get())};
        REQUIRE(std::abs(det_old) > 1e-12);

        const std::size_t moved{sample % n};
        particles.pos_x_get()[moved] =
            wrap_coordinate(particles.pos_x_get()[moved] + move_dist(rng), box_length);
        particles.pos_y_get()[moved] =
            wrap_coordinate(particles.pos_y_get()[moved] + move_dist(rng), box_length);
        particles.pos_z_get()[moved] =
            wrap_coordinate(particles.pos_z_get()[moved] + move_dist(rng), box_length);

        slater.update_trig_cache(moved, particles);
        const double* const new_row{slater.build_row(moved)};
        const double ratio{slater.determinant_ratio(moved, new_row)};

        SlaterPlaneWave rebuilt{particles, box_length};
        const double rebuilt_log_det{rebuilt.log_abs_det(particles)};

        INFO("Fresh rebuild must be finite for randomized determinant-ratio testing.");
        CAPTURE(sample, moved, rebuilt_log_det);
        REQUIRE(std::isfinite(rebuilt_log_det));

        const double det_new{phy_determinant3x3(rebuilt.determinant_get())};
        const double exact_ratio{det_new / det_old};

        INFO("Fast determinant ratio must agree with exact determinant ratio over many random configurations.");
        CAPTURE(sample, moved, ratio, exact_ratio, det_old, det_new);
        REQUIRE(std::abs(ratio - exact_ratio) <= 1e-9);
    }
}

TEST_CASE("SlaterPlaneWave reports a singular determinant for duplicate particle rows",
          "[rigorous][slater]") {
    constexpr std::size_t n{3U};
    constexpr double box_length{10.0};

    Particles particles{n};
    particles.pos_x_get()[0] = 1.25;
    particles.pos_y_get()[0] = 2.50;
    particles.pos_z_get()[0] = 3.75;

    particles.pos_x_get()[1] = particles.pos_x_get()[0];
    particles.pos_y_get()[1] = particles.pos_y_get()[0];
    particles.pos_z_get()[1] = particles.pos_z_get()[0];

    particles.pos_x_get()[2] = 4.10;
    particles.pos_y_get()[2] = 0.90;
    particles.pos_z_get()[2] = 1.70;

    SlaterPlaneWave slater{particles, box_length};
    const double log_abs_det{slater.log_abs_det(particles)};
    const double det{phy_determinant3x3(slater.determinant_get())};

    INFO("Duplicate particle positions should create duplicate Slater rows and a singular matrix.");
    CAPTURE(log_abs_det, det);
    REQUIRE_FALSE(std::isfinite(log_abs_det));
    REQUIRE(std::abs(det) <= 1e-12);
}

TEST_CASE("Repeated accepted Slater updates preserve predictive determinant ratios against fresh rebuilds",
          "[rigorous][slater]") {
    constexpr std::size_t n{7U};
    constexpr double box_length{10.0};
    constexpr std::size_t steps{64U};

    Particles particles{n};
    set_stable_closed_shell_positions(particles);

    SlaterPlaneWave maintained{particles, box_length};
    const double initial_log_det{maintained.log_abs_det(particles)};

    INFO("Initial Slater matrix for the multi-step drift test must be well-defined.");
    CAPTURE(initial_log_det);
    REQUIRE(std::isfinite(initial_log_det));

    for (std::size_t step = 0; step < steps; ++step) {
        const std::size_t moved{step % n};

        const double dx{0.003 * static_cast<double>((step % 5U) + 1U)};
        const double dy{-0.0025 * static_cast<double>((step % 7U) + 1U)};
        const double dz{0.002 * static_cast<double>((step % 3U) + 1U)};

        particles.pos_x_get()[moved] =
            wrap_coordinate(particles.pos_x_get()[moved] + dx, box_length);
        particles.pos_y_get()[moved] =
            wrap_coordinate(particles.pos_y_get()[moved] + dy, box_length);
        particles.pos_z_get()[moved] =
            wrap_coordinate(particles.pos_z_get()[moved] + dz, box_length);

        maintained.update_trig_cache(moved, particles);
        const double* const accepted_row{maintained.build_row(moved)};
        const double accepted_ratio{maintained.determinant_ratio(moved, accepted_row)};

        INFO("Accepted-step determinant ratio became non-finite or too close to singular.");
        CAPTURE(step, moved, accepted_ratio);
        REQUIRE(std::isfinite(accepted_ratio));
        REQUIRE(std::abs(accepted_ratio) > 1e-10);

        maintained.accept_move(moved, accepted_row, accepted_ratio);

        SlaterPlaneWave rebuilt{particles, box_length};
        const double rebuilt_log_det{rebuilt.log_abs_det(particles)};

        INFO("Fresh rebuild produced a non-finite log determinant during the multi-step predictive-consistency test.");
        CAPTURE(step, moved, rebuilt_log_det);
        REQUIRE(std::isfinite(rebuilt_log_det));

        const double maintained_residual{slater_identity_residual(maintained)};
        const double rebuilt_residual{slater_identity_residual(rebuilt)};

        INFO("Both maintained and rebuilt Slater states must remain valid inverses of their determinant matrices.");
        CAPTURE(step, moved, maintained_residual, rebuilt_residual);
        REQUIRE(maintained_residual <= 1e-8);
        REQUIRE(rebuilt_residual <= 1e-10);

        double max_det_diff{};
        for (std::size_t idx = 0; idx < n * n; ++idx) {
            const double det_diff{
                std::abs(maintained.determinant_get()[idx] - rebuilt.determinant_get()[idx])};
            max_det_diff = std::max(max_det_diff, det_diff);
        }

        INFO("Maintained and rebuilt determinant matrices should match after each accepted update.");
        CAPTURE(step, moved, max_det_diff);
        REQUIRE(max_det_diff <= 1e-11);

        const std::size_t probe{(moved + 1U) % n};
        Particles probe_particles{n};
        phy_copy_positions(particles, probe_particles);

        const double probe_dx{0.0017 * static_cast<double>((step % 4U) + 1U)};
        const double probe_dy{-0.0011 * static_cast<double>((step % 6U) + 1U)};
        const double probe_dz{0.0013 * static_cast<double>((step % 5U) + 1U)};

        probe_particles.pos_x_get()[probe] =
            wrap_coordinate(probe_particles.pos_x_get()[probe] + probe_dx, box_length);
        probe_particles.pos_y_get()[probe] =
            wrap_coordinate(probe_particles.pos_y_get()[probe] + probe_dy, box_length);
        probe_particles.pos_z_get()[probe] =
            wrap_coordinate(probe_particles.pos_z_get()[probe] + probe_dz, box_length);

        maintained.save_trig_row(probe);
        rebuilt.save_trig_row(probe);

        maintained.update_trig_cache(probe, probe_particles);
        rebuilt.update_trig_cache(probe, probe_particles);

        const double* const maintained_probe_row{maintained.build_row(probe)};
        const double* const rebuilt_probe_row{rebuilt.build_row(probe)};

        const double maintained_probe_ratio{
            maintained.determinant_ratio(probe, maintained_probe_row)};
        const double rebuilt_probe_ratio{
            rebuilt.determinant_ratio(probe, rebuilt_probe_row)};

        maintained.restore_trig_row(probe);
        rebuilt.restore_trig_row(probe);

        INFO("Maintained and rebuilt Slater states should predict the same future determinant ratio.");
        CAPTURE(step, moved, probe, maintained_probe_ratio, rebuilt_probe_ratio);
        REQUIRE(std::isfinite(maintained_probe_ratio));
        REQUIRE(std::isfinite(rebuilt_probe_ratio));
        REQUIRE(std::abs(maintained_probe_ratio - rebuilt_probe_ratio) <= 1e-8);
    }
}

TEST_CASE("WaveFunction log_psi and derivatives are invariant under integer box translations",
          "[rigorous][wavefunction]") {
    constexpr std::size_t n{7U};
    constexpr double box_length{9.5};

    Particles particles{n};
    set_stable_closed_shell_positions(particles);

    WaveFunction baseline{particles, box_length};
    const double baseline_log_psi{baseline.evaluate_log_psi(particles)};
    baseline.evaluate_derivatives(particles);

    std::vector<double> baseline_grad_x(n);
    std::vector<double> baseline_grad_y(n);
    std::vector<double> baseline_grad_z(n);
    std::vector<double> baseline_lap(n);

    for (std::size_t i = 0; i < n; ++i) {
        baseline_grad_x[i] = particles.grad_log_psi_x_get()[i];
        baseline_grad_y[i] = particles.grad_log_psi_y_get()[i];
        baseline_grad_z[i] = particles.grad_log_psi_z_get()[i];
        baseline_lap[i] = particles.lap_log_psi_get()[i];
    }

    Particles shifted{n};
    phy_copy_positions(particles, shifted);

    shifted.pos_x_get()[3] += box_length;
    shifted.pos_y_get()[3] -= 2.0 * box_length;
    shifted.pos_z_get()[3] += 3.0 * box_length;

    // After shifting:
    for (std::size_t i = 0; i < n; ++i) {
        shifted.pos_x_get()[i] -= box_length * std::floor(shifted.pos_x_get()[i] / box_length);
        shifted.pos_y_get()[i] -= box_length * std::floor(shifted.pos_y_get()[i] / box_length);
        shifted.pos_z_get()[i] -= box_length * std::floor(shifted.pos_z_get()[i] / box_length);
    }

    WaveFunction translated{shifted, box_length};
    const double translated_log_psi{translated.evaluate_log_psi(shifted)};
    translated.evaluate_derivatives(shifted);

    INFO("WaveFunction must be periodic under integer box translations.");
    CAPTURE(baseline_log_psi, translated_log_psi);
    REQUIRE(std::abs(baseline_log_psi - translated_log_psi) <= 1e-12);

    for (std::size_t i = 0; i < n; ++i) {
        CAPTURE(i);
        REQUIRE(std::abs(baseline_grad_x[i] - shifted.grad_log_psi_x_get()[i]) <= 1e-10);
        REQUIRE(std::abs(baseline_grad_y[i] - shifted.grad_log_psi_y_get()[i]) <= 1e-10);
        REQUIRE(std::abs(baseline_grad_z[i] - shifted.grad_log_psi_z_get()[i]) <= 1e-10);
        REQUIRE(std::abs(baseline_lap[i] - shifted.lap_log_psi_get()[i]) <= 1e-7);
    }
}

TEST_CASE("EnergyTracker matches an exact Ewald reference on many random small configurations",
          "[rigorous][energy]") {
    constexpr std::size_t n{3U};
    constexpr double box_length{8.25};
    constexpr std::size_t samples{25U};

    std::mt19937_64 rng{987654321ULL};
    std::uniform_real_distribution<double> dist{0.0, box_length};

    for (std::size_t sample = 0; sample < samples; ++sample) {
        Particles particles{n};

        for (std::size_t i = 0; i < n; ++i) {
            particles.pos_x_get()[i] = dist(rng);
            particles.pos_y_get()[i] = dist(rng);
            particles.pos_z_get()[i] = dist(rng);
        }

        EnergyTracker tracker{box_length, static_cast<double>(n)};
        tracker.initialize_structure_factors(particles);
        tracker.initialize_reciprocal_energy();
        tracker.initialize_real_energy(particles);

        const double tracker_total{tracker.eval_total_energy(particles)};
        const double exact_total{exact_total_potential(particles, box_length)};

        INFO("EnergyTracker total potential should match the exact Ewald reference on many random small configurations.");
        CAPTURE(sample, tracker_total, exact_total);
        REQUIRE(std::abs(tracker_total - exact_total) <= 8e-6);
    }
}

TEST_CASE("EnergyTracker remains close to the exact Ewald reference across many cached updates",
          "[rigorous][energy]") {
    constexpr std::size_t n{3U};
    constexpr double box_length{8.75};
    constexpr std::size_t steps{48U};

    Particles particles{n};
    particles.pos_x_get()[0] = 0.55;
    particles.pos_y_get()[0] = 1.20;
    particles.pos_z_get()[0] = 2.10;

    particles.pos_x_get()[1] = 3.15;
    particles.pos_y_get()[1] = 2.45;
    particles.pos_z_get()[1] = 1.05;

    particles.pos_x_get()[2] = 5.60;
    particles.pos_y_get()[2] = 0.85;
    particles.pos_z_get()[2] = 4.45;

    EnergyTracker tracker{box_length, static_cast<double>(n)};
    tracker.initialize_structure_factors(particles);
    tracker.initialize_reciprocal_energy();
    tracker.initialize_real_energy(particles);

    for (std::size_t step = 0; step < steps; ++step) {
        const std::size_t moved{step % n};

        const double old_x{particles.pos_x_get()[moved]};
        const double old_y{particles.pos_y_get()[moved]};
        const double old_z{particles.pos_z_get()[moved]};

        const double dx{0.017 * static_cast<double>((step % 4U) + 1U)};
        const double dy{-0.012 * static_cast<double>((step % 5U) + 1U)};
        const double dz{0.010 * static_cast<double>((step % 3U) + 1U)};

        particles.pos_x_get()[moved] = wrap_coordinate(old_x + dx, box_length);
        particles.pos_y_get()[moved] = wrap_coordinate(old_y + dy, box_length);
        particles.pos_z_get()[moved] = wrap_coordinate(old_z + dz, box_length);

        tracker.update_structure_factors(old_x, old_y, old_z, particles.pos_x_get()[moved],
                                         particles.pos_y_get()[moved],
                                         particles.pos_z_get()[moved]);
        tracker.update_real_energy(moved, old_x, old_y, old_z, particles);

        const double tracker_total{tracker.eval_total_energy(particles)};
        const double exact_total{exact_total_potential(particles, box_length)};

        INFO("Cached Ewald updates drifted too far from an exact recomputation.");
        CAPTURE(step, moved, tracker_total, exact_total);
        REQUIRE(std::abs(tracker_total - exact_total) <= 8e-6);
    }
}