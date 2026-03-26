#include "test_utilities.hpp"

#include "energy_tracking/energy_tracking.hpp"
#include "wavefunction/wavefunction.hpp"

#include <cstddef>
#include <numbers>
#include <random>

// N=1: only k=0 orbital, T_exact = 0, uniform wavefunction
TEST_CASE("Free gas N=1: kinetic energy is exactly zero", "[validation]") {
    constexpr std::size_t N{1U};
    constexpr double L{5.0};

    Particles particles{N};
    SlaterPlaneWave slater{particles, L};
    const double T_EXACT{exact_kinetic_energy(slater)};
    require_near(T_EXACT, 0.0);

    WaveFunction wf{particles, L, 0.0, 1.0}; // a = 0 -> Jastrow off

    std::mt19937_64 rng{42};
    std::uniform_real_distribution<double> uniform{0.0, L};

    for (int sample = 0; sample < 5; ++sample) {
        particles.pos_x_get()[0] = uniform(rng);
        particles.pos_y_get()[0] = uniform(rng);
        particles.pos_z_get()[0] = uniform(rng);

        wf.evaluate_log_psi(particles);
        wf.evaluate_derivatives(particles);

        const double T_local{local_kinetic_energy(particles)};

        require_near(T_local, 0.0);
    }
}

// N=7 closed shell: T_exact = 3(2pi/L)^2
TEST_CASE("Free gas N=7: local kinetic energy matches exact value at every sample", "[validation]") {
    constexpr std::size_t N{7U};
    constexpr double L{6.0};

    Particles particles{N};
    SlaterPlaneWave slater{particles, L};
    const double T_EXACT{exact_kinetic_energy(slater)};

    const double TWO_PI_OVER_L{2.0 * std::numbers::pi / L};
    const double T_ANALYTICAL{3.0 * TWO_PI_OVER_L * TWO_PI_OVER_L};
    require_near(T_EXACT, T_ANALYTICAL);

    WaveFunction wf{particles, L, 0.0, 1.0};

    std::mt19937_64 rng{314159};
    std::uniform_real_distribution<double> uniform{0.0, L};

    for (int sample = 0; sample < 10; ++sample) {
        for (std::size_t i = 0; i < N; ++i) {
            particles.pos_x_get()[i] = uniform(rng);
            particles.pos_y_get()[i] = uniform(rng);
            particles.pos_z_get()[i] = uniform(rng);
        }

        wf.evaluate_log_psi(particles);
        wf.evaluate_derivatives(particles);

        const double T_local{local_kinetic_energy(particles)};

        require_near(T_local, T_EXACT, 1e-8);
    }
}

// N=19 closed shell: T_exact = 15(2pi/L)^2
TEST_CASE("Free gas N=19: local kinetic energy matches exact value at every sample", "[validation]") {
    constexpr std::size_t N{19U};
    constexpr double L{7.0};

    Particles particles{N};
    SlaterPlaneWave slater{particles, L};
    const double T_EXACT{exact_kinetic_energy(slater)};

    const double K_SQ{(2.0 * std::numbers::pi / L) * (2.0 * std::numbers::pi / L)};
    const double T_ANALYTICAL{15.0 * K_SQ};
    require_near(T_EXACT, T_ANALYTICAL);

    WaveFunction wf{particles, L, 0.0, 1.0};

    std::mt19937_64 rng{271828};
    std::uniform_real_distribution<double> uniform{0.0, L};

    for (int sample = 0; sample < 10; ++sample) {
        for (std::size_t i = 0; i < N; ++i) {
            particles.pos_x_get()[i] = uniform(rng);
            particles.pos_y_get()[i] = uniform(rng);
            particles.pos_z_get()[i] = uniform(rng);
        }

        wf.evaluate_log_psi(particles);
        wf.evaluate_derivatives(particles);

        const double T_local{local_kinetic_energy(particles)};

        require_near(T_local, T_EXACT, 1e-7);
    }
}

// Jastrow off -> zero variance, local KE identical at every config
TEST_CASE("Free gas zero-variance: local kinetic energy is configuration-independent", "[validation]") {
    constexpr std::size_t N{7U};
    constexpr double L{5.5};

    Particles particles{N};
    WaveFunction wf{particles, L, 0.0, 1.0};

    std::mt19937_64 rng{999};
    std::uniform_real_distribution<double> uniform{0.0, L};

    constexpr int NUM_SAMPLES{50};
    double first_T{};

    for (int sample = 0; sample < NUM_SAMPLES; ++sample) {
        for (std::size_t i = 0; i < N; ++i) {
            particles.pos_x_get()[i] = uniform(rng);
            particles.pos_y_get()[i] = uniform(rng);
            particles.pos_z_get()[i] = uniform(rng);
        }

        wf.evaluate_log_psi(particles);
        wf.evaluate_derivatives(particles);

        const double T_local{local_kinetic_energy(particles)};

        if (sample == 0) {
            first_T = T_local;
        } else {
            require_near(T_local, first_T, 1e-8);
        }
    }
}

// EnergyTracker kinetic term should match local_kinetic_energy with Jastrow off
TEST_CASE("Free gas: EnergyTracker kinetic term matches manual computation", "[validation]") {
    constexpr std::size_t N{7U};
    constexpr double L{6.0};

    Particles particles{N};
    WaveFunction wf{particles, L, 0.0, 1.0};
    EnergyTracker tracker{L, static_cast<double>(N)};

    std::mt19937_64 rng{65537};
    std::uniform_real_distribution<double> uniform{0.0, L};

    for (std::size_t i = 0; i < N; ++i) {
        particles.pos_x_get()[i] = uniform(rng);
        particles.pos_y_get()[i] = uniform(rng);
        particles.pos_z_get()[i] = uniform(rng);
    }

    wf.evaluate_log_psi(particles);
    wf.evaluate_derivatives(particles);

    tracker.initialize_structure_factors(particles);

    const double T_manual{local_kinetic_energy(particles)};
    const SlaterPlaneWave& slater{wf.slater_plane_wave_get()};
    const double T_EXACT{exact_kinetic_energy(slater)};

    require_near(T_manual, T_EXACT, 1e-8);
}

// Verify closed-shell orbital counts: N = 1, 7, 19, 27
TEST_CASE("Shell filling produces correct closed-shell orbital counts", "[validation]") {
    constexpr double L{10.0};

    SECTION("N = 1") {
        Particles p{1U};
        const SlaterPlaneWave slater{p, L};
        REQUIRE(slater.num_unique_k_get() == 1U);
        REQUIRE(slater.num_orbitals_get() == 1U);
    }

    SECTION("N = 7") {
        Particles p{7U};
        const SlaterPlaneWave slater{p, L};
        REQUIRE(slater.num_unique_k_get() == 4U);
        REQUIRE(slater.num_orbitals_get() == 7U);
    }

    SECTION("N = 19") {
        Particles p{19U};
        const SlaterPlaneWave slater{p, L};
        REQUIRE(slater.num_unique_k_get() == 10U);
        REQUIRE(slater.num_orbitals_get() == 19U);
    }

    SECTION("N = 27") {
        Particles p{27U};
        const SlaterPlaneWave slater{p, L};
        REQUIRE(slater.num_unique_k_get() == 14U);
        REQUIRE(slater.num_orbitals_get() == 27U);
    }
}

// N=16 partial shell: zero-variance property still holds
TEST_CASE("Free gas partial shell N=16: zero-variance property still holds", "[validation]") {
    constexpr std::size_t N{16U};
    constexpr double L{6.5};

    Particles particles{N};
    WaveFunction wf{particles, L, 0.0, 1.0};
    const SlaterPlaneWave& slater{wf.slater_plane_wave_get()};
    const double T_EXACT{exact_kinetic_energy(slater)};

    std::mt19937_64 rng{1337};
    std::uniform_real_distribution<double> uniform{0.0, L};

    constexpr int NUM_SAMPLES{20};

    for (int sample = 0; sample < NUM_SAMPLES; ++sample) {
        for (std::size_t i = 0; i < N; ++i) {
            particles.pos_x_get()[i] = uniform(rng);
            particles.pos_y_get()[i] = uniform(rng);
            particles.pos_z_get()[i] = uniform(rng);
        }

        wf.evaluate_log_psi(particles);
        wf.evaluate_derivatives(particles);

        const double T_local{local_kinetic_energy(particles)};

        require_near(T_local, T_EXACT, 1e-7);
    }
}