#pragma once

#include "../particles/particles.hpp"
#include "../pbc/pbc.hpp"

#if defined(__GNUC__) || defined(__clang__)
#define RESTRICT __restrict__
#elif defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

class JastrowPade {
private:
    double a_{};
    double b_{};
public:
    // eqn (28) on paper
    // u(r) = (a*r) / (1 + b*r)
    explicit JastrowPade(double a = 0.5, double b = 1) noexcept
    : a_{a}
    , b_{b}
    { }

    // Getters:
    [[nodiscard]] double a() const { return a_; }
    [[nodiscard]] double b() const { return b_; }

    // eqn (27) on paper
    // J(R) = sum_{i<j} u(r_ij)
    [[nodiscard]] double value(
        const Particles& particles,
        const PeriodicBoundaryCondition& pbc
    ) const noexcept;

    /*
     eqn (29 & 30) on paper
     add jastrow contributions into provided derivative buffers:
     grad X/Y/Z is ∇_i J, lap is ∇_i^2 J (per particle i)
    */
    void addDerivatives(
        const Particles& particles,
        const PeriodicBoundaryCondition& pbc,
        double* RESTRICT gradX,
        double* RESTRICT gradY,
        double* RESTRICT gradZ,
        double* RESTRICT lap
    ) const noexcept;
};
