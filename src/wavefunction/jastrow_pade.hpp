#pragma once
#include "particles/particles.hpp"
#include "pbc/pbc.hpp"

class JastrowPade {
private:
    double a_{};
    double b_{};
public:
    // eqn (28) on paper 
    // u(r) = (a*r) / (1 + b*r) 
    explicit JastrowPade(double a, double b) noexcept : a_(a), b_(b) {}

    // eqn (27) on paper
    // J(R) = sum_{i<j} u(r_ij)
    [[nodiscard]]double value(const Particles<>& particles, const periodicBoundaryCondition& pbc) const noexcept;

    /*
     eqn (29 & 30) on paper
     add jastrow contributions into provided derivative buffers:
     grad X/Y/Z is ∇_i J, lap is ∇_i^2 J (per particle i)
    */
    void addDerivatives(const Particles<>& particles, const periodicBoundaryCondition& pbc, double* gradX, double* gradY, double* gradZ, double* lap) const noexcept;
};