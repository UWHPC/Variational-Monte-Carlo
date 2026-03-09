#pragma once

#include "../particles/particles.hpp"

class JastrowPade {
private:
    // Box length:
    double box_length_;

    // Jastrow-Pade parameters:
    double a_;
    double b_;

public:
    // eqn (28) on paper
    // u(r) = (a*r) / (1 + b*r)
    explicit JastrowPade(double box_length, double a = 0.25, double b = 1) noexcept
        : box_length_{box_length}, a_{a}, b_{b} {}

    // Getters:
    [[nodiscard]] double a_get() const { return a_; }
    [[nodiscard]] double b_get() const { return b_; }

    // eqn (27) on paper
    // J(R) = sum_{i<j} u(r_ij)
    [[nodiscard]] double value(const Particles& particles) const noexcept;

    /*
     eqn (29 & 30) on paper
     add jastrow contributions into provided derivative buffers:
     grad X/Y/Z is ∇_i J, lap is ∇_i^2 J (per particle i)
    */
    void add_derivatives(const Particles& particles, double* RESTRICT grad_x, double* RESTRICT grad_y,
                         double* RESTRICT grad_z, double* RESTRICT laplacian) const noexcept;
};
