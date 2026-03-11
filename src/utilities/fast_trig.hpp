#pragma once

// 6th-order Taylor polynomial approximations for sin and cos.

inline double fast_cos(double x) noexcept {
    const double x2 = x * x;
    const double x4 = x2 * x2;
    const double x6 = x4 * x2;
    return 1.0 - x2 * (1.0 / 2.0) + x4 * (1.0 / 24.0) - x6 * (1.0 / 720.0);
}

inline double fast_sin(double x) noexcept {
    const double x2 = x * x;
    const double x3 = x2 * x;
    const double x5 = x3 * x2;
    return x - x3 * (1.0 / 6.0) + x5 * (1.0 / 120.0);
}
