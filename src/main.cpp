#include <exception>
#include <iostream>
#include <stdexcept>

#include "config/config.hpp"
#include "particles/particles.hpp"

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif

int main() {
    // Example:
    Config cfg{};
    Particles electrons{cfg.numParticles};

    // Never use the particles.ptr() in loops
    double* RESTRICT px{electrons.posX()};
    double const* RESTRICT px{electrons.posX()};

    return 0;
}
