#include <exception>
#include <iostream>
#include <stdexcept>

#include "config/config.hpp"
#include "particles/particles.hpp"
#include "pbc/pbc.hpp"

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif

int main() {

    return 0;
}
