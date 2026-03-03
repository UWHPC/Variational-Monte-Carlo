    #include "jastrow_pade.hpp"
    #include <cmath>
    #include <cstddef>

    double JastrowPade::value(const Particles<>& particles, const periodicBoundaryCondition& pbc) const noexcept
    {
        // NOTE: assumes gradX/gradY/gradZ/lap are zero-initialized by caller
        const std::size_t N{ particles.numParticles() };

        const double* x{ particles.posX() };
        const double* y{ particles.posY() };
        const double* z{ particles.posZ() };

        double J{ 0.0 };

        for (std::size_t i{ 0 }; i < N; ++i) {
            for (std::size_t j{ i + 1 }; j < N; ++j) {
                double dx{ 0.0 }, dy{ 0.0 }, dz{ 0.0 };
                pbc.displacement(x[i], y[i], z[i], 
                                x[j], y[j], z[j], 
                                dx, dy, dz);
                const double r2{ dx*dx + dy*dy + dz*dz };
                if ( r2 < 1e-24) continue; // according to chat this is defensive: avoid pathological coincident particles
                const double r{ std::sqrt(r2) };

                // u(r) = a*r / (1 + b*r)
                const double denom{ 1.0 + b_ * r };
                J += (a_ * r) / denom;
            }
        }

        return J;
    }

    void JastrowPade::addDerivatives(const Particles<>& particles, const periodicBoundaryCondition& pbc, double* gradX, double* gradY, double* gradZ, double* lap) const noexcept
    {
        const std::size_t N{ particles.numParticles() };

        const double* x{ particles.posX() };
        const double* y{ particles.posY() };
        const double* z{ particles.posZ() };

        for (std::size_t i{ 0 }; i < N; ++i) {
            for (std::size_t j{ i + 1 }; j < N; ++j ) {
                double dx{ 0.0 }, dy{ 0.0 }, dz{ 0.0 };
                pbc.displacement(x[i], y[i], z[i], 
                                x[j], y[j], z[j], 
                                dx, dy, dz);
                
                const double r2{ dx*dx + dy*dy + dz*dz };
                if ( r2 < 1e-24) continue; // 1/r will blow up if r2 < 1e-24.
                const double r{ std::sqrt(r2) };
                const double inv_r{ 1.0 / r};

                // u(r) = a*r / (1 + b*r)
                // u'(r) = a / (1 + b*r)^2
                // u''(r) = -2ab / (1 + b*r)^3
                const double denom{ 1.0 + b_ * r};
                const double denom2{ denom * denom };
                const double denom3{ denom2 * denom };

                const double uprime{ a_ / denom2 };
                const double usecond{ -2.0 * a_ * b_ / denom3 };

                // ∇_i u(r_ij) = u'(r) * (r_vec / r)
                const double fx{ uprime * dx * inv_r };
                const double fy{ uprime * dy * inv_r };
                const double fz{ uprime * dz * inv_r };

                gradX[i] += fx; gradY[i] += fy; gradZ[i] += fz;
                gradX[j] -= fx; gradY[j] -= fy; gradZ[j] -= fz;

                // ∇^2 u(r) = u''(r) + (2/r) u'(r)
                const double lap_pair{ usecond + 2.0 * uprime * inv_r };
                
                lap[i] += lap_pair;
                lap[j] += lap_pair;
            }
        }
        
    }