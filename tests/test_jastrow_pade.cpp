#include "test_utilities.hpp"

#include "jastrow_pade/jastrow_pade.cuh"

#include <cstddef>
#include <vector>

namespace {

#ifdef FP_64
constexpr real_t FD_STEP{1e-5_r};
constexpr real_t FD_GRAD_TOLERANCE{1e-7_r};
constexpr real_t FD_LAP_TOLERANCE{2e-4_r};
constexpr real_t GRAD_SUM_TOLERANCE{1e-12_r};
#else
constexpr real_t FD_STEP{2e-2_r};
constexpr real_t FD_GRAD_TOLERANCE{5e-4_r};
constexpr real_t FD_LAP_TOLERANCE{2e-2_r};
constexpr real_t GRAD_SUM_TOLERANCE{1e-6_r};
#endif

real_t valueAtOffset(
  const JastrowPade& jastrow,
  const Particles& reference,
  std::size_t particle,
  real_t dx, real_t dy, real_t dz
) {
  Particles shifted{copy_particle_positions(reference)};
  shifted.pos().x_[particle] += dx;
  shifted.pos().y_[particle] += dy;
  shifted.pos().z_[particle] += dz;
  return jastrow.value(shifted);
}

} // namespace

TEST_CASE("Jastrow value uses minimum-image pair distances", "[jastrow]") {
  const JastrowPade jastrow{10.0_r, 0.25_r, 1.0_r};
  Particles particles{2U};

  particles.pos().x_[0] = 0.1_r;
  particles.pos().y_[0] = 0.0_r;
  particles.pos().z_[0] = 0.0_r;

  particles.pos().x_[1] = 9.9_r;
  particles.pos().y_[1] = 0.0_r;
  particles.pos().z_[1] = 0.0_r;

  const real_t r{0.2_r};
  const real_t expected{(0.25_r * r) / (1.0_r + r)};
  require_near(jastrow.value(particles), expected);
}

TEST_CASE("Jastrow value skips degenerate pairs", "[jastrow]") {
  const JastrowPade jastrow{10.0_r, 0.25_r, 1.0_r};
  Particles particles{2U};

  particles.pos().x_[0] = 1.0_r;
  particles.pos().y_[0] = 2.0_r;
  particles.pos().z_[0] = 3.0_r;

  particles.pos().x_[1] = 1.0_r;
  particles.pos().y_[1] = 2.0_r;
  particles.pos().z_[1] = 3.0_r;

  require_near(jastrow.value(particles), 0.0_r);
}

TEST_CASE("Jastrow derivatives match the analytic two-particle result", "[jastrow]") {
  const JastrowPade jastrow{100.0_r, 0.25_r, 1.0_r};
  Particles particles{2U};

  particles.pos().x_[0] = 0.0_r;
  particles.pos().y_[0] = 0.0_r;
  particles.pos().z_[0] = 0.0_r;

  particles.pos().x_[1] = 1.0_r;
  particles.pos().y_[1] = 0.0_r;
  particles.pos().z_[1] = 0.0_r;

  const std::size_t stride{particles.p_stride()};
  std::vector<real_t> gradX(stride, 0.0_r);
  std::vector<real_t> gradY(stride, 0.0_r);
  std::vector<real_t> gradZ(stride, 0.0_r);
  std::vector<real_t> lap(stride, 0.0_r);

  jastrow.add_derivatives(particles, gradX.data(), gradY.data(), gradZ.data(), lap.data());

  require_near(gradX[0], -0.0625_r);
  require_near(gradX[1], 0.0625_r);
  require_near(gradX[0] + gradX[1], 0.0_r);

  require_near(gradY[0], 0.0_r);
  require_near(gradY[1], 0.0_r);
  require_near(gradZ[0], 0.0_r);
  require_near(gradZ[1], 0.0_r);

  // ∇²u = u'' + (2/r)u' = -0.0625 + 2(0.0625) = 0.0625
  require_near(lap[0], 0.0625_r);
  require_near(lap[1], 0.0625_r);

  for (std::size_t i = 2; i < stride; ++i) {
    require_near(gradX[i], 0.0_r);
    require_near(gradY[i], 0.0_r);
    require_near(gradZ[i], 0.0_r);
    require_near(lap[i], 0.0_r);
  }
}

TEST_CASE("Jastrow derivatives are unchanged for degenerate pairs", "[jastrow]") {
  const JastrowPade jastrow{10.0_r, 0.25_r, 1.0_r};
  Particles particles{2U};

  particles.pos().x_[0] = 4.0_r;
  particles.pos().y_[0] = 5.0_r;
  particles.pos().z_[0] = 6.0_r;

  particles.pos().x_[1] = 4.0_r;
  particles.pos().y_[1] = 5.0_r;
  particles.pos().z_[1] = 6.0_r;

  const std::size_t stride{particles.p_stride()};
  std::vector<real_t> gradX(stride, 3.0_r);
  std::vector<real_t> gradY(stride, -2.0_r);
  std::vector<real_t> gradZ(stride, 1.5_r);
  std::vector<real_t> lap(stride, 7.0_r);

  jastrow.add_derivatives(particles, gradX.data(), gradY.data(), gradZ.data(), lap.data());

  for (std::size_t i = 0; i < stride; ++i) {
    require_near(gradX[i], 3.0_r);
    require_near(gradY[i], -2.0_r);
    require_near(gradZ[i], 1.5_r);
    require_near(lap[i], 7.0_r);
  }
}

TEST_CASE("Jastrow derivatives match finite-difference gradients and Laplacians", "[jastrow]") {
  const JastrowPade jastrow{20.0_r, 0.6_r, 0.9_r};
  Particles particles{3U};

  particles.pos().x_[0] = 1.1_r;
  particles.pos().y_[0] = 2.2_r;
  particles.pos().z_[0] = 0.7_r;

  particles.pos().x_[1] = 3.8_r;
  particles.pos().y_[1] = 1.4_r;
  particles.pos().z_[1] = 2.5_r;

  particles.pos().x_[2] = 0.2_r;
  particles.pos().y_[2] = 4.1_r;
  particles.pos().z_[2] = 3.3_r;

  const std::size_t stride{particles.p_stride()};
  std::vector<real_t> gradX(stride, 0.0_r);
  std::vector<real_t> gradY(stride, 0.0_r);
  std::vector<real_t> gradZ(stride, 0.0_r);
  std::vector<real_t> lap(stride, 0.0_r);
  jastrow.add_derivatives(particles, gradX.data(), gradY.data(), gradZ.data(), lap.data());

  const real_t h{FD_STEP};
  const real_t valueCenter{jastrow.value(particles)};

  const real_t dJdx{
    (
      valueAtOffset(jastrow, particles, 0U, h, 0.0_r, 0.0_r) -
      valueAtOffset(jastrow, particles, 0U, -h, 0.0_r, 0.0_r)
    ) / (2.0_r * h)
  };
  const real_t dJdy{
    (
      valueAtOffset(jastrow, particles, 0U, 0.0_r, h, 0.0_r) -
      valueAtOffset(jastrow, particles, 0U, 0.0_r, -h, 0.0_r)
    ) / (2.0_r * h)
  };
  const real_t dJdz{
    (
      valueAtOffset(jastrow, particles, 0U, 0.0_r, 0.0_r, h) -
      valueAtOffset(jastrow, particles, 0U, 0.0_r, 0.0_r, -h)
    ) / (2.0_r * h)
  };

  const real_t d2Jdx2{
    (
      valueAtOffset(jastrow, particles, 0U, h, 0.0_r, 0.0_r) - 2.0_r * valueCenter +
      valueAtOffset(jastrow, particles, 0U, -h, 0.0_r, 0.0_r)
    ) / (h * h)
  };
  const real_t d2Jdy2{
    (
      valueAtOffset(jastrow, particles, 0U, 0.0_r, h, 0.0_r) - 2.0_r * valueCenter +
      valueAtOffset(jastrow, particles, 0U, 0.0_r, -h, 0.0_r)
    ) / (h * h)
  };
  const real_t d2Jdz2{
    (
      valueAtOffset(jastrow, particles, 0U, 0.0_r, 0.0_r, h) - 2.0_r * valueCenter +
      valueAtOffset(jastrow, particles, 0U, 0.0_r, 0.0_r, -h)
    ) / (h * h)
  };

  require_near(gradX[0], dJdx, FD_GRAD_TOLERANCE);
  require_near(gradY[0], dJdy, FD_GRAD_TOLERANCE);
  require_near(gradZ[0], dJdz, FD_GRAD_TOLERANCE);
  require_near(lap[0], d2Jdx2 + d2Jdy2 + d2Jdz2, FD_LAP_TOLERANCE);

  require_near(gradX[0] + gradX[1] + gradX[2], 0.0_r, GRAD_SUM_TOLERANCE);
  require_near(gradY[0] + gradY[1] + gradY[2], 0.0_r, GRAD_SUM_TOLERANCE);
  require_near(gradZ[0] + gradZ[1] + gradZ[2], 0.0_r, GRAD_SUM_TOLERANCE);
}