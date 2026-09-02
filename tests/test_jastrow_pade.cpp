#include "test_utilities.hpp"

#include "jastrow_pade/jastrow_pade.hpp"

#include <cstddef>
#include <vector>

namespace {

#ifdef FP_64
constexpr fp_t FD_STEP{1e-5_fp};
constexpr fp_t FD_GRAD_TOLERANCE{1e-7_fp};
constexpr fp_t FD_LAP_TOLERANCE{2e-4_fp};
constexpr fp_t GRAD_SUM_TOLERANCE{1e-12_fp};
#else
constexpr fp_t FD_STEP{2e-2_fp};
constexpr fp_t FD_GRAD_TOLERANCE{5e-4_fp};
constexpr fp_t FD_LAP_TOLERANCE{2e-2_fp};
constexpr fp_t GRAD_SUM_TOLERANCE{1e-6_fp};
#endif

fp_t valueAtOffset(
  const JastrowPade& jastrow,
  const Particles& reference,
  std::size_t particle,
  fp_t dx, fp_t dy, fp_t dz
) {
  Particles shifted{copy_particle_positions(reference)};
  shifted.pos().x_[particle] += dx;
  shifted.pos().y_[particle] += dy;
  shifted.pos().z_[particle] += dz;
  return jastrow.value(shifted.view());
}

} // namespace

TEST_CASE("Jastrow value uses minimum-image pair distances", "[jastrow]") {
  const JastrowPade jastrow{10.0_fp, 0.25_fp, 1.0_fp};
  Particles particles{2U};

  particles.pos().x_[0] = 0.1_fp;
  particles.pos().y_[0] = 0.0_fp;
  particles.pos().z_[0] = 0.0_fp;

  particles.pos().x_[1] = 9.9_fp;
  particles.pos().y_[1] = 0.0_fp;
  particles.pos().z_[1] = 0.0_fp;

  const fp_t r{0.2_fp};
  const fp_t expected{(0.25_fp * r) / (1.0_fp + r)};
  require_near(jastrow.value(particles.view()), expected);
}

TEST_CASE("Jastrow value skips degenerate pairs", "[jastrow]") {
  const JastrowPade jastrow{10.0_fp, 0.25_fp, 1.0_fp};
  Particles particles{2U};

  particles.pos().x_[0] = 1.0_fp;
  particles.pos().y_[0] = 2.0_fp;
  particles.pos().z_[0] = 3.0_fp;

  particles.pos().x_[1] = 1.0_fp;
  particles.pos().y_[1] = 2.0_fp;
  particles.pos().z_[1] = 3.0_fp;

  require_near(jastrow.value(particles.view()), 0.0_fp);
}

TEST_CASE("Jastrow derivatives match the analytic two-particle result", "[jastrow]") {
  const JastrowPade jastrow{100.0_fp, 0.25_fp, 1.0_fp};
  Particles particles{2U};

  particles.pos().x_[0] = 0.0_fp;
  particles.pos().y_[0] = 0.0_fp;
  particles.pos().z_[0] = 0.0_fp;

  particles.pos().x_[1] = 1.0_fp;
  particles.pos().y_[1] = 0.0_fp;
  particles.pos().z_[1] = 0.0_fp;

  const std::size_t stride{particles.p_stride()};
  std::vector<fp_t> gradX(stride, 0.0_fp);
  std::vector<fp_t> gradY(stride, 0.0_fp);
  std::vector<fp_t> gradZ(stride, 0.0_fp);
  std::vector<fp_t> lap(stride, 0.0_fp);

  jastrow.add_derivatives(particles.view(), gradX.data(), gradY.data(), gradZ.data(), lap.data());

  require_near(gradX[0], -0.0625_fp);
  require_near(gradX[1], 0.0625_fp);
  require_near(gradX[0] + gradX[1], 0.0_fp);

  require_near(gradY[0], 0.0_fp);
  require_near(gradY[1], 0.0_fp);
  require_near(gradZ[0], 0.0_fp);
  require_near(gradZ[1], 0.0_fp);

  // ∇²u = u'' + (2/r)u' = -0.0625 + 2(0.0625) = 0.0625
  require_near(lap[0], 0.0625_fp);
  require_near(lap[1], 0.0625_fp);

  for (std::size_t i = 2; i < stride; ++i) {
    require_near(gradX[i], 0.0_fp);
    require_near(gradY[i], 0.0_fp);
    require_near(gradZ[i], 0.0_fp);
    require_near(lap[i], 0.0_fp);
  }
}

TEST_CASE("Jastrow derivatives are unchanged for degenerate pairs", "[jastrow]") {
  const JastrowPade jastrow{10.0_fp, 0.25_fp, 1.0_fp};
  Particles particles{2U};

  particles.pos().x_[0] = 4.0_fp;
  particles.pos().y_[0] = 5.0_fp;
  particles.pos().z_[0] = 6.0_fp;

  particles.pos().x_[1] = 4.0_fp;
  particles.pos().y_[1] = 5.0_fp;
  particles.pos().z_[1] = 6.0_fp;

  const std::size_t stride{particles.p_stride()};
  std::vector<fp_t> gradX(stride, 3.0_fp);
  std::vector<fp_t> gradY(stride, -2.0_fp);
  std::vector<fp_t> gradZ(stride, 1.5_fp);
  std::vector<fp_t> lap(stride, 7.0_fp);

  jastrow.add_derivatives(particles.view(), gradX.data(), gradY.data(), gradZ.data(), lap.data());

  for (std::size_t i = 0; i < stride; ++i) {
    require_near(gradX[i], 3.0_fp);
    require_near(gradY[i], -2.0_fp);
    require_near(gradZ[i], 1.5_fp);
    require_near(lap[i], 7.0_fp);
  }
}

TEST_CASE("Jastrow derivatives match finite-difference gradients and Laplacians", "[jastrow]") {
  const JastrowPade jastrow{20.0_fp, 0.6_fp, 0.9_fp};
  Particles particles{3U};

  particles.pos().x_[0] = 1.1_fp;
  particles.pos().y_[0] = 2.2_fp;
  particles.pos().z_[0] = 0.7_fp;

  particles.pos().x_[1] = 3.8_fp;
  particles.pos().y_[1] = 1.4_fp;
  particles.pos().z_[1] = 2.5_fp;

  particles.pos().x_[2] = 0.2_fp;
  particles.pos().y_[2] = 4.1_fp;
  particles.pos().z_[2] = 3.3_fp;

  const std::size_t stride{particles.p_stride()};
  std::vector<fp_t> gradX(stride, 0.0_fp);
  std::vector<fp_t> gradY(stride, 0.0_fp);
  std::vector<fp_t> gradZ(stride, 0.0_fp);
  std::vector<fp_t> lap(stride, 0.0_fp);
  jastrow.add_derivatives(particles.view(), gradX.data(), gradY.data(), gradZ.data(), lap.data());

  const fp_t h{FD_STEP};
  const fp_t valueCenter{jastrow.value(particles.view())};

  const fp_t dJdx{
    (
      valueAtOffset(jastrow, particles, 0U, h, 0.0_fp, 0.0_fp) -
      valueAtOffset(jastrow, particles, 0U, -h, 0.0_fp, 0.0_fp)
    ) / (2.0_fp * h)
  };
  const fp_t dJdy{
    (
      valueAtOffset(jastrow, particles, 0U, 0.0_fp, h, 0.0_fp) -
      valueAtOffset(jastrow, particles, 0U, 0.0_fp, -h, 0.0_fp)
    ) / (2.0_fp * h)
  };
  const fp_t dJdz{
    (
      valueAtOffset(jastrow, particles, 0U, 0.0_fp, 0.0_fp, h) -
      valueAtOffset(jastrow, particles, 0U, 0.0_fp, 0.0_fp, -h)
    ) / (2.0_fp * h)
  };

  const fp_t d2Jdx2{
    (
      valueAtOffset(jastrow, particles, 0U, h, 0.0_fp, 0.0_fp) - 2.0_fp * valueCenter +
      valueAtOffset(jastrow, particles, 0U, -h, 0.0_fp, 0.0_fp)
    ) / (h * h)
  };
  const fp_t d2Jdy2{
    (
      valueAtOffset(jastrow, particles, 0U, 0.0_fp, h, 0.0_fp) - 2.0_fp * valueCenter +
      valueAtOffset(jastrow, particles, 0U, 0.0_fp, -h, 0.0_fp)
    ) / (h * h)
  };
  const fp_t d2Jdz2{
    (
      valueAtOffset(jastrow, particles, 0U, 0.0_fp, 0.0_fp, h) - 2.0_fp * valueCenter +
      valueAtOffset(jastrow, particles, 0U, 0.0_fp, 0.0_fp, -h)
    ) / (h * h)
  };

  require_near(gradX[0], dJdx, FD_GRAD_TOLERANCE);
  require_near(gradY[0], dJdy, FD_GRAD_TOLERANCE);
  require_near(gradZ[0], dJdz, FD_GRAD_TOLERANCE);
  require_near(lap[0], d2Jdx2 + d2Jdy2 + d2Jdz2, FD_LAP_TOLERANCE);

  require_near(gradX[0] + gradX[1] + gradX[2], 0.0_fp, GRAD_SUM_TOLERANCE);
  require_near(gradY[0] + gradY[1] + gradY[2], 0.0_fp, GRAD_SUM_TOLERANCE);
  require_near(gradZ[0] + gradZ[1] + gradZ[2], 0.0_fp, GRAD_SUM_TOLERANCE);
}
