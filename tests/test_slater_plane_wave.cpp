#include "test_utilities.hpp"

#include <catch2/catch_message.hpp>

#include "slater_plane_wave/slater_plane_wave.hpp"

#include <cmath>
#include <numbers>
#include <random>
#include <tuple>
#include <vector>

namespace {
#ifdef FP_64
constexpr real_t SLATER_PRECISION_SCALE{1.0_r};
#else
constexpr real_t SLATER_PRECISION_SCALE{1e6_r};
#endif
} // namespace

TEST_CASE("SlaterPlaneWave constructor initializes correctly", "[slater]") {
  constexpr std::size_t N{3U};
  Particles particles{N};
  SlaterPlaneWave slater{particles, 5.0_r};

  REQUIRE(slater.num_orbitals() == N);
  REQUIRE(slater.box_length() == 5.0_r);
  REQUIRE(slater.matrix_row_stride() >= N);
  REQUIRE(slater.matrix_size() == slater.matrix_row_stride() * N);
}

TEST_CASE("log_abs_det handles the N=1 constant orbital case", "[slater]") {
  // N=1: orbital 0 is k=0, cos → D = cos(0) = 1
  Particles particles{1U};
  SlaterPlaneWave slater{particles, 10.0_r};

  particles.pos().x_[0] = 3.25_r;
  particles.pos().y_[0] = 1.50_r;
  particles.pos().z_[0] = 7.75_r;

  const real_t logDet{slater.log_abs_det(particles)};

  require_near(logDet, 0.0_r);
  require_near(slater.determinant()[0], 1.0_r);
  require_near(slater.inv_determinant()[0], 1.0_r);
}

TEST_CASE("log_abs_det computes an inverse satisfying D*invD = I", "[slater]") {
  // N=3: orbital 0 = cos(0)=1, orbital 1 = cos(k1·r), orbital 2 = sin(k1·r)
  constexpr std::size_t N{3U};
  Particles particles{N};
  SlaterPlaneWave slater{particles, 11.0_r};

  particles.pos().x_[0] = 0.3_r;
  particles.pos().y_[0] = 0.4_r;
  particles.pos().z_[0] = 0.5_r;

  particles.pos().x_[1] = 1.7_r;
  particles.pos().y_[1] = 0.2_r;
  particles.pos().z_[1] = 2.1_r;

  particles.pos().x_[2] = 2.2_r;
  particles.pos().y_[2] = 1.8_r;
  particles.pos().z_[2] = 0.9_r;

  const real_t logDet{slater.log_abs_det(particles)};
  REQUIRE(std::isfinite(logDet));

  for (std::size_t row = 0; row < N; ++row) {
    for (std::size_t col = 0; col < N; ++col) {
      real_t value{};
      for (std::size_t k = 0; k < N; ++k) {
        value +=
          slater.determinant()[matrix_index(row, k, slater.matrix_row_stride())] *
          slater.inv_determinant()[matrix_index(col, k, slater.matrix_row_stride())];
      }
      const real_t expected{row == col ? 1.0_r : 0.0_r};
      require_near(value, expected, 1e-9_r * SLATER_PRECISION_SCALE);
    }
  }
}

#ifdef XPU_CUDA
TEST_CASE("log_abs_det reuses CUDA scratch across repeated N=512 calls", "[cuda-scratch]") {
  constexpr std::size_t N{512U};
  constexpr real_t L{20.0_r};
  Particles particles{N};
  SlaterPlaneWave slater{particles, L};

  std::mt19937_64 rng{0x66CADAULL};
  std::uniform_real_distribution<real_t> coordinate{0.0_r, L};
  for (std::size_t particle = 0; particle < N; ++particle) {
    particles.pos().x_[particle] = coordinate(rng);
    particles.pos().y_[particle] = coordinate(rng);
    particles.pos().z_[particle] = coordinate(rng);
  }

  const real_t first{slater.log_abs_det(particles)};
  REQUIRE(std::isfinite(first));
  REQUIRE(std::abs(first) > 1e-3_r);

  const real_t tolerance{
    10.0_r * DEFAULT_TOLERANCE * std::max(1.0_r, std::abs(first))
  };

  for (std::size_t call = 0; call < 4U; ++call) {
    CAPTURE(call, tolerance);
    const real_t repeated{slater.log_abs_det(particles)};
    REQUIRE(std::isfinite(repeated));
    require_near(repeated, first, tolerance);
  }
}
#endif

TEST_CASE("SlaterPlaneWave zero-initializes the full trig cache span", "[slater]") {
  constexpr std::size_t N{7U};
  Particles particles{N};
  SlaterPlaneWave slater{particles, 9.0_r};

  const std::size_t numK{slater.num_unique_k()};
  const std::size_t ROW_STRIDE{slater.trig_row_stride()};

  INFO("Checking full trig-cache initialization across all particle/k entries.");
  CAPTURE(N, numK, ROW_STRIDE);

  for (std::size_t p = 0; p < N; ++p) {
    for (std::size_t k = 0; k < numK; ++k) {
      const std::size_t idx{p * ROW_STRIDE + k};
      CAPTURE(p, k, idx);
      REQUIRE(slater.sin_cache()[idx] == 0.0_r);
      REQUIRE(slater.cos_cache()[idx] == 0.0_r);
    }
  }
}

TEST_CASE("determinant_ratio matches exact determinant ratio for a moved row", "[slater]") {
  constexpr std::size_t N{3U};
  constexpr real_t L{12.0_r};
  Particles particles{N};

  particles.pos().x_[0] = 0.8_r;
  particles.pos().y_[0] = 1.1_r;
  particles.pos().z_[0] = 0.3_r;

  particles.pos().x_[1] = 2.7_r;
  particles.pos().y_[1] = 0.4_r;
  particles.pos().z_[1] = 1.9_r;

  particles.pos().x_[2] = 4.2_r;
  particles.pos().y_[2] = 3.5_r;
  particles.pos().z_[2] = 2.6_r;

  SlaterPlaneWave slater{particles, L};
  const real_t logDetOld{slater.log_abs_det(particles)};
  REQUIRE(std::isfinite(logDetOld));

  const real_t detOld{determinant_3x3(slater.determinant(), slater.matrix_row_stride())};
  REQUIRE(std::abs(detOld) > (1e-12_r * SLATER_PRECISION_SCALE));

  constexpr std::size_t moved{1U};
  particles.pos().x_[moved] += 0.35_r;
  particles.pos().y_[moved] -= 0.20_r;
  particles.pos().z_[moved] += 0.15_r;

  slater.update_trig_cache(moved, particles);
  const real_t* const newRow{slater.build_row(moved)};
  const real_t ratio{slater.determinant_ratio(moved, newRow)};

  SlaterPlaneWave exactSlater{particles, L};
  const real_t logDetNew{exactSlater.log_abs_det(particles)};
  REQUIRE(std::isfinite(logDetNew));

  const real_t detNew{
    determinant_3x3(
      exactSlater.determinant(),
      exactSlater.matrix_row_stride()
    )
  };
  const real_t exactRatio{detNew / detOld};

  INFO("Trial-move determinant ratio should equal det(D_new) / det(D_old).");
  CAPTURE(detOld, detNew, ratio, exactRatio, moved);
  require_near(ratio, exactRatio, 1e-10_r * SLATER_PRECISION_SCALE);
}

TEST_CASE("accept_move matches a fresh full rebuild after an accepted row update", "[slater]") {
  constexpr std::size_t N{3U};
  constexpr real_t L{10.5_r};
  Particles particles{N};

  particles.pos().x_[0] = 0.4_r;
  particles.pos().y_[0] = 1.5_r;
  particles.pos().z_[0] = 2.7_r;

  particles.pos().x_[1] = 3.1_r;
  particles.pos().y_[1] = 2.2_r;
  particles.pos().z_[1] = 0.9_r;

  particles.pos().x_[2] = 4.8_r;
  particles.pos().y_[2] = 0.7_r;
  particles.pos().z_[2] = 3.3_r;

  SlaterPlaneWave updated{particles, L};
  const real_t logDetInitial{updated.log_abs_det(particles)};
  REQUIRE(std::isfinite(logDetInitial));

  constexpr std::size_t moved{2U};
  particles.pos().x_[moved] -= 0.28_r;
  particles.pos().y_[moved] += 0.19_r;
  particles.pos().z_[moved] += 0.41_r;

  updated.update_trig_cache(moved, particles);
  const real_t* const newRow{updated.build_row(moved)};
  const real_t ratio{updated.determinant_ratio(moved, newRow)};

  INFO("Accepted-move update should preserve determinant/inverse consistency.");
  CAPTURE(moved, ratio);
  REQUIRE(std::isfinite(ratio));
  REQUIRE(std::abs(ratio) > (1e-10_r * SLATER_PRECISION_SCALE));

  updated.accept_move(moved, newRow, ratio);

  SlaterPlaneWave rebuilt{particles, L};
  const real_t logDetRebuilt{rebuilt.log_abs_det(particles)};
  REQUIRE(std::isfinite(logDetRebuilt));

  const std::size_t S{updated.matrix_row_stride()};
  for (std::size_t row = 0; row < N; ++row) {
    for (std::size_t col = 0; col < N; ++col) {
      const std::size_t idx{row * S + col};
      CAPTURE(row, col);
      require_near(updated.determinant()[idx], rebuilt.determinant()[idx], 1e-12_r * SLATER_PRECISION_SCALE);
      require_near(updated.inv_determinant()[idx], rebuilt.inv_determinant()[idx], 1e-9_r * SLATER_PRECISION_SCALE);
    }
  }
}

TEST_CASE("N=3 determinant matrix uses cos/sin basis correctly", "[slater]") {
  // N=3: orbital 0 = cos(0)=1, orbital 1 = cos(k1·r), orbital 2 = sin(k1·r)
  // k1 is the first nonzero canonical k-vector
  constexpr std::size_t N{3U};
  constexpr real_t L{10.0_r};
  Particles particles{N};
  SlaterPlaneWave slater{particles, L};

  particles.pos().x_[0] = 1.0_r;
  particles.pos().y_[0] = 2.0_r;
  particles.pos().z_[0] = 3.0_r;

  particles.pos().x_[1] = 4.0_r;
  particles.pos().y_[1] = 5.0_r;
  particles.pos().z_[1] = 6.0_r;

  particles.pos().x_[2] = 7.0_r;
  particles.pos().y_[2] = 8.0_r;
  particles.pos().z_[2] = 9.0_r;

  slater.log_abs_det(particles);

  const auto& k_index{slater.orbital_k_index()};
  const auto& o_type{slater.orbital_type()};

  // Orbital 0 should be cos type with k=0
  REQUIRE(o_type[0] == 0);
  require_near(slater.k_vector().x_[k_index[0]], 0.0_r);
  require_near(slater.k_vector().y_[k_index[0]], 0.0_r);
  require_near(slater.k_vector().z_[k_index[0]], 0.0_r);

  // First column should all be 1.0 (cos(0))
  for (std::size_t i = 0; i < N; ++i) {
    require_near(slater.determinant()[matrix_index(i, 0, slater.matrix_row_stride())], 1.0_r);
  }

  // Orbital 1 should be cos, orbital 2 should be sin, same k-vector
  REQUIRE(o_type[1] == 0);
  REQUIRE(o_type[2] == 1);
  REQUIRE(k_index[1] == k_index[2]);

  // Verify D entries for orbitals 1 and 2
  const std::size_t ki{k_index[1]};
  const real_t kx{slater.k_vector().x_[ki]};
  const real_t ky{slater.k_vector().y_[ki]};
  const real_t kz{slater.k_vector().z_[ki]};

  for (std::size_t i = 0; i < N; ++i) {
    const real_t k_dot_r{
      kx * particles.pos().x_[i] +
      ky * particles.pos().y_[i] +
      kz * particles.pos().z_[i]
    };
    require_near(
      slater.determinant()[matrix_index(i, 1, slater.matrix_row_stride())],
      std::cos(k_dot_r)
    );
    require_near(
      slater.determinant()[matrix_index(i, 2, slater.matrix_row_stride())],
      std::sin(k_dot_r)
    );
  }
}

TEST_CASE("N=7 determinant is nonzero with cos/sin basis", "[slater]") {
  // N=7 is a closed shell: 1 (k=0) + 3 pairs × 2 = 7
  constexpr std::size_t N{7U};
  Particles particles{N};
  SlaterPlaneWave slater{particles, 10.0_r};

  // Spread particles around the box
  for (std::size_t i = 0; i < N; ++i) {
    particles.pos().x_[i] = 1.0_r + static_cast<real_t>(i) * 1.1_r;
    particles.pos().y_[i] = 0.5_r + static_cast<real_t>(i) * 0.7_r;
    particles.pos().z_[i] = 0.3_r + static_cast<real_t>(i) * 1.3_r;
  }

  const real_t logDet{slater.log_abs_det(particles)};
  REQUIRE(std::isfinite(logDet));

  // Verify D * D^{-1} = I
  for (std::size_t row = 0; row < N; ++row) {
    for (std::size_t col = 0; col < N; ++col) {
      real_t value{};
      for (std::size_t k = 0; k < N; ++k) {
        value +=
          slater.determinant()[matrix_index(row, k, slater.matrix_row_stride())] *
          slater.inv_determinant()[matrix_index(col, k, slater.matrix_row_stride())];
      }
      const real_t expected{row == col ? 1.0_r : 0.0_r};
      require_near(value, expected, 1e-9_r * SLATER_PRECISION_SCALE);
    }
  }
}

TEST_CASE("Slater derivatives match finite-difference for N=3 cos/sin basis", "[slater]") {
  constexpr std::size_t N{3U};
  constexpr real_t L{10.0_r};
  Particles particles{N};
  SlaterPlaneWave slater{particles, L};

  particles.pos().x_[0] = 1.1_r;
  particles.pos().y_[0] = 2.3_r;
  particles.pos().z_[0] = 0.7_r;

  particles.pos().x_[1] = 4.2_r;
  particles.pos().y_[1] = 1.8_r;
  particles.pos().z_[1] = 3.5_r;

  particles.pos().x_[2] = 7.6_r;
  particles.pos().y_[2] = 5.1_r;
  particles.pos().z_[2] = 8.9_r;

  // Compute analytic derivatives
  slater.log_abs_det(particles);
  const std::size_t stride{particles.p_stride()};
  std::vector<real_t> gradX(stride, 0.0_r);
  std::vector<real_t> gradY(stride, 0.0_r);
  std::vector<real_t> gradZ(stride, 0.0_r);
  std::vector<real_t> lap(stride, 0.0_r);
  slater.add_derivatives(gradX.data(), gradY.data(), gradZ.data(), lap.data());

  // Finite-difference check for particle 0
#ifdef FP_64
  const real_t h{1e-5_r};
#else
  const real_t h{1e-3_r};
#endif
  const real_t center{slater.log_abs_det(particles)};

  auto shift_and_eval = [&](std::size_t p, real_t dx, real_t dy, real_t dz) {
    Particles shifted{N};
    for (std::size_t i = 0; i < N; ++i) {
      shifted.pos().x_[i] = particles.pos().x_[i];
      shifted.pos().y_[i] = particles.pos().y_[i];
      shifted.pos().z_[i] = particles.pos().z_[i];
    }
    shifted.pos().x_[p] += dx;
    shifted.pos().y_[p] += dy;
    shifted.pos().z_[p] += dz;
    return slater.log_abs_det(shifted);
  };

  for (std::size_t p = 0; p < N; ++p) {
    const real_t fd_gx{
      (
        shift_and_eval(p, h, 0.0_r, 0.0_r) -
        shift_and_eval(p, -h, 0.0_r, 0.0_r)
      ) / (2.0_r * h)
    };
    const real_t fd_gy{
      (
        shift_and_eval(p, 0.0_r, h, 0.0_r) -
        shift_and_eval(p, 0.0_r, -h, 0.0_r)
      ) / (2.0_r * h)
    };
    const real_t fd_gz{
      (
        shift_and_eval(p, 0.0_r, 0.0_r, h) -
        shift_and_eval(p, 0.0_r, 0.0_r, -h)
      ) / (2.0_r * h)
    };

    const real_t fd_lx{
      (
        shift_and_eval(p, h, 0.0_r, 0.0_r) - 2.0_r * center +
        shift_and_eval(p, -h, 0.0_r, 0.0_r)
      ) / (h * h)
    };
    const real_t fd_ly{
      (
        shift_and_eval(p, 0.0_r, h, 0.0_r) - 2.0_r * center +
        shift_and_eval(p, 0.0_r, -h, 0.0_r)
      ) / (h * h)
    };
    const real_t fd_lz{
      (
        shift_and_eval(p, 0.0_r, 0.0_r, h) - 2.0_r * center +
        shift_and_eval(p, 0.0_r, 0.0_r, -h)
      ) / (h * h)
    };

    require_near(gradX[p], fd_gx, 1e-6_r * SLATER_PRECISION_SCALE);
    require_near(gradY[p], fd_gy, 1e-6_r * SLATER_PRECISION_SCALE);
    require_near(gradZ[p], fd_gz, 1e-6_r * SLATER_PRECISION_SCALE);
    require_near(lap[p], fd_lx + fd_ly + fd_lz, 5e-4_r * SLATER_PRECISION_SCALE);
  }
}

TEST_CASE("Shell filling produces (0,0,0) as the first n-vector", "[slater]") {
  Particles p{1U};
  SlaterPlaneWave slater{p, 5.0_r};

  REQUIRE(slater.n_vector().x_[0] == 0);
  REQUIRE(slater.n_vector().y_[0] == 0);
  REQUIRE(slater.n_vector().z_[0] == 0);
}

TEST_CASE("Shell filling for N=7 uses canonical n-vectors with 4 unique k-vectors", "[slater]") {
  // N=7: 1 (k=0) + 3 nonzero k-vectors × 2 (cos,sin) = 7
  Particles p{7U};
  SlaterPlaneWave slater{p, 10.0_r};

  REQUIRE(slater.num_unique_k() == 4U);

  const int* n_x{slater.n_vector().x_};
  const int* n_y{slater.n_vector().y_};
  const int* n_z{slater.n_vector().z_};

  // k-index 0 should be (0,0,0)
  REQUIRE(n_x[0] == 0);
  REQUIRE(n_y[0] == 0);
  REQUIRE(n_z[0] == 0);

  // k-indices 1,2,3 should have |n|^2 = 1 and be canonical (first nonzero > 0)
  for (std::size_t i = 1; i < 4; ++i) {
    const int mag_sq{n_x[i] * n_x[i] + n_y[i] * n_y[i] + n_z[i] * n_z[i]};
    REQUIRE(mag_sq == 1);

    // Canonical: first nonzero component is positive
    if (n_x[i] != 0) {
      REQUIRE(n_x[i] > 0);
    } else if (n_y[i] != 0) {
      REQUIRE(n_y[i] > 0);
    } else {
      REQUIRE(n_z[i] > 0);
    }
  }
}

TEST_CASE("Shell filling orbital types alternate cos/sin for nonzero k", "[slater]") {
  Particles p{7U};
  SlaterPlaneWave slater{p, 10.0_r};

  const auto& o_type{slater.orbital_type()};
  const auto& k_index{slater.orbital_k_index()};

  // Orbital 0: cos (k=0)
  REQUIRE(o_type[0] == 0);

  // Orbitals 1-6: pairs of (cos, sin)
  for (std::size_t i = 1; i < 7; i += 2) {
    REQUIRE(o_type[i] == 0);               // cos
    REQUIRE(o_type[i + 1] == 1);           // sin
    REQUIRE(k_index[i] == k_index[i + 1]); // same k-vector
  }
}

TEST_CASE("Shell filling n-vectors are sorted by magnitude then lexicographically", "[slater]") {
  Particles p{7U};
  SlaterPlaneWave slater{p, 10.0_r};

  const std::size_t num_k{slater.num_unique_k()};
  const int* n_x{slater.n_vector().x_};
  const int* n_y{slater.n_vector().y_};
  const int* n_z{slater.n_vector().z_};

  for (std::size_t i = 0; i + 1 < num_k; ++i) {
    const int mag_sq_a{n_x[i] * n_x[i] + n_y[i] * n_y[i] + n_z[i] * n_z[i]};
    const int mag_sq_b{n_x[i + 1] * n_x[i + 1] + n_y[i + 1] * n_y[i + 1] + n_z[i + 1] * n_z[i + 1]};

    const bool magnitude_ok{mag_sq_a <= mag_sq_b};
    REQUIRE(magnitude_ok);

    if (mag_sq_a == mag_sq_b) {
      const bool lex_ok{
        std::tie(n_x[i], n_y[i], n_z[i]) <=
        std::tie(n_x[i + 1], n_y[i + 1], n_z[i + 1])
      };
      REQUIRE(lex_ok);
    }
  }
}

TEST_CASE("Shell filling k-vectors match 2pi/L times n-vectors", "[slater]") {
  constexpr std::size_t N{7U};
  constexpr real_t L{8.0_r};
  Particles p{N};
  SlaterPlaneWave slater{p, L};

  const real_t TWO_PI_OVER_L{2.0_r * std::numbers::pi_v<real_t> / L};
  const std::size_t num_k{slater.num_unique_k()};

  for (std::size_t i = 0; i < num_k; ++i) {
    require_near(
      slater.k_vector().x_[i],
      TWO_PI_OVER_L * static_cast<real_t>(slater.n_vector().x_[i])
    );
    require_near(
      slater.k_vector().y_[i],
      TWO_PI_OVER_L * static_cast<real_t>(slater.n_vector().y_[i])
    );
    require_near(
      slater.k_vector().z_[i],
      TWO_PI_OVER_L * static_cast<real_t>(slater.n_vector().z_[i])
    );
  }
}
