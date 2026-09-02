#include <catch2/catch_test_macros.hpp>

#include "particles/particles.hpp"
#include "slater_plane_wave/slater_plane_wave.hpp"

#include <array>
#include <cmath>

TEST_CASE("Batched Slater walkers match independent calculations", "[slater][walkers]") {
  constexpr auto num_particles{3uz};
  constexpr auto num_walkers{2uz};
  constexpr auto box_length{11.0_fp};
  constexpr std::array walker_positions{
    std::array{
      std::array{0.3_fp, 0.4_fp, 0.5_fp},
      std::array{1.7_fp, 0.2_fp, 2.1_fp},
      std::array{2.2_fp, 1.8_fp, 0.9_fp}
    },
    std::array{
      std::array{0.8_fp, 1.1_fp, 0.3_fp},
      std::array{2.7_fp, 0.4_fp, 1.9_fp},
      std::array{4.2_fp, 3.5_fp, 2.6_fp}
    }
  };

  Particles particles{num_particles, num_walkers};
  const auto set_positions{
    [&](Particles& target, std::size_t target_walker, std::size_t source_walker) {
      auto pos{target.pos(target_walker)};
      for (auto axis{idx(Axis::X)}; axis < idx(Axis::NUM); ++axis) {
        for (auto particle{0uz}; particle < num_particles; ++particle) {
          pos[axis][particle] = walker_positions[source_walker][particle][axis];
        }
      }
    }
  };

  for (auto walker{0uz}; walker < num_walkers; ++walker) {
    set_positions(particles, walker, walker);
  }

  SlaterPlaneWave batched_slater{particles, box_length};
  std::array<fp_t, num_walkers> batched_log_dets{};

  for (auto walker{0uz}; walker < num_walkers; ++walker) {
    Particles independent_particles{num_particles};
    set_positions(independent_particles, 0uz, walker);
    SlaterPlaneWave independent_slater{independent_particles, box_length};

    batched_log_dets[walker] = batched_slater.log_abs_det(particles.view(walker), walker);
    const auto independent_log_det{
      independent_slater.log_abs_det(independent_particles.view())
    };

    CAPTURE(walker, batched_log_dets[walker], independent_log_det);
    REQUIRE(std::isfinite(batched_log_dets[walker]));
    REQUIRE(std::abs(batched_log_dets[walker] - independent_log_det) <= 1e-12_fp);
  }

  REQUIRE(std::abs(batched_log_dets[0uz] - batched_log_dets[1uz]) > 1e-12_fp);
}
