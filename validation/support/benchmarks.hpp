#pragma once

#include "utilities/macros.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <string_view>

namespace validation {

struct LiteratureBenchmark {
  fp_t r_s;
  fp_t comparison_tolerance;
  fp_t step_size;
  std::uint64_t seed;
};

inline constexpr auto literature_particle_count{19uz};
inline constexpr auto literature_walker_count{16uz};
inline constexpr auto literature_warmup_sweeps{200uz};
inline constexpr auto literature_measure_sweeps{1000uz};
inline constexpr auto literature_block_size{100uz};

inline constexpr std::string_view literature_citation{
  "J. P. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981), "
  "Appendix C and Table XII, doi:10.1103/PhysRevB.23.5048"
};

inline constexpr std::string_view literature_provenance{
  "Fully polarized three-dimensional homogeneous electron gas; "
  "thermodynamic-limit total energy reconstructed from the PZ81 correlation "
  "parametrization and analytical Hartree-Fock kinetic and exchange energies"
};

inline constexpr std::array literature_benchmarks{
  LiteratureBenchmark{2.0_fp, 0.07_fp, 0.3_fp, 2019uz},
  LiteratureBenchmark{5.0_fp, 0.01_fp, 0.6_fp, 5019uz},
  LiteratureBenchmark{10.0_fp, 0.005_fp, 1.0_fp, 10019uz}
};

[[nodiscard]] inline fp_t pz81_total_energy(fp_t r_s) {
  constexpr auto two{2.0_fp};
  constexpr auto three{3.0_fp};
  constexpr auto four{4.0_fp};
  constexpr auto nine{9.0_fp};
  const auto fermi_wave_vector{
    std::cbrt(nine * std::numbers::pi_v<fp_t> / two) / r_s
  };
  const auto kinetic{
    three * fermi_wave_vector * fermi_wave_vector / 10.0_fp
  };
  const auto exchange{
    -three * fermi_wave_vector / (four * std::numbers::pi_v<fp_t>)
  };

  fp_t correlation{};
  if (r_s >= 1.0_fp) {
    constexpr auto gamma{-0.0843_fp};
    constexpr auto beta_1{1.3981_fp};
    constexpr auto beta_2{0.2611_fp};
    correlation = gamma / (
      1.0_fp + beta_1 * std::sqrt(r_s) + beta_2 * r_s
    );
  } else {
    constexpr auto a{0.01555_fp};
    constexpr auto b{-0.0269_fp};
    constexpr auto c{0.0007_fp};
    constexpr auto d{-0.0048_fp};
    const auto logarithm{std::log(r_s)};
    correlation = a * logarithm + b +
      c * r_s * logarithm + d * r_s;
  }

  return kinetic + exchange + correlation;
}

[[nodiscard]] inline fp_t box_length_from_r_s(
  fp_t r_s,
  std::size_t particle_count
) {
  constexpr auto four{4.0_fp};
  constexpr auto three{3.0_fp};
  const auto volume{
    four * std::numbers::pi_v<fp_t> * r_s * r_s * r_s *
    scast<fp_t>(particle_count) / three
  };
  return std::cbrt(volume);
}

} // namespace validation
