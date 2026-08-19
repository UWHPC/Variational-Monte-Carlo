#pragma once

#include "../utilities/components.hpp"
#include "../utilities/macros.hpp"

#include <xpu/soa.hpp>

#include <cstddef>

namespace simulation {

struct StepResult {
  bool accepted;
  std::size_t moved_particle;
  xpu::array<real_t, idx(Axis::NUM)> old_pos;
  xpu::array<real_t, idx(Axis::NUM)> new_pos;
  real_t log_psi_delta;
  real_t real_energy_delta;
  real_t reciprocal_energy;
};

} // namespace simulation
