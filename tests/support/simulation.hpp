#pragma once

#include "config/config.hpp"
#include "output_writer/output_writer.hpp"

#include <optional>
#include <utility>
#include <vector>

class RecordingOutputWriter final : public OutputWriter {
public:
  std::optional<InitData> init{};
  std::vector<FrameData> frames{};
  std::optional<DoneData> done{};

  void write_init(const InitData& data) override {
    init = data;
  }

  void write_frame(const FrameData& data) override {
    frames.push_back(data);
  }

  void write_done(const DoneData& data) override {
    done = data;
  }
};

inline Config make_config(
  std::size_t num_particles,
  fp_t box_length,
  std::size_t warmup_sweeps,
  std::size_t measure_sweeps,
  fp_t step_size,
  std::uint64_t master_seed,
  std::size_t block_size,
  std::size_t num_walkers = 1uz
) {
  Config config{};
  config.num_particles = num_particles;
  config.num_walkers = num_walkers;
  config.warmup_sweeps = warmup_sweeps;
  config.measure_sweeps = measure_sweeps;
  config.box_length = box_length;
  config.step_size = step_size;
  config.master_seed = master_seed;
  config.block_size = block_size;
  config.is_master_thread = false;
  config.warmup_steps = num_particles * num_walkers * warmup_sweeps;
  config.measure_steps = num_particles * num_walkers * measure_sweeps;
  return config;
}
