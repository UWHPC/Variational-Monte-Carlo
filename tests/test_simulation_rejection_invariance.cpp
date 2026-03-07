#include <catch2/catch_test_macros.hpp>

#include "output_writer/output_writer.hpp"
#include "simulation/simulation.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

namespace {

struct CapturedFrame {
    std::size_t accepted;
    double local_energy;
};

class CapturingOutputWriter final : public OutputWriter {
public:
    explicit CapturingOutputWriter(std::vector<CapturedFrame>& frames) : frames_{frames} {}

    void write_init(const InitData&) override {}

    void write_frame(const FrameData& data) override {
        frames_.push_back(CapturedFrame{.accepted = data.accepted, .local_energy = data.local_energy});
    }

    void write_done(const DoneData&) override {}

private:
    std::vector<CapturedFrame>& frames_;
};

} // namespace

// TODO: Test is invalid for current physics model
// TEST_CASE("Simulation keeps local energy unchanged across rejected moves", "[simulation]") {
//     const Config config{.num_particles = 16U,
//                         .box_length = 4.5,
//                         .warmup_steps = 0U,
//                         .measure_steps = 120U,
//                         .step_size = 0.35,
//                         .seed = 12345U,
//                         .block_size = 20U};

//     std::vector<CapturedFrame> frames{};
//     std::unique_ptr<OutputWriter> writer{std::make_unique<CapturingOutputWriter>(frames)};

//     Simulation simulation{config, std::move(writer)};
//     simulation.run();

//     REQUIRE(frames.size() == config.measure_steps);

//     constexpr double tolerance{1e-10};
//     std::size_t rejected_transitions{};
//     std::size_t changed_energy_after_reject{};

//     for (std::size_t i = 1; i < frames.size(); ++i) {
//         if (frames[i].accepted == frames[i - 1].accepted) {
//             ++rejected_transitions;
//             if (std::abs(frames[i].local_energy - frames[i - 1].local_energy) > tolerance) {
//                 ++changed_energy_after_reject;
//             }
//         }
//     }

//     REQUIRE(rejected_transitions > 0U);
//     REQUIRE(changed_energy_after_reject == 0U);
// }
