#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

struct InitData {
    std::string run_id;

    std::size_t num_particles;
    double box_length;
    std::size_t warmup_steps;
    std::size_t measure_steps;
    double step_size;
    uint64_t seed;
    std::size_t block_size;
};

struct FrameData {
    std::size_t step;
    std::size_t accepted;
    std::size_t proposed;
    double acceptance_rate;
    double local_energy;
    double mean_energy;
    std::optional<double> standard_error;
    std::vector<double> positions; // size = 3 * num_particles
};

struct DoneData {
    std::size_t total_accepted;
    std::size_t total_proposed;
    double final_acceptance_rate;
    double final_mean_energy;
    std::optional<double> final_standard_error;
};

enum class OutputFormat { JSON, CSV, BIN };

class OutputWriter {
public:
    virtual ~OutputWriter() = default;

    virtual void write_init(const InitData& data) = 0;
    virtual void write_frame(const FrameData& data) = 0;
    virtual void write_done(const DoneData& data) = 0;
};

class JsonOutputWriter final : public OutputWriter {
public:
    explicit JsonOutputWriter(std::ostream& out);

    void write_init(const InitData& data) override;
    void write_frame(const FrameData& data) override;
    void write_done(const DoneData& data) override;

private:
    std::ostream& out_;
};

class CsvOutputWriter final : public OutputWriter {
public:
    explicit CsvOutputWriter(std::ostream& out);

    void write_init(const InitData& data) override;
    void write_frame(const FrameData& data) override;
    void write_done(const DoneData& data) override;

private:
    std::ostream& out_;
};

class BinOutputWriter final : public OutputWriter {
public:
    explicit BinOutputWriter(std::ostream& out);

    void write_init(const InitData& data) override;
    void write_frame(const FrameData& data) override;
    void write_done(const DoneData& data) override;

private:
    std::ostream& out_;
};

std::unique_ptr<OutputWriter> make_output_writer(OutputFormat format, std::ostream& out);
