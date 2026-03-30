#include "output_writer.hpp"

#include <limits>
#include <memory>
#include <stdexcept>
#include <string_view>

namespace {

void write_json_string(std::ostream& out, std::string_view value) {
    constexpr char kHexDigits[]{"0123456789abcdef"};

    out << '"';
    for (const char raw_ch : value) {
        const unsigned char ch{static_cast<unsigned char>(raw_ch)};
        switch (ch) {
        case '\"':
            out << "\\\"";
            break;
        case '\\':
            out << "\\\\";
            break;
        case '\b':
            out << "\\b";
            break;
        case '\f':
            out << "\\f";
            break;
        case '\n':
            out << "\\n";
            break;
        case '\r':
            out << "\\r";
            break;
        case '\t':
            out << "\\t";
            break;
        default:
            if (ch < 0x20U) {
                const std::size_t hi{static_cast<std::size_t>((ch >> 4U) & 0x0FU)};
                const std::size_t lo{static_cast<std::size_t>(ch & 0x0FU)};
                out << "\\u00" << kHexDigits[hi] << kHexDigits[lo];
            } else {
                out << raw_ch;
            }
            break;
        }
    }
    out << '"';
}

void write_positions(std::ostream& out, const std::vector<double>& positions) {
    out << '[';
    for (std::size_t i{}; i < positions.size(); ++i) {
        if (i != 0U) {
            out << ',';
        }
        out << positions[i];
    }
    out << ']';
}

} // namespace

JsonOutputWriter::JsonOutputWriter(std::ostream& out) : out_{out} {
    out_.precision(std::numeric_limits<double>::max_digits10);
}

void JsonOutputWriter::write_init(const InitData& data) {
    out_ << "{\"type\":\"init\",\"runId\":";
    write_json_string(out_, data.run_id);
    out_ << ",\"numParticles\":" << data.num_particles << ",\"boxLength\":" << data.box_length
         << ",\"warmupSteps\":" << data.warmup_steps << ",\"measureSteps\":" << data.measure_steps
         << ",\"stepSize\":" << data.step_size << ",\"blockSize\":" << data.block_size
         << ",\"seed\":" << data.seed << "}\n";
}

void JsonOutputWriter::write_frame(const FrameData& data) {
    out_ << "{\"type\":\"frame\",\"step\":" << data.step << ",\"accepted\":" << data.accepted
         << ",\"proposed\":" << data.proposed << ",\"acceptanceRate\":" << data.acceptance_rate
         << ",\"localEnergy\":" << data.local_energy << ",\"meanEnergy\":" << data.mean_energy;

    if (data.standard_error.has_value()) {
        out_ << ",\"standardErrorAvailable\":true,\"standardError\":" << *data.standard_error;
    } else {
        out_ << ",\"standardErrorAvailable\":false,\"standardError\":0";
    }

    out_ << ",\"positions\":";
    write_positions(out_, data.positions);
    out_ << "}\n";
}

void JsonOutputWriter::write_done(const DoneData& data) {
    out_ << "{\"type\":\"done\",\"totalAccepted\":" << data.total_accepted
         << ",\"totalProposed\":" << data.total_proposed
         << ",\"finalAcceptanceRate\":" << data.final_acceptance_rate
         << ",\"finalMeanEnergy\":" << data.final_mean_energy;

    if (data.final_standard_error.has_value()) {
        out_ << ",\"finalStandardErrorAvailable\":true,\"finalStandardError\":"
             << *data.final_standard_error;
    } else {
        out_ << ",\"finalStandardErrorAvailable\":false,\"finalStandardError\":0";
    }

    out_ << "}\n";
}

CsvOutputWriter::CsvOutputWriter(std::ostream& out) : out_{out} {}

void CsvOutputWriter::write_init(const InitData&) {
    static_cast<void>(out_);
    throw std::logic_error{"CsvOutputWriter is not implemented yet"};
}

void CsvOutputWriter::write_frame(const FrameData&) {
    static_cast<void>(out_);
    throw std::logic_error{"CsvOutputWriter is not implemented yet"};
}

void CsvOutputWriter::write_done(const DoneData&) {
    static_cast<void>(out_);
    throw std::logic_error{"CsvOutputWriter is not implemented yet"};
}

BinOutputWriter::BinOutputWriter(std::ostream& out) : out_{out} {}

void BinOutputWriter::write_init(const InitData& data) {
    // Header: num_particles (u64), box_length (f64), measure_steps (u64)
    const uint64_t np{static_cast<uint64_t>(data.num_particles)};
    const double bl{data.box_length};
    const uint64_t ms{static_cast<uint64_t>(data.measure_steps)};
    out_.write(reinterpret_cast<const char*>(&np), sizeof(np));
    out_.write(reinterpret_cast<const char*>(&bl), sizeof(bl));
    out_.write(reinterpret_cast<const char*>(&ms), sizeof(ms));
    out_.flush();
}

void BinOutputWriter::write_frame(const FrameData& data) {
    // Per frame: local_energy, mean_energy, standard_error, acceptance_rate, positions[3*N]
    const double se{data.standard_error.value_or(0.0)};
    out_.write(reinterpret_cast<const char*>(&data.local_energy), sizeof(double));
    out_.write(reinterpret_cast<const char*>(&data.mean_energy), sizeof(double));
    out_.write(reinterpret_cast<const char*>(&se), sizeof(double));
    out_.write(reinterpret_cast<const char*>(&data.acceptance_rate), sizeof(double));
    out_.write(reinterpret_cast<const char*>(data.positions.data()),
               static_cast<std::streamsize>(data.positions.size() * sizeof(double)));
}

void BinOutputWriter::write_done(const DoneData&) {
    out_.flush();
}

std::unique_ptr<OutputWriter> make_output_writer(OutputFormat format, std::ostream& out) {
    switch (format) {
    case OutputFormat::JSON:
        return std::make_unique<JsonOutputWriter>(out);
    case OutputFormat::CSV:
        return std::make_unique<CsvOutputWriter>(out);
    case OutputFormat::BIN:
        return std::make_unique<BinOutputWriter>(out);
    }
    throw std::invalid_argument("Unsupported output format");
}
