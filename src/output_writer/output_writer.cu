#include "output_writer.hpp"

#include <memory>
#include <stdexcept>

CsvOutputWriter::CsvOutputWriter(std::ostream& out)
: out_{out}
{ }

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

BinOutputWriter::BinOutputWriter(std::ostream& out)
: out_{out}
{ }

void BinOutputWriter::write_init(const InitData& data) {
  const uint64_t np{static_cast<uint64_t>(data.num_particles)};
  const fp_t bl{data.box_length};
  const uint64_t ms{static_cast<uint64_t>(data.measure_steps)};
  out_.write(reinterpret_cast<const char*>(&np), sizeof(np));
  out_.write(reinterpret_cast<const char*>(&bl), sizeof(bl));
  out_.write(reinterpret_cast<const char*>(&ms), sizeof(ms));
  out_.flush();
}

void BinOutputWriter::write_frame(const FrameData& data) {
  const fp_t se{data.standard_error.value_or(0.0_fp)};
  out_.write(reinterpret_cast<const char*>(&data.local_energy), sizeof(fp_t));
  out_.write(reinterpret_cast<const char*>(&data.mean_energy), sizeof(fp_t));
  out_.write(reinterpret_cast<const char*>(&se), sizeof(fp_t));
  out_.write(reinterpret_cast<const char*>(&data.acceptance_rate), sizeof(fp_t));
  out_.write(
    reinterpret_cast<const char*>(data.positions.data()),
    static_cast<std::streamsize>(data.positions.size() * sizeof(fp_t))
  );
}

void BinOutputWriter::write_done(const DoneData&) { out_.flush(); }

std::unique_ptr<OutputWriter> make_output_writer(OutputFormat format, std::ostream& out) {
  switch (format) {
  case OutputFormat::CSV:
    return std::make_unique<CsvOutputWriter>(out);
  case OutputFormat::BIN:
    return std::make_unique<BinOutputWriter>(out);
  }
  throw std::invalid_argument("Unsupported output format");
}
