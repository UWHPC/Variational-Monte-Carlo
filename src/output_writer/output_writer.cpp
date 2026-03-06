#include "output_writer.hpp"

#include <memory>
#include <stdexcept>

std::unique_ptr<OutputWriter> make_output_writer(OutputFormat format, std::ostream& out) {
    switch (format) {
    case OutputFormat::JSON:
        return std::make_unique<JsonOutputWriter>(out);
    case OutputFormat::CSV:
        return std::make_unique<CsvOutputWriter>(out);
    }
    throw std::invalid_argument("Unsupported output format");
}
