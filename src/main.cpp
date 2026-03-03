#include "config/config.hpp"
#include <cassert>

int main(int argc, char** argv) {
    try {
        const Config config{parseArgs(argc, argv)};
        assert(config.numParticles >= 1);
        assert(config.boxLength > 0);
        assert(config.stepSize > 0);
        assert(config.warmupSteps >= 0);
        assert(config.measureSteps >= 1);
        assert(config.blockSize >= 1);
        // WARN: check blocking uncertainty assert
        printConfig(config);
        return 0;
    } catch (const HelpRequested&) {
        printUsage(argv[0]);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Argument parsing error: " << ex.what() << '\n';
        printUsage(argv[0]);
        return 1;
    }
}
