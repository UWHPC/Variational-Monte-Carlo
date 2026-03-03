#include "config/config.hpp"
#include "particles/particles.hpp"
#include "pbc/pbc.hpp"
#include "blocking_analysis/blocking_analysis.hpp"
#include "simulation/simulation.hpp"
#include "wavefunction/wavefunction.hpp"

int main(int argc, char** argv) {
    try {
        const Config config{parseArgs(argc, argv)};
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