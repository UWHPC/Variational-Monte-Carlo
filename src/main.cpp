#include "config/config.hpp"
#include "simulation/simulation.hpp"

int main(int argc, char** argv) {
    try {
        const Config config{parse_args(argc, argv)};
        Simulation sim{config};
        sim.run();
        return 0;
    } catch (const HelpRequested&) {
        // TODO: make the help section better (and different than usage)
        print_usage(argv[0]);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Argument parsing error: " << ex.what() << '\n';
        print_usage(argv[0]);
        return 1;
    }
}
