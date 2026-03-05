#include "config/config.hpp"

int main(int argc, char** argv) {
    try {
        const Config config{parseArgs(argc, argv)};
        print_config(config);
        return 0;
    } catch (const HelpRequested&) {
        print_usage(argv[0]);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Argument parsing error: " << ex.what() << '\n';
        print_usage(argv[0]);
        return 1;
    }
}
