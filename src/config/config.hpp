#pragma once

#include <cstddef>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

// // Closed shells (up to 10k):
// N = 1, 7, 19, 27, 33, 57, 81, 93, 123, 147, 171, 179,
// 203, 251, 257, 305, 341, 365, 389, 437, 461, 485, 515, 587,
// 619, 691, 739, 751, 799, 847, 895, 925, 949, 1021, 1045, 1141,
// 1189, 1213, 1237, 1309, 1357, 1365, 1419, 1503, 1551, 1575, 1647, 1743,
// 1791, 1839, 1863, 1935, 2007, 2103, 2109, 2205, 2301, 2325, 2373, 2469,
// 2517, 2553, 2601, 2721, 2777, 2801, 2897, 2945, 2969, 3071, 3119, 3191,
// 3239, 3287, 3407, 3431, 3575, 3695, 3743, 3791, 3887, 3911, 3959, 4067,
// 4139, 4169, 4337, 4385, 4457, 4553, 4625, 4697, 4729, 4801, 4945, 5041,
// 5137, 5185, 5257, 5377, 5449, 5497, 5575, 5695, 5743, 5887, 6031, 6043,
// 6187, 6235, 6355, 6403, 6451, 6619, 6667, 6763, 6859, 6931, 6979, 7075,
// 7123, 7153, 7249, 7441, 7497, 7521, 7689, 7809, 7881, 8025, 8121, 8217,
// 8289, 8385, 8409, 8601, 8709, 8733, 8829, 8925, 9045, 9093, 9171, 9315,
// 9435, 9459, 9627, 9771, 9795, 9843, 9939, 10059

class Config {
public:
    // Independent parameters defaults; overridable in config.cfg:
    std::size_t num_threads{1}; // Number of threads used
    std::size_t num_particles{7U}; // Number of particles
    std::size_t warmup_sweeps{100U}; // Number of sweeps for warmup
    std::size_t measure_sweeps{100U}; // Number of sweeps for measuring
    double box_length{10.0};         // Length of box (grid)
    std::size_t block_size{100U};    // Size of block
    uint64_t master_seed{42U};      // Random seed
    bool is_master_thread{};   // Is the Master Thread

    // Derived (computed from above):
    std::size_t warmup_steps;  // Warmup steps
    std::size_t measure_steps; // Measure steps
    double step_size;          // proposal half-width s

    // Default constructor uses hardcoded defaults:
    Config() { compute_derived(); }

    // Load from file, falling back to defaults for missing keys:
    [[nodiscard]] static Config from_file( std::string const &path ) {
        Config cfg{};

        std::ifstream file{ path };
        if ( !file.is_open() ) {
            std::cout << "[Config] File not found: " << path
                      << "; using defaults.\n";
            return cfg;
        }

        auto map{ parse( file ) };

        read( map, "Num_Threads", cfg.num_threads );
        read( map, "Num_Particles", cfg.num_particles );
        read( map, "Warmup_Sweeps", cfg.warmup_sweeps );
        read( map, "Measure_Sweeps", cfg.measure_sweeps );
        read( map, "Box_Length", cfg.box_length );
        read( map, "Block_Size", cfg.block_size );
        read( map, "Master_Seed", cfg.master_seed );
        read( map, "Is_Master_Thread", cfg.is_master_thread );

        if ( !map.empty() ) {
            std::cout << "[Config] Warning! unknown keys ignored:";
            for ( auto const &[key, val] : map ) {
                std::cout << " " << key;
            }
            std::cout << "\n";
        }

        cfg.compute_derived();
        return cfg;
    }

private:
    void compute_derived() {
        this->warmup_steps = num_particles * warmup_sweeps;
        this->measure_steps = num_particles * measure_sweeps;
        this->step_size = box_length / 10.0;
    }

    // Because C++ is so ugly:
    using Map = std::unordered_map<std::string, std::string>;

    // Parser:
    [[nodiscard]] static Map parse( std::ifstream &file ) {
        Map map{};
        std::string line{};

        while ( std::getline( file, line ) ) {
            // Strip comments:
            if ( auto pos{ line.find( '#' ) }; pos != std::string::npos ) {
                line.erase( pos );
            }
            // Skip blank lines:
            if ( line.find_first_not_of( " \t\r" ) == std::string::npos ) {
                continue;
            }
            auto eq{ line.find( '=' ) };
            if ( eq == std::string::npos ) {
                continue;
            }
            std::string key{ trim( line.substr( 0, eq ) ) };
            std::string val{ trim( line.substr( eq + 1 ) ) };

            if ( !key.empty() && !val.empty() ) {
                map[key] = val;
            }
        }
        return map;
    }

    [[nodiscard]] static std::string trim( std::string s ) {
        auto const start{ s.find_first_not_of( " \t\r" ) };
        if ( start == std::string::npos ) return {};
        auto const end{ s.find_last_not_of( " \t\r" ) };
        
        return s.substr( start, end - start + 1 );
    }

    // Type-specific readers; erase consumed keys so leftovers can be warned:

    static void read( Map &map, std::string const &key, double &out ) {
        if ( auto it{ map.find( key ) }; it != map.end() ) {
            out = std::stod( it->second );
            map.erase( it );
        }
    }

    static void read( Map &map, std::string const &key, std::size_t &out ) {
        if ( auto it{ map.find( key ) }; it != map.end() ) {
            out = static_cast<std::size_t>( std::stoull( it->second ) );
            map.erase( it );
        }
    }

    #if UINTPTR_MAX != UINT64_MAX
    static void read( Map &map, std::string const &key, uint64_t &out ) {
        if ( auto it{ map.find( key ) }; it != map.end() ) {
            out = std::stoull( it->second );
            map.erase( it );
        }
    }
    #endif

    static void read( Map &map, std::string const &key, int &out ) {
        if ( auto it{ map.find( key ) }; it != map.end() ) {
            out = std::stoi( it->second );
            map.erase( it );
        }
    }

    static void read( Map &map, std::string const &key, bool &out ) {
        if ( auto it{ map.find( key ) }; it != map.end() ) {
            std::string val{ it->second };
            // Lowercase for comparison:
            for ( auto &ch : val ) ch = static_cast<char>( std::tolower( ch ) );
            out = ( val == "true" || val == "1" );
            map.erase( it );
        }
    }
};