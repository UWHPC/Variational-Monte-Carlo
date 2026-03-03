cmake_minimum_required(VERSION 3.20)
project(Variational-Monte-Carlo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Compiler warnings and optimizations ---
if(MSVC)
    add_compile_options(/W4 /permissive- /arch:AVX2)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion)
endif()

# --- Library target for all source files ---
add_library(vmc_lib STATIC
  src/pbc/pbc.cpp
  src/jastrow_pade/jastrow_pade.cpp
  src/wavefunction/wavefunction.cpp
  src/slater_plane_wave/slater_plane_wave.cpp
  src/simulation/simulation.cpp
)
target_include_directories(vmc_lib PUBLIC src)

# --- Main executable ---
add_executable(vmc src/main.cpp)
target_link_libraries(vmc PRIVATE vmc_lib)

# --- Testing ---
#include(CTest)
#enable_testing()

# include(FetchContent)

# fetchcontent_declare(
#   Catch2
#   GIT_REPOSITORY https://github.com/catchorg/Catch2.git
#   GIT_TAG        v3.8.1
# )
# fetchcontent_makeavailable(Catch2)

# add_executable(vmc_tests
#   tests/test_tests.cpp
# )
# target_link_libraries(vmc_tests PRIVATE Catch2::Catch2WithMain vmc_lib)

# include(Catch)
# catch_discover_tests(vmc_tests)
