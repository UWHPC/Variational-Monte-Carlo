cmake_minimum_required(VERSION 3.20)

if(
    NOT DEFINED CMAKE_EXPORT_COMPILE_COMMANDS
    OR CMAKE_EXPORT_COMPILE_COMMANDS STREQUAL ""
)
    set(
        CMAKE_EXPORT_COMPILE_COMMANDS
        ON
        CACHE BOOL "Generate compile_commands.json for tooling" FORCE
    )
endif()

project(Variational-Monte-Carlo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Dependencies ---
find_package(OpenMP REQUIRED COMPONENTS CXX)

# --- Build type default (Release unless overridden) ---
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)

# --- Compiler warnings and optimizations ---
if(MSVC)
    set(VMC_WARNING_OPTIONS /W4 /permissive- /arch:AVX2)
    set(VMC_RELEASE_OPTIONS /O2)
    set(VMC_DEBUG_OPTIONS   /Od /Zi)
    set(VMC_RELWITHDEBINFO_OPT_1 /O2)
    set(VMC_RELWITHDEBINFO_OPT_2 /Oy-)
    set(VMC_RELWITHDEBINFO_OPT_3 "")
else()
    set(VMC_WARNING_OPTIONS
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wconversion
        -Wsign-conversion
    )
    set(VMC_RELEASE_OPTIONS -O3 -march=native -ffast-math -fno-finite-math-only)
    set(VMC_DEBUG_OPTIONS   -g -O0)
    set(VMC_RELWITHDEBINFO_OPT_1 -O3)
    set(VMC_RELWITHDEBINFO_OPT_2 -march=native)
    set(VMC_RELWITHDEBINFO_OPT_3 -fno-omit-frame-pointer -ffast-math)
endif()

# --- Library target for all source files ---
add_library(vmc_lib STATIC
    src/jastrow_pade/jastrow_pade.cpp
    src/wavefunction/wavefunction.cpp
    src/slater_plane_wave/slater_plane_wave.cpp
    src/simulation/simulation.cpp
    src/blocking_analysis/blocking_analysis.cpp
    src/output_writer/output_writer.cpp
    src/energy_tracking/energy_tracking.cpp
)

target_include_directories(vmc_lib PUBLIC src)
target_link_libraries(vmc_lib PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(vmc_lib PRIVATE
    ${VMC_WARNING_OPTIONS}
    $<$<CONFIG:Release>:${VMC_RELEASE_OPTIONS}>
    $<$<CONFIG:Debug>:${VMC_DEBUG_OPTIONS}>
    $<$<CONFIG:RelWithDebInfo>:${VMC_RELWITHDEBINFO_OPT_1}>
    $<$<CONFIG:RelWithDebInfo>:${VMC_RELWITHDEBINFO_OPT_2}>
    $<$<CONFIG:RelWithDebInfo>:${VMC_RELWITHDEBINFO_OPT_3}>
)

# --- Main executable ---
add_executable(vmc
    src/main.cpp
)

target_link_libraries(vmc PRIVATE vmc_lib)
target_compile_options(vmc PRIVATE
    ${VMC_WARNING_OPTIONS}
    $<$<CONFIG:Release>:${VMC_RELEASE_OPTIONS}>
    $<$<CONFIG:Debug>:${VMC_DEBUG_OPTIONS}>
    $<$<CONFIG:RelWithDebInfo>:${VMC_RELWITHDEBINFO_OPT_1}>
    $<$<CONFIG:RelWithDebInfo>:${VMC_RELWITHDEBINFO_OPT_2}>
    $<$<CONFIG:RelWithDebInfo>:${VMC_RELWITHDEBINFO_OPT_3}>
)

# --- Testing (opt-in via BUILD_TESTING argument) ---
include(CTest)

if(BUILD_TESTING)

    include(FetchContent)

    fetchcontent_declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.8.1
    )

    fetchcontent_makeavailable(Catch2)

    add_executable(vmc_tests
        tests/test_tests.cpp
    )

    target_link_libraries(vmc_tests PRIVATE
        Catch2::Catch2WithMain
        vmc_lib
    )
    target_compile_options(vmc_tests PRIVATE
        ${VMC_WARNING_OPTIONS}
        $<$<CONFIG:Release>:${VMC_RELEASE_OPTIONS}>
        $<$<CONFIG:Debug>:${VMC_DEBUG_OPTIONS}>
        $<$<CONFIG:RelWithDebInfo>:${VMC_RELWITHDEBINFO_OPT_1}>
        $<$<CONFIG:RelWithDebInfo>:${VMC_RELWITHDEBINFO_OPT_2}>
        $<$<CONFIG:RelWithDebInfo>:${VMC_RELWITHDEBINFO_OPT_3}>
    )

    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    include(Catch)
    catch_discover_tests(vmc_tests)

endif()