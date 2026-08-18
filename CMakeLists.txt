cmake_minimum_required(VERSION 3.25)

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
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Toolchain requirements ---
if(WIN32)
  message(FATAL_ERROR
    "Native Windows builds are not supported. Build under WSL2 or Linux.")
endif()

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  message(FATAL_ERROR
    "GCC is the only supported host compiler (found ${CMAKE_CXX_COMPILER_ID}). "
    "Configure with -DCMAKE_CXX_COMPILER=g++-15.")
endif()

if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15)
  message(FATAL_ERROR
    "GCC 15 or newer is required (found ${CMAKE_CXX_COMPILER_VERSION}). "
    "Configure with -DCMAKE_CXX_COMPILER=g++-15.")
endif()

option(VMC_ENABLE_CUDA "Enable CUDA acceleration if a CUDA compiler is available" ON)

# --- CUDA detection ---
set(VMC_CUDA OFF)
if(VMC_ENABLE_CUDA)
  include(CheckLanguage)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    if(NOT CMAKE_CUDA_HOST_COMPILER)
      set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
    endif()

    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 23)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_EXTENSIONS OFF)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      set(CMAKE_CUDA_ARCHITECTURES native)
    endif()
    set(VMC_CUDA ON)
    message(STATUS
      "VMC: CUDA ${CMAKE_CUDA_COMPILER_VERSION}, arch ${CMAKE_CUDA_ARCHITECTURES}, "
      "device host compiler ${CMAKE_CUDA_HOST_COMPILER}")
  endif()
endif()

# --- Floating point precision ---
option(FP_64 "Use double precision (real_t = double) instead of single precision (real_t = float)" ON)
option(VMC_FAST_MATH "Fast approximate math intrinsics in device code (inert: pending an xpu fast-math tier); forces single precision (FP_64=OFF)" OFF)

if(VMC_FAST_MATH AND FP_64)
  message(STATUS "VMC_FAST_MATH is ON: forcing single precision (FP_64=OFF).")
  set(FP_64 OFF)
endif()

# --- Dependencies ---
include(FetchContent)

set(XPU_ENABLE_CUDA ${VMC_CUDA} CACHE BOOL "" FORCE)
set(XPU_ENABLE_LINALG ON CACHE BOOL "" FORCE)
set(XPU_ARCH_NATIVE OFF CACHE BOOL "" FORCE)

fetchcontent_declare(
  xpu
  GIT_REPOSITORY https://github.com/UWHPC/xpu.git
  GIT_TAG aa950a6
)
fetchcontent_makeavailable(xpu)

find_package(OpenMP QUIET COMPONENTS CXX)

if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found (CXX ${OpenMP_CXX_VERSION}); using it for '#pragma omp simd'.")
else()
  message(STATUS "OpenMP not found; building without it (we only use '#pragma omp simd').")
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-fopenmp-simd" VMC_HAS_FOPENMP_SIMD)
endif()

# --- Build type default (Release unless overridden) ---
if(NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# --- Compiler warnings and optimizations ---
set(VMC_WARNING_OPTIONS
  -Wall
  -Wextra
  -Wpedantic
  -Wshadow
  -Wconversion
  -Wsign-conversion
)
set(VMC_RELEASE_OPTIONS -O3 -march=native)
set(VMC_DEBUG_OPTIONS   -g -O0)
set(VMC_RELWITHDEBINFO_OPT_1 -O3)
set(VMC_RELWITHDEBINFO_OPT_2 -march=native)
set(VMC_RELWITHDEBINFO_OPT_3 -fno-omit-frame-pointer)

if(VMC_CUDA)
  set(VMC_CUDA_XCOMPILE "-Xcompiler=-Wall,-Wextra,-march=native")
  set(VMC_CUDA_RELEASE_OPTIONS -O3 -Xcompiler=-O3)
  set(VMC_CUDA_DEBUG_OPTIONS   -g -G -O0 -Xcompiler=-O0)
endif()

# --- Source lists ---
set(VMC_CUDA_SOURCES
  src/jastrow_pade/jastrow_pade.cu
  src/wavefunction/wavefunction.cu
  src/slater_plane_wave/slater_plane_wave.cu
  src/simulation/simulation.cu
  src/energy_tracking/energy_tracking.cu
  src/optimizer/jastrow_optimizer.cu
  src/blocking_analysis/blocking_analysis.cu
  src/output_writer/output_writer.cu
)

set(VMC_MAIN_SOURCE src/main.cu)
set(VMC_TEST_SOURCE tests/test_tests.cu)

set(VMC_CUDA_EXTENSION_SOURCES
  ${VMC_CUDA_SOURCES}
  ${VMC_MAIN_SOURCE}
  ${VMC_TEST_SOURCE}
)

# --- Library target for all source files ---
add_library(vmc_lib STATIC
  ${VMC_CUDA_SOURCES}
)

target_include_directories(vmc_lib PUBLIC src)

target_link_libraries(vmc_lib PUBLIC
  xpu::xpu
  xpu::linalg
)

if(NOT VMC_CUDA AND VMC_CUDA_EXTENSION_SOURCES)
  set_source_files_properties(
    ${VMC_CUDA_EXTENSION_SOURCES}
    PROPERTIES LANGUAGE CXX
  )

  set_source_files_properties(
    ${VMC_CUDA_EXTENSION_SOURCES}
    PROPERTIES COMPILE_OPTIONS "-x;c++"
  )
endif()

if(FP_64)
  target_compile_definitions(vmc_lib PUBLIC FP_64)
endif()

if(VMC_FAST_MATH)
  target_compile_definitions(vmc_lib PUBLIC VMC_FAST_MATH)
endif()

if(OpenMP_CXX_FOUND)
  target_link_libraries(vmc_lib PUBLIC OpenMP::OpenMP_CXX)
elseif(VMC_HAS_FOPENMP_SIMD)
  target_compile_options(vmc_lib PUBLIC -fopenmp-simd)
endif()
set(VMC_COMPILE_OPTIONS
  $<$<COMPILE_LANGUAGE:CXX>:${VMC_WARNING_OPTIONS}>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:${VMC_RELEASE_OPTIONS}>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>:${VMC_DEBUG_OPTIONS}>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RelWithDebInfo>>:${VMC_RELWITHDEBINFO_OPT_1}>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RelWithDebInfo>>:${VMC_RELWITHDEBINFO_OPT_2}>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RelWithDebInfo>>:${VMC_RELWITHDEBINFO_OPT_3}>
  $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<BOOL:${VMC_CUDA}>>:${VMC_CUDA_XCOMPILE}>
  $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<BOOL:${VMC_CUDA}>,$<CONFIG:Release>>:${VMC_CUDA_RELEASE_OPTIONS}>
  $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<BOOL:${VMC_CUDA}>,$<CONFIG:Debug>>:${VMC_CUDA_DEBUG_OPTIONS}>
)

target_compile_options(vmc_lib PRIVATE ${VMC_COMPILE_OPTIONS})

if(VMC_CUDA)
  set_target_properties(vmc_lib PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
  )
endif()

# --- Main executable ---
add_executable(vmc
  ${VMC_MAIN_SOURCE}
)

target_link_libraries(vmc PRIVATE vmc_lib)
target_compile_options(vmc PRIVATE ${VMC_COMPILE_OPTIONS})

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
    ${VMC_TEST_SOURCE}
  )

  target_link_libraries(vmc_tests PRIVATE
    Catch2::Catch2WithMain
    vmc_lib
  )
  target_compile_options(vmc_tests PRIVATE ${VMC_COMPILE_OPTIONS})

  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
  include(Catch)
  catch_discover_tests(vmc_tests)

endif()
