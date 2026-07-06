cmake_minimum_required(VERSION 3.24)

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

option(VMC_ENABLE_CUDA "Enable CUDA acceleration if a CUDA compiler is available" ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- CUDA detection ---
set(VMC_CUDA OFF)
if(VMC_ENABLE_CUDA)
  include(CheckLanguage)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      set(CMAKE_CUDA_ARCHITECTURES native)
    endif()
    set(VMC_CUDA ON)
  endif()
endif()

# --- Floating point precision ---
option(FP_64 "Use double precision (real_t = double) instead of single precision (real_t = float)" ON)
option(VMC_FAST_MATH "Use fast approximate math intrinsics (__sincosf, __expf, ...) in CUDA device code; forces single precision (FP_64=OFF)" OFF)

if(VMC_FAST_MATH AND FP_64)
  message(STATUS "VMC_FAST_MATH is ON: forcing single precision (FP_64=OFF).")
  set(FP_64 OFF)
endif()

# --- Dependencies ---
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
  set(VMC_RELEASE_OPTIONS -O3 -march=native)
  set(VMC_DEBUG_OPTIONS   -g -O0)
  set(VMC_RELWITHDEBINFO_OPT_1 -O3)
  set(VMC_RELWITHDEBINFO_OPT_2 -march=native)
  set(VMC_RELWITHDEBINFO_OPT_3 -fno-omit-frame-pointer)
endif()

if(VMC_CUDA)
  if(MSVC)
    set(VMC_CUDA_XCOMPILE "-Xcompiler=/W4,/permissive-")
    set(VMC_CUDA_RELEASE_OPTIONS -Xcompiler=/O2)
    set(VMC_CUDA_DEBUG_OPTIONS   -Xcompiler=/Od,/Zi -G)
  else()
    set(VMC_CUDA_XCOMPILE "-Xcompiler=-Wall,-Wextra,-march=native")
    set(VMC_CUDA_RELEASE_OPTIONS -O3 -Xcompiler=-O3)
    set(VMC_CUDA_DEBUG_OPTIONS   -g -G -O0 -Xcompiler=-O0)
  endif()
endif()

# --- Source lists ---
set(VMC_CUDA_SOURCES
  src/jastrow_pade/jastrow_pade.cu
  src/wavefunction/wavefunction.cu
  src/slater_plane_wave/log_abs_det.cu
  src/slater_plane_wave/slater_plane_wave.cu
  src/slater_plane_wave/update_restore_trig.cu
  src/simulation/simulation.cu
  src/energy_tracking/energy_tracking.cu
  src/optimizer/jastrow_optimizer.cu
)

set(VMC_CXX_SOURCES
  src/blocking_analysis/blocking_analysis.cpp
  src/output_writer/output_writer.cpp
)

set(VMC_CUDA_EXTENSION_SOURCES
  ${VMC_CUDA_SOURCES}
)

# --- Library target for all source files ---
add_library(vmc_lib STATIC
  ${VMC_CXX_SOURCES}
  ${VMC_CUDA_SOURCES}
)

target_include_directories(vmc_lib PUBLIC src)

if(NOT VMC_CUDA AND VMC_CUDA_EXTENSION_SOURCES)
  set_source_files_properties(
    ${VMC_CUDA_EXTENSION_SOURCES}
    PROPERTIES LANGUAGE CXX
  )

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set_source_files_properties(
      ${VMC_CUDA_EXTENSION_SOURCES}
      PROPERTIES COMPILE_OPTIONS "-x;c++"
    )
  endif()
endif()

if(FP_64)
  target_compile_definitions(vmc_lib PUBLIC FP_64)
endif()

if(VMC_FAST_MATH)
  target_compile_definitions(vmc_lib PUBLIC VMC_FAST_MATH)
endif()

if(OpenMP_CXX_FOUND)
  target_link_libraries(vmc_lib PUBLIC OpenMP::OpenMP_CXX)
  if(MSVC)
    target_compile_options(vmc_lib PUBLIC $<$<COMPILE_LANGUAGE:CXX>:/openmp:experimental>)
  endif()
elseif(VMC_HAS_FOPENMP_SIMD)
  target_compile_options(vmc_lib PUBLIC -fopenmp-simd)
endif()
target_compile_options(vmc_lib PRIVATE
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

if(VMC_CUDA)
  find_package(CUDAToolkit REQUIRED)

  set_target_properties(vmc_lib PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
  )

  target_compile_definitions(vmc_lib PUBLIC VMC_CUDA_BACKEND)

  target_link_libraries(vmc_lib PUBLIC
    CUDA::cudart
    CUDA::cusolver
    CUDA::cublas
  )
endif()

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