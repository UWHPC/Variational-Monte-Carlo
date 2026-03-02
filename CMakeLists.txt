cmake_minimum_required(VERSION 3.20)
project(Variational-Monte-Carlo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)
enable_testing()

include(FetchContent)

fetchcontent_declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.8.1
)
fetchcontent_makeavailable(Catch2)

add_executable(vmc_tests
  tests/test_tests.cpp
)
target_link_libraries(vmc_tests PRIVATE Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(vmc_tests)


