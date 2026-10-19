cmake_minimum_required(VERSION 3.20)
project(MonteCarloPricer VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # lets editors / clangd understand the project

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Core library
add_library(mcpricer
  src/black_scholes.cpp
  src/monte_carlo.cpp
  src/implied_vol.cpp
)
target_include_directories(mcpricer PUBLIC include)
target_compile_options(mcpricer PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
)

# Demo executable
add_executable(mcpricer_demo src/main.cpp)
target_link_libraries(mcpricer_demo PRIVATE mcpricer)

# Tests
enable_testing()
add_executable(mcpricer_tests tests/test_pricer.cpp)
target_link_libraries(mcpricer_tests PRIVATE mcpricer)
add_test(NAME mcpricer_tests COMMAND mcpricer_tests)
