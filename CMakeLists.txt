# Sticking to a modern, stable CMake version.
cmake_minimum_required(VERSION 3.20)

# This is us! Our project's name is TracEon.
project(TracEon)

# We're using C++20 for all the nice modern features.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Source File Definitions ---
# By defining our source files in variables, we only have to list them once.
set(CORE_SOURCES
        Cache.cpp
        IEncodingStrategy.cpp
        SmartStrategy.cpp
)

set(TEST_SOURCES
        tests/main.cpp
        tests/CacheTests.cpp
        tests/FastqTests.cpp
        tests/SmartStrategyTests.cpp
)

# --- Library Target ---
# Here's our main library, 'traceon_core'. All our core logic lives here.
add_library(traceon_core STATIC ${CORE_SOURCES})
target_include_directories(traceon_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# --- Executable Targets ---
# A simple app for manual testing.
add_executable(trac-eon-test main.cpp)
target_link_libraries(trac-eon-test PRIVATE traceon_core)

# Benchmark executables.
add_executable(benchmark_disk benchmark_disk.cpp)
add_executable(benchmark_traceon benchmark_traceon.cpp)
target_link_libraries(benchmark_traceon PRIVATE traceon_core)

# --- Test Setup using Catch2 ---
# This CMake trick automatically downloads Catch2 from GitHub for us.
include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.6.0 # Pinning to a specific version for stability
)
FetchContent_MakeAvailable(Catch2)

# Now let's set up everything for our unit tests.
enable_testing()

# We'll create a separate program just for running our tests...
add_executable(unit_tests ${TEST_SOURCES})

# ...and our test program needs to be linked against our library and Catch2.
target_link_libraries(unit_tests PRIVATE traceon_core Catch2::Catch2WithMain)

# This officially registers our 'unit_tests' program with CTest.
include(CTest)
add_test(NAME TracEonTests COMMAND unit_tests)