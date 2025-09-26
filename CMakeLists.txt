# Let's start by requiring a reasonably modern version of CMake.
cmake_minimum_required(VERSION 3.2)

# This is us! Our project's name is TracEon.
project(TracEon)

# We agreed on C++20 for this project, so let's make sure it's enforced.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Here's the heart of our project: the 'traceon_core' static library.
# All our main C++ source files will get added here eventually.
add_library(traceon_core STATIC
        Cache.cpp  # <-- ADD THIS LINE
        # e.g., Cache.cpp will go here
)
# Make sure our library knows where to find its own header files.
target_include_directories(traceon_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# This is just a little app for any manual testing or demos we want to do.
add_executable(trac-eon-test main.cpp
        main.cpp
        Cache.cpp
        Cache.h
        IEncodingStrategy.h
        RecordTypes.h)
target_link_libraries(trac-eon-test PRIVATE traceon_core)

# --- Test Setup ---

# This is a neat CMake trick to automatically download Catch2 from GitHub for us.
# No need to manually add the library to our project, which is nice.
include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.10.0 # Pinning to a specific version for stability
)
FetchContent_MakeAvailable(Catch2)

# Now let's set up everything for our unit tests.
enable_testing()

# We'll create a separate program just for running our tests...
add_executable(unit_tests
        tests/main.cpp
        tests/CacheTests.cpp
        # Our test files, like tests/CacheTests.cpp, will be added here
)

# ...and of course, our test program needs to be linked against our library and Catch2.
target_link_libraries(unit_tests PRIVATE traceon_core Catch2::Catch2WithMain)

# This officially registers our 'unit_tests' program with CTest, CMake's test runner.
include(CTest)
add_test(NAME TracEonTests COMMAND unit_tests)