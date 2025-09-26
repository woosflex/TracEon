# Sticking to a modern, stable CMake version.
cmake_minimum_required(VERSION 3.20)

# This is us! Our project's name is TracEon.
project(TracEon)

# We agreed on C++20 for this project, so let's make sure it's enforced.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Here's our main library, 'traceon_core'. All our core logic lives here.
add_library(traceon_core STATIC
        Cache.cpp
        PlainTextStrategy.cpp
)
# Make sure our library knows where to find its own header files.
target_include_directories(traceon_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# This is our little app for manual testing. It only needs main.cpp,
# because all the other code is linked in from the traceon_core library.
add_executable(trac-eon-test main.cpp)
target_link_libraries(trac-eon-test PRIVATE traceon_core)

# --- Test Setup ---

# This is a neat CMake trick to grab Catch2 for us directly from GitHub.
include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.6.0 # Let's stick with this tested version for now
)
FetchContent_MakeAvailable(Catch2)

# Now let's set up everything for our unit tests.
enable_testing()

# We'll create a separate program just for running our tests.
add_executable(unit_tests
        tests/main.cpp
        tests/CacheTests.cpp
)

# ...and our test program needs to be linked against our library and Catch2.
target_link_libraries(unit_tests PRIVATE traceon_core Catch2::Catch2WithMain)

# This officially registers our 'unit_tests' program with CTest.
include(CTest)
add_test(NAME TracEonTests COMMAND unit_tests)