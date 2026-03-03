# Sticking to a modern, stable CMake version.
cmake_minimum_required(VERSION 3.20)

# This is us! Our project's name is TracEon.
project(TracEon)

# We're using C++20 for all the nice modern features.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- ZLIB Support ---
add_subdirectory(third_party/zlib)

# --- LZ4 Support ---
add_subdirectory(third_party/lz4/build/cmake)

# --- Source File Definitions (UPDATED PATHS) ---
set(CORE_SOURCES
        src/Cache.cpp
        src/IEncodingStrategy.cpp
        src/SmartStrategy.cpp
        # MapDefs.h is a header, not a source, so strictly speaking not needed here
        # but kept if you want IDEs to index it easily via add_library
)

set(TEST_SOURCES
        tests/main.cpp
        tests/CacheTests.cpp
        tests/FastqTests.cpp
        tests/SmartStrategyTests.cpp
)

# --- Library Target ---
add_library(traceon_core STATIC ${CORE_SOURCES})

# Include directories (UPDATED PATHS)
# traceon_core users need 'include/'
# traceon_core implementation needs 'third_party/lz4/lib' and 'third_party/'
target_include_directories(traceon_core PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lz4/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/zlib
    ${CMAKE_CURRENT_SOURCE_DIR} # For third_party/robin_hood.h inclusion
)

# Link traceon_core to both zlib and lz4
target_link_libraries(traceon_core PUBLIC zlib_static lz4_static)

# ── Architecture-aware runtime dispatch macros (Fat Binary) ───────────────
# Detect host processor and enable compile-time macros that allow the code to
# compile AVX2/NEON variants while deferring runtime selection to the
# dispatcher (no global -mavx2 is applied; avoids SIGILL on older CPUs).
include(CheckCXXCompilerFlag)
string(TOLOWER "${CMAKE_HOST_SYSTEM_PROCESSOR}" TRACEON_HOST_PROC)
message(STATUS "TracEon: host processor = ${TRACEON_HOST_PROC}")
if(TRACEON_HOST_PROC MATCHES "x86_64|amd64")
    # Check if the compiler understands -mavx2; if so, define TRACEON_HAS_AVX2
    check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_MAVX2)
    if(COMPILER_SUPPORTS_MAVX2)
        target_compile_definitions(traceon_core PUBLIC TRACEON_HAS_AVX2)
        message(STATUS "TracEon: Compiler supports AVX2; defining TRACEON_HAS_AVX2")
    endif()
elseif(TRACEON_HOST_PROC MATCHES "aarch64|arm64")
    target_compile_definitions(traceon_core PUBLIC TRACEON_HAS_NEON)
    message(STATUS "TracEon: ARM64 target detected; defining TRACEON_HAS_NEON")
endif()
# NEON on AArch64 requires no extra compiler flags; runtime dispatch will call NEON path.

# --- Executable Targets ---
add_executable(trac-eon-test main.cpp)
target_link_libraries(trac-eon-test PRIVATE traceon_core)

# --- Benchmark Driver ---
add_executable(traceon_driver traceon_driver.cpp)
target_link_libraries(traceon_driver PRIVATE traceon_core)

# --- Example Usage (NEW) ---
add_executable(example_usage examples/simple_usage.cpp)
target_link_libraries(example_usage PRIVATE traceon_core)

# --- Test Setup using Catch2 ---
include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.6.0 
)
FetchContent_MakeAvailable(Catch2)

enable_testing()

add_executable(unit_tests ${TEST_SOURCES})
target_link_libraries(unit_tests PRIVATE traceon_core Catch2::Catch2WithMain)

include(CTest)
add_test(NAME TracEonTests COMMAND unit_tests)