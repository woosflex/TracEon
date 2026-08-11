# Sticking to a modern, stable CMake version.
cmake_minimum_required(VERSION 3.20)

# This is us! Our project's name is TracEon.
project(TracEon)

# We're using C++20 for all the nice modern features.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- ZLIB Support (zlib-ng via FetchContent, zlib-compat mode) ---
include(FetchContent)
FetchContent_Declare(
        zlib-ng
        GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git
        GIT_TAG        2.2.2
)
set(ZLIB_COMPAT ON)
set(WITH_NATIVE_INSTRUCTIONS ON)
FetchContent_MakeAvailable(zlib-ng)

# --- ankerl::unordered_dense (header-only, MIT license) ---
FetchContent_Declare(
        ankerl_unordered_dense
        GIT_REPOSITORY https://github.com/martinus/unordered_dense.git
        GIT_TAG        v4.4.0
)
FetchContent_MakeAvailable(ankerl_unordered_dense)

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
        tests/MapDefsTests.cpp
        # SmartStrategyTests.cpp split into per-domain TUs (128 cases / 4195
        # assertions preserved across the combined unit_tests target).
        tests/test_parser_fastq.cpp
        tests/test_parser_fasta.cpp
        tests/test_gzip.cpp
        tests/test_binary_cache.cpp
        tests/test_lifecycle.cpp
        tests/test_api_misc.cpp
)

# --- Library Target ---
add_library(traceon_core STATIC ${CORE_SOURCES})

# Include directories (UPDATED PATHS)
# traceon_core users need 'include/'
# traceon_core implementation needs 'third_party/lz4/lib' and 'third_party/'
target_include_directories(traceon_core PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lz4/lib
    ${CMAKE_CURRENT_SOURCE_DIR} # For third_party/robin_hood.h inclusion
    ${ankerl_unordered_dense_SOURCE_DIR}/include # ankerl::unordered_dense
)

# GCC 16.1.1+ removed implicit <cstdint> includes; robin_hood.h needs UINT64_C
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 16)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include cstdint")
endif()

# Link traceon_core to both zlib and lz4
target_link_libraries(traceon_core PUBLIC zlibstatic lz4_static)

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

# ── Per-domain split test executables ────────────────────────────────────────
# SmartStrategyTests.cpp was split into per-domain translation units so each
# test executable is small enough to start under AddressSanitizer on hosts
# where the ~42 MB combined unit_tests binary fails at ASan init
# (InternalMmapVector ENOMEM — binary-size/ASLR interaction). Each domain exe
# links Catch2::Catch2WithMain (the single TU providing CATCH_CONFIG_MAIN).
#
# The [oom] test (setrlimit RLIMIT_AS) aborts under ASan shadow-mapping; run
# these under ASan with the [oom] tag filtered:  <exe> "~[oom]"
set(SPLIT_TEST_TARGETS
        test_parser_fastq
        test_parser_fasta
        test_gzip
        test_binary_cache
        test_lifecycle
        test_api_misc
)
foreach(t ${SPLIT_TEST_TARGETS})
    add_executable(${t} tests/${t}.cpp)
    target_link_libraries(${t} PRIVATE traceon_core Catch2::Catch2WithMain)
    add_test(NAME ${t} COMMAND ${t})
endforeach()

# ── K-mer experimental targets (OUT of the core umbrella) ──────────────────
# traceon_kmer is deliberately NOT part of CORE_SOURCES, does not appear in
# the public umbrella header (TracEon.h), and is not linked into
# traceon_core. It is an opt-in, separately built static library covering
# the C shim for the minimap2-style minimizer table
# (src/kmerindex_c_api.cpp) plus the header-only k-mer encoders/indexes
# (KmerEncoding.h, KmerReferenceIndex.h).
add_library(traceon_kmer STATIC src/kmerindex_c_api.cpp)
# Consumers need include/ for kmerindex_c_api.h + the header-only k-mer
# headers; ankerl is an implementation detail of the shim only.
target_include_directories(traceon_kmer PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR} # For third_party/robin_hood.h inclusion
)
target_include_directories(traceon_kmer PRIVATE
    ${ankerl_unordered_dense_SOURCE_DIR}/include # ankerl::unordered_dense
)

# ── K-mer unit tests (Catch2) ───────────────────────────────────────────────
add_executable(kmer_tests tests/kmer_tests.cpp)
target_link_libraries(kmer_tests PRIVATE traceon_kmer Catch2::Catch2WithMain)
add_test(NAME KmerTests COMMAND kmer_tests)

# ── C compilation/linkage proof (the .c file is compiled as C99, proving
# the kmerindex_c_api.h boundary is real C, not C++ masquerading) ──────────
add_executable(kmer_c_compile_test tests/kmer_c_api_compile_test.c)
set_source_files_properties(tests/kmer_c_api_compile_test.c PROPERTIES LANGUAGE C)
set_property(TARGET kmer_c_compile_test PROPERTY C_STANDARD 99)
set_property(TARGET kmer_c_compile_test PROPERTY C_STANDARD_REQUIRED ON)
target_link_libraries(kmer_c_compile_test PRIVATE traceon_kmer)
# The executable has only a C source, but links a C++ static library, so
# drive the link with the C++ driver to pull in libstdc++.
set_target_properties(kmer_c_compile_test PROPERTIES LINKER_LANGUAGE CXX)
add_test(NAME KmerCCompileTest COMMAND kmer_c_compile_test)