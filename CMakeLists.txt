# Sticking to a modern, stable CMake version.
cmake_minimum_required(VERSION 3.20)

# This is us! Our project's name is TracEon.
project(TracEon VERSION 3.0.0 LANGUAGES C CXX)
include(GNUInstallDirs)
include(CTest)

# We're using C++20 for all the nice modern features.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Threads are used by the remote server (traceon_remote) and the benchmarks.
find_package(Threads REQUIRED)

# --- ZLIB Support (zlib-ng via FetchContent, zlib-compat mode) ---
include(FetchContent)
FetchContent_Declare(
        zlib-ng
        GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git
        # zlib-ng 2.2.2, pinned to an immutable commit for reproducibility.
        GIT_TAG        c939498c7f2a1d645833a98251365e2c814fd736
)
set(ZLIB_COMPAT ON CACHE BOOL "Build zlib-ng with the zlib-compatible API" FORCE)
option(WITH_NATIVE_INSTRUCTIONS "Tune zlib-ng for this build host (non-portable)" OFF)
set(ZLIB_ENABLE_TESTS OFF CACHE BOOL "Disable bundled zlib-ng tests" FORCE)
set(ZLIBNG_ENABLE_TESTS OFF CACHE BOOL "Disable bundled zlib-ng tests" FORCE)
set(WITH_GTEST OFF CACHE BOOL "Disable zlib-ng GoogleTest targets" FORCE)
set(SKIP_INSTALL_ALL ON) # Install the pinned static dependency privately with TracEon.
FetchContent_MakeAvailable(zlib-ng)
unset(SKIP_INSTALL_ALL)

# --- ankerl::unordered_dense (header-only, MIT license) ---
FetchContent_Declare(
        ankerl_unordered_dense
        GIT_REPOSITORY https://github.com/martinus/unordered_dense.git
        # Immutable commit used by the audited candidate (API-compatible with v4.4.0).
        GIT_TAG        231e48c9426bd21c273669e5fdcd042c146975cf
)
FetchContent_MakeAvailable(ankerl_unordered_dense)

# --- LZ4 Support ---
set(LZ4_BUILD_CLI OFF CACHE BOOL "Disable bundled LZ4 command-line tools" FORCE)
add_subdirectory(third_party/lz4/build/cmake)

# --- Source File Definitions (UPDATED PATHS) ---
set(CORE_SOURCES
        src/Cache.cpp
        src/SmartStrategy.cpp
)

# ── Remote read access (v2.3.0 "Harpe") ──────────────────────────────────────
# POSIX-only by design: the server/client use sockets and are never linked
# into traceon_core (core stays free of socket code — the remote lib links
# the core, not vice versa).
option(TRACEON_BUILD_SERVER "Build TracEon remote read access (POSIX sockets)" ON)
if(TRACEON_BUILD_SERVER AND NOT WIN32)
    find_package(OpenSSL REQUIRED COMPONENTS Crypto)
    set(REMOTE_SOURCES
            src/TraceonServer.cpp
            src/TraceonClient.cpp
    )
    add_library(traceon_remote STATIC ${REMOTE_SOURCES})
    target_include_directories(traceon_remote PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/TracEon>)
    target_compile_features(traceon_remote PUBLIC cxx_std_20)
    set_target_properties(traceon_remote PROPERTIES EXPORT_NAME remote POSITION_INDEPENDENT_CODE ON)
    add_library(TracEon::remote ALIAS traceon_remote)
    target_link_libraries(traceon_remote PUBLIC traceon_core Threads::Threads)
    target_link_libraries(traceon_remote PRIVATE OpenSSL::Crypto)
    # The socket helpers in TraceonProto.h need the POSIX feature-test macros
    # at compile time (GCC 16's -include cstdint pulls libc headers in before
    # any in-file #define would take effect).
    target_compile_definitions(traceon_remote PRIVATE _GNU_SOURCE)
else()
    message(STATUS "TracEon: remote server/client disabled (TRACEON_BUILD_SERVER=OFF or WIN32)")
endif()

set(TEST_SOURCES
        tests/CacheTests.cpp
        tests/FastqTests.cpp
        tests/MapDefsTests.cpp
)

# --- Library Target ---
add_library(traceon_core STATIC ${CORE_SOURCES})

# Include directories (UPDATED PATHS)
# traceon_core users need 'include/'
# traceon_core implementation needs 'third_party/lz4/lib' and 'third_party/'
target_include_directories(traceon_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${ankerl_unordered_dense_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/TracEon>)
target_compile_features(traceon_core PUBLIC cxx_std_20)
set_target_properties(traceon_core PROPERTIES EXPORT_NAME core POSITION_INDEPENDENT_CODE ON)
add_library(TracEon::core ALIAS traceon_core)

# GCC 16.1.1+ removed implicit <cstdint> includes; robin_hood.h needs UINT64_C
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 16)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include cstdint")
endif()

# Link traceon_core to both zlib and lz4
target_link_libraries(traceon_core PUBLIC zlibstatic lz4_static Threads::Threads)

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

# ── Remote access tools (v2.3.0) ────────────────────────────────────────────
if(TRACEON_BUILD_SERVER AND NOT WIN32)
    # Benchmark: --mode=local (in-process getView baseline) | --mode=remote
    # (TCP client-server) | --mode=serve (embedded server; used by the Docker
    # testbed and manual verification).
    add_executable(remote_bench benchmarks/remote_bench.cpp)
    target_link_libraries(remote_bench PRIVATE traceon_remote Threads::Threads)
    target_compile_definitions(remote_bench PRIVATE _GNU_SOURCE)

    add_executable(traceon tools/traceon.cpp)
    target_link_libraries(traceon PRIVATE traceon_remote)
endif()

# Cache construction also works in core-only builds.
add_executable(traceon_make_cache tools/traceon_make_cache.cpp)
target_link_libraries(traceon_make_cache PRIVATE traceon_core)

# --- Benchmark Driver ---
add_executable(traceon_driver traceon_driver.cpp)
target_link_libraries(traceon_driver PRIVATE traceon_core)

# Opt-in FASTA subset adapter for the mapping feasibility experiment.
add_executable(traceon_subset EXCLUDE_FROM_ALL benchmarks/traceon_subset.cpp)
target_link_libraries(traceon_subset PRIVATE traceon_core)
if(TARGET traceon_remote)
    target_link_libraries(traceon_subset PRIVATE traceon_remote)
    target_compile_definitions(traceon_subset PRIVATE TRACEON_SUBSET_REMOTE)
endif()

# ── Read-scaling benchmark (v2.3.0 "Harpe", slice 2) ────────────────────────
# Standalone benchmark target: lock-free getView() read throughput vs thread
# count (1/2/4/8/14) against one in-memory immutable cache built via the
# addEntry() path. Core-only (no sockets), so it is NOT guarded by
# TRACEON_BUILD_SERVER. Deliberately NOT registered with CTest and NOT part
# of unit_tests (same policy as remote_bench) so the default suite stays
# fast; it is a benchmark target run manually/CI-on-demand:
#
#   ./build/read_scaling [--records=N --keys=N --count=N --warmup=N
#                         --threads=1,2,4,8,14 --reps=R --seed=S]
add_executable(read_scaling benchmarks/read_scaling.cpp)
target_link_libraries(read_scaling PRIVATE traceon_core Threads::Threads)

# --- Example Usage (NEW) ---
add_executable(example_usage examples/simple_usage.cpp)
target_link_libraries(example_usage PRIVATE traceon_core)

# ── K-mer experimental targets (OUT of the core umbrella) ──────────────────
# traceon_kmer is deliberately NOT part of CORE_SOURCES, does not appear in
# the public umbrella header (TracEon.h), and is not linked into
# traceon_core. It is an opt-in, separately built static library covering
# the C shim for the minimap2-style minimizer table
# (src/kmerindex_c_api.cpp) plus the header-only k-mer encoder
# (KmerEncoding.h).
add_library(traceon_kmer STATIC src/kmerindex_c_api.cpp)
# Consumers need include/ for kmerindex_c_api.h + the header-only k-mer
# headers; ankerl is an implementation detail of the shim only.
target_include_directories(traceon_kmer PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/TracEon>
)
target_include_directories(traceon_kmer PRIVATE
    ${ankerl_unordered_dense_SOURCE_DIR}/include # ankerl::unordered_dense
)

target_compile_features(traceon_kmer PUBLIC cxx_std_20)
set_target_properties(traceon_kmer PROPERTIES EXPORT_NAME kmer POSITION_INDEPENDENT_CODE ON)
add_library(TracEon::kmer ALIAS traceon_kmer)

if(BUILD_TESTING)
# --- Test Setup using Catch2 ---
include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        # Immutable commit used by the audited candidate (Catch2 3.6 series).
        GIT_TAG        4e8d92bf02f7d1c8006a0e7a5ecabd8e62d98502
)
FetchContent_MakeAvailable(Catch2)

enable_testing()

add_executable(unit_tests ${TEST_SOURCES})
target_link_libraries(unit_tests PRIVATE traceon_core Catch2::Catch2WithMain)

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
target_compile_definitions(test_binary_cache PRIVATE
    TRACEON_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

# Remote domain tests need the remote lib (server/client classes).
if(TRACEON_BUILD_SERVER AND NOT WIN32)
    add_executable(test_remote tests/test_remote.cpp)
    target_link_libraries(test_remote PRIVATE traceon_remote traceon_core Catch2::Catch2WithMain)
    target_compile_definitions(test_remote PRIVATE _GNU_SOURCE)
    add_test(NAME test_remote COMMAND test_remote)
endif()

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

# Focused regressions from the correctness audit.
add_executable(test_save_allocations tests/test_save_allocations.cpp)
target_link_libraries(test_save_allocations PRIVATE traceon_core)
add_test(NAME test_save_allocations COMMAND test_save_allocations)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    foreach(check indexed_mapping_study summarize_paired)
        add_test(NAME ${check} COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_${check}.py)
    endforeach()
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set_target_properties(traceon_subset PROPERTIES EXCLUDE_FROM_ALL FALSE)
        add_test(NAME mapping_baseline COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mapping_baseline.py)
        add_test(NAME subset_mapping_pilot COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_subset_mapping_pilot.py
            $<TARGET_FILE:traceon_make_cache> $<TARGET_FILE:traceon_subset>)
        if(TARGET remote_bench)
            add_test(NAME shared_mapping_pilot COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_shared_mapping_pilot.py
                $<TARGET_FILE:traceon_make_cache> $<TARGET_FILE:traceon_subset>
                $<TARGET_FILE:remote_bench>)
        endif()
    endif()
    if(TARGET traceon)
        add_test(NAME service_cli COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_service_cli.py $<TARGET_FILE:traceon>)
        add_test(NAME service_cli_mapped COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_service_cli.py $<TARGET_FILE:traceon> --mapped)
        add_test(NAME remote_profile COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_remote_profile.py
            $<TARGET_FILE:traceon> $<TARGET_FILE:remote_bench>)
    endif()
    set(AUDIT_CHECK_ARGS --build-dir ${CMAKE_CURRENT_BINARY_DIR})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND NOT CMAKE_CXX_FLAGS MATCHES "sanitize")
        add_library(atomic_file_faults SHARED tests/atomic_file_faults.c)
        target_link_libraries(atomic_file_faults PRIVATE ${CMAKE_DL_LIBS})
        if(TARGET traceon_make_cache)
            add_test(NAME atomic_save COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_atomic_save.py
                $<TARGET_FILE:traceon_make_cache> $<TARGET_FILE:traceon_driver>
                $<TARGET_FILE:atomic_file_faults>)
            add_test(NAME atomic_save_mapped COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_atomic_save.py
                $<TARGET_FILE:traceon_make_cache> $<TARGET_FILE:traceon_driver>
                $<TARGET_FILE:atomic_file_faults> --mapped)
        endif()
        add_library(thread_failure_inject SHARED tests/thread_failure_inject.c)
        target_link_libraries(thread_failure_inject PRIVATE ${CMAKE_DL_LIBS})
        list(APPEND AUDIT_CHECK_ARGS --thread-shim $<TARGET_FILE:thread_failure_inject>)
    endif()
    add_test(NAME audit_checks COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/audit_checks.py ${AUDIT_CHECK_ARGS})
endif()

endif() # BUILD_TESTING

include(cmake/Install.cmake)

# ── Fuzzing targets (libFuzzer / OSS-Fuzz style) ─────────────────────────────
# Off by default: fuzz/*.cpp are NOT compiled at all on normal builds (they
# would compile fine with GCC too — each carries a portable file-fed main
# guarded by #ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION — but keeping
# them out of the default build avoids dead targets and sanitizer flag noise
# in the regular build/test flow).
#
# When ON, clang is REQUIRED (-fsanitize=fuzzer is a clang-only driver) and
# each fuzz target is built with -fsanitize=fuzzer,address,undefined. The
# fuzzers are deliberately NOT registered with CTest: CI runs them with an
# explicit time budget (see .github/workflows/fuzz.yml) instead.
#
# Whole-program sanitizer coverage (canonical OSS-Fuzz pattern): the core/
# kmer libraries are instrumented with -fsanitize=fuzzer-no-link,address,
# undefined. fuzzer-no-link adds libFuzzer's coverage callbacks and the
# sanitizer instrumentation but links NO runtime (neither libFuzzer nor
# ASan/UBSan) into the static libs, so nothing propagates to consumers.
# Each fuzz target then adds its own single -fsanitize=fuzzer,address,
# undefined on compile AND link — the only runtime spec on the final link
# line, so no duplicate sanitizer runtimes. (clang 18/19 rejects two ASan
# runtime specs with 'incompatible ASan runtimes'; clang 22 tolerated the
# duplicate, which is why the old CMAKE_EXE_LINKER_FLAGS approach worked
# locally on Fedora but failed on ubuntu-latest.) The recommended configure
# is exactly what .github/workflows/fuzz.yml does:
#
#   cmake -B build-fuzz -DTRACEON_BUILD_FUZZERS=ON \
#         -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=RelWithDebInfo
option(TRACEON_BUILD_FUZZERS "Build libFuzzer fuzz targets (requires clang)" OFF)

if(TRACEON_BUILD_FUZZERS)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(FATAL_ERROR
            "TRACEON_BUILD_FUZZERS=ON requires clang (CMAKE_CXX_COMPILER_ID is "
            "'${CMAKE_CXX_COMPILER_ID}'). Configure with -DCMAKE_CXX_COMPILER=clang++")
    endif()

    # Instrument the libraries so sanitizers see into SmartStrategy.cpp,
    # kmerindex_c_api.cpp, and the header-only parsers/encoders. PRIVATE +
    # fuzzer-no-link: coverage/sanitizer instrumentation only — no runtime
    # is linked into the static libs and nothing propagates to consumers
    # (the fuzz targets supply the single -fsanitize=fuzzer,address,undefined
    # spec on their own compile+link lines).
    target_compile_options(traceon_core PRIVATE -fsanitize=fuzzer-no-link,address,undefined -fno-omit-frame-pointer)
    target_compile_options(traceon_kmer PRIVATE -fsanitize=fuzzer-no-link,address,undefined -fno-omit-frame-pointer)

    set(TRACEON_FUZZ_TARGETS
        fuzz_v4_loader
        fuzz_gzip_loader
        fuzz_fastq
        fuzz_fasta
        fuzz_kmer_encode
        fuzz_kmer_api
    )
    foreach(fuzz_target ${TRACEON_FUZZ_TARGETS})
        add_executable(${fuzz_target} fuzz/${fuzz_target}.cpp)
        target_include_directories(${fuzz_target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fuzz)
        target_link_libraries(${fuzz_target} PRIVATE traceon_core traceon_kmer)
        # libFuzzer does NOT define FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
        # (verified on clang 22: -fsanitize=fuzzer adds no such macro); pass it
        # explicitly so fuzz_common.h suppresses its file-fed main() and lets
        # libFuzzer's driver main win (OSS-Fuzz convention). GCC file-fed
        # builds (outside this configure) still get the guarded main.
        target_compile_definitions(${fuzz_target} PRIVATE FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION)
        target_compile_options(${fuzz_target} PRIVATE
            -fsanitize=fuzzer,address,undefined
            -fno-omit-frame-pointer)
        target_link_options(${fuzz_target} PRIVATE -fsanitize=fuzzer,address,undefined)
        # NOT registered with CTest (see comment above).
    endforeach()
endif()
