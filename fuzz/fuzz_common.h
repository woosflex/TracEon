#ifndef TRACEON_FUZZ_COMMON_H
#define TRACEON_FUZZ_COMMON_H

// fuzz_common.h — shared helpers for the TracEon libFuzzer targets
// (OSS-Fuzz style, clang-only when built with TRACEON_BUILD_FUZZERS=ON).
//
// 1. Portable file-fed main. libFuzzer (clang -fsanitize=fuzzer) links its
//    own driver (main); note that clang does NOT define
//    FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION on its own — the fuzz CMake
//    build passes it explicitly (OSS-Fuzz convention) so this file-fed main
//    is suppressed. GCC / plain-clang builds — where that macro is NOT
//    defined — get the guarded main below instead: it reads every argv file
//    and feeds its bytes to LLVMFuzzerTestOneInput. The exit code is 0 as
//    long as every file was read and fed, regardless of what the target
//    does with the bytes (each target must never let an exception escape).
//    This lets the exact same .cpp compile and run on GCC-only hosts like
//    this one.
//
// 2. Temp-file helper. The loaders under test take a FILE PATH
//    (loadBinary/loadGzipFile/loadFile/mmap_open), but libFuzzer feeds
//    in-memory buffers, so the harness writes each input to a per-process
//    scratch dir under /tmp before invoking the API. A per-pid dir keeps
//    -fork=N / parallel CI jobs from colliding on file names.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#ifndef _WIN32
#include <unistd.h>
#include <sys/stat.h>
#else
#include <process.h>
#include <direct.h>
#endif

namespace traceon_fuzz {

// Per-process scratch dir: "/tmp/traceon_fuzz_<pid>". Created lazily on the
// first write. Never cleaned up (a fuzz process is short-lived; leftover
// files are harmless and /tmp is ephemeral).
inline std::string& tmp_dir() {
    static std::string dir;
    return dir;
}

inline void ensure_tmp_dir() {
    std::string& d = tmp_dir();
    if (!d.empty()) return;
    char buf[64];
#ifdef _WIN32
    std::snprintf(buf, sizeof(buf), "%s\\traceon_fuzz_%ld", std::getenv("TEMP") ? std::getenv("TEMP") : ".", static_cast<long>(_getpid()));
#else
    std::snprintf(buf, sizeof(buf), "/tmp/traceon_fuzz_%ld", static_cast<long>(getpid()));
#endif
    d = buf;
#ifdef _WIN32
    _mkdir(d.c_str());
#else
    ::mkdir(d.c_str(), 0700);
#endif
}

// Write `data` into a fresh file in the scratch dir. Returns the path.
// `suffix` (e.g. ".gz") lets targets control the temp name for
// extension-based format detection in loadFile().
// Aborts only on harness-level failure (unwritable temp dir) — never on
// fuzzer input.
inline std::string write_input_to_tmp(const uint8_t* data, size_t size,
                                      const char* suffix = "") {
    ensure_tmp_dir();
    static unsigned long counter = 0;
    std::string path = tmp_dir() + "/in" + std::to_string(counter++) + suffix;
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) std::abort(); // harness bug: cannot create scratch file
    if (size > 0 && std::fwrite(data, 1, size, f) != size) {
        std::fclose(f);
        std::abort(); // harness bug: short write
    }
    std::fclose(f);
    return path;
}

} // namespace traceon_fuzz

// Declared here so the non-libFuzzer main below can call it; each fuzz_*.cpp
// defines it with `extern "C"` linkage — REQUIRED because libFuzzer's driver
// declares LLVMFuzzerTestOneInput as extern "C". Without it the definition is
// C++-mangled and does not resolve the driver's symbol, so linking fails with
// 'undefined reference to LLVMFuzzerTestOneInput'. GCC's guarded main below
// calls the same extern "C" entity directly.
extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size);

#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
// File-fed main for GCC / non-libFuzzer builds. Reads each argv file and
// feeds it to the target. Exit 0 on all inputs; nonzero only for a
// file-open/read failure (a real harness/environment problem, not a
// finding — targets catch everything themselves).
int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <input-file> [input-file ...]\n", argv[0]);
        return 1;
    }
    int rc = 0;
    for (int i = 1; i < argc; ++i) {
        FILE* f = std::fopen(argv[i], "rb");
        if (!f) {
            std::fprintf(stderr, "%s: cannot open\n", argv[i]);
            rc = 2;
            continue;
        }
        if (std::fseek(f, 0, SEEK_END) != 0) { std::fclose(f); rc = 3; continue; }
        const long sz = std::ftell(f);
        if (sz < 0) { std::fclose(f); rc = 3; continue; }
        std::rewind(f);
        std::vector<uint8_t> buf(static_cast<size_t>(sz));
        if (sz > 0 && std::fread(buf.data(), 1, buf.size(), f) != buf.size()) {
            std::fprintf(stderr, "%s: read error\n", argv[i]);
            std::fclose(f);
            rc = 3;
            continue;
        }
        std::fclose(f);
        LLVMFuzzerTestOneInput(buf.data(), buf.size());
    }
    return rc;
}
#endif // FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION

#endif // TRACEON_FUZZ_COMMON_H
