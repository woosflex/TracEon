// fuzz_gzip_loader.cpp — libFuzzer target for the GZIP load path
// (SmartStrategy::loadGzipFile / Cache::loadFile-with-gzip-detection).
//
// Contract under test: arbitrary bytes must either decompress+parse or be
// rejected with std::runtime_error (truncated streams, trailing garbage,
// corrupt members, OOM guards) — never a crash, an OOB read, an infinite
// loop, or a leak.
//
// The temp file is written with a ".gz" suffix so loadFile() takes its
// extension-based routing (in addition to the magic-byte check) when the
// input happens to start with the gzip magic.

#include "SmartStrategy.h"
#include "Cache.h"
#include "fuzz_common.h"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace {

void exercise_reads(TracEon::SmartStrategy& s) {
    volatile size_t sink = 0;
    const auto keys = s.getAllKeys();
    for (const auto& k : keys) {
        const std::string_view v = s.getView(k);
        for (char c : v) sink += static_cast<unsigned char>(c);
    }
    (void)sink;
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    // Harness-level amplification guard (return 0, NOT a finding): never feed
    // more than 8 MiB of compressed input. A gzip member can amplify up to
    // ISIZE (~4 GiB) from a tiny input; the loader pre-sizes its arena from
    // the input (compressed_size x 3, OOM-guarded) and has its own bad_alloc
    // catch, but under ASan a failed huge malloc aborts the process BEFORE
    // bad_alloc propagates (same limitation as the [oom] test). The CI budget
    // is max_len=262144, so this 8 MiB cap is already effectively true — it is
    // explicit here so any local run without max_len stays bounded. Nothing is
    // decompressed in the harness; the loader remains the authoritative path.
    if (size > 8ull * 1024 * 1024) return 0;
    const std::string path = traceon_fuzz::write_input_to_tmp(data, size, ".gz");
    // Alternate entry points: the explicit loadGzipFile() and the
    // auto-detecting loadFile() (routes by extension, then by magic bytes).
    const bool via_load_file = (size > 0) && ((data[0] & 1) != 0);
    try {
        if (via_load_file) {
            TracEon::Cache cache;
            cache.loadFile(path);
        } else {
            TracEon::SmartStrategy s;
            s.loadGzipFile(path);
            exercise_reads(s);
        }
    } catch (const std::exception&) {
        // Expected rejection for malformed/truncated gzip input.
    }
    return 0;
}
