// fuzz_v4_loader.cpp — libFuzzer target for the `.traceon` v4 binary
// loader (the saveBinary() format: "TRO\x04" + 22-byte header + streamed
// LZ4 Frame + whole-payload CRC32C).
//
// Contract under test: arbitrary bytes fed to SmartStrategy::loadBinary()
// / Cache::restore() must be rejected cleanly (std::runtime_error) or load
// successfully — never a crash, an out-of-bounds read/write, an infinite
// loop, or a leak. A crash here is a finding.
//
// The read path is exercised on any successfully-loaded cache (getAllKeys /
// getView over the arena-backed string_views) and, for small caches, a
// saveBinary() round-trip — so serialization bugs surface too.

#include "SmartStrategy.h"
#include "Cache.h"
#include "fuzz_common.h"

#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

namespace {

// Bounded read-path exercise on a loaded cache. Views are consumed byte-by-
// byte so ASan/UBSan catch any out-of-bounds arena string_view.
void exercise_reads(TracEon::SmartStrategy& s) {
    volatile size_t sink = 0;
    const auto keys = s.getAllKeys();
    for (const auto& k : keys) {
        const std::string_view v = s.getView(k);
        for (char c : v) sink += static_cast<unsigned char>(c);
        sink += s.getSequence(k).size();
        sink += s.getQuality(k).size();
    }
    (void)sink;
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    const std::string path = traceon_fuzz::write_input_to_tmp(data, size);
    // Alternate between the two public entry points so both get coverage:
    // SmartStrategy::loadBinary (core) and Cache::restore (public API).
    const bool via_cache = (size > 0) && ((data[0] & 1) != 0);
    try {
        if (via_cache) {
            TracEon::Cache cache;
            cache.restore(path);
            if (cache.size() > 0) {
                // Public read path with keys derived from the input bytes.
                volatile size_t sink = 0;
                for (size_t off = 0; off + 8 < size && off < 256; off += 9) {
                    const std::string key(reinterpret_cast<const char*>(data + off), 8);
                    sink += cache.getView(key).size();
                    sink += cache.hasSequence(key) ? 1 : 0;
                }
                (void)sink;
            }
        } else {
            TracEon::SmartStrategy s;
            s.loadBinary(path);
            exercise_reads(s);
            // Round-trip: a small successfully-loaded cache must re-serialize.
            if (s.getFileCacheSize() < 10'000) {
                const std::string out = path + ".save";
                s.saveBinary(out);
            }
        }
    } catch (const std::exception&) {
        // Expected: malformed/implausible/truncated/corrupt caches must be
        // rejected with an exception, never a crash.
    }
    return 0;
}
