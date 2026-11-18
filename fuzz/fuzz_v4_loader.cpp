// fuzz_v4_loader.cpp — libFuzzer target for `.traceon` v4 migration and
// v5 development-format loading (version-aware headers + LZ4 Frame + CRC32C).
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

// ── Harness-level amplification guard (return 0, never a finding) ──────────
// A fuzz input is capped at max_len=262144 bytes in CI, so a v4/v5 header
// claiming > 64 MiB of logical payload (or frame) can never be legitimate:
// the loader cannot synthesize > 64 MiB of uncompressed data from <= 256 KiB
// of input. Skipping such inputs prevents pathological multi-GB allocations
// entirely. The loader's own OOM guard (bad_alloc catch around the arena
// resize) is the designed defense, but under ASan a failed huge malloc aborts
// the process BEFORE bad_alloc can propagate (same limitation as the [oom]
// unit test) — hence this harness-side cap. The CI workflow additionally
// exports ASAN_OPTIONS=allocator_may_return_null=1 so any residual huge
// allocation fails gracefully into the bad_alloc path instead of aborting.
//
// v4: logical/frame lengths at offsets 6/14.
// v5: format byte at 6, then logical/frame lengths at offsets 7/15.
static constexpr uint64_t kMaxDeclaredLen = 64ull * 1024 * 1024; // 64 MiB

// Explicit little-endian decode (host may be big-endian; memcpy-free is
// fine here since we only ever read within the fuzz buffer).
inline uint64_t read_le64(const uint8_t* p) noexcept {
    uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v |= static_cast<uint64_t>(p[i]) << (8 * i);
    return v;
}

// True => skip the input. Only the header prefix is parsed (magic + the two
// u64 length fields); nothing else is touched. Bounds are exact: logical_len
// offsets are selected only after a complete TRO version prefix; shorter or
// unknown-version inputs are left to the loader's normal rejection paths.
inline bool declares_implausible_lengths(const uint8_t* data, size_t size) {
    if (size < 4 || std::memcmp(data, "TRO", 3) != 0) return false;
    const uint8_t version = data[3];
    const size_t logical_offset = version == 4 ? 6 : version == 5 ? 7 : 0;
    const size_t frame_offset = version == 4 ? 14 : version == 5 ? 15 : 0;
    if (logical_offset == 0 || size < logical_offset + 8) return false;
    if (read_le64(data + logical_offset) > kMaxDeclaredLen) return true;
    if (size >= frame_offset + 8 && read_le64(data + frame_offset) > kMaxDeclaredLen) return true;
    return false;
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    // Amplification guard: skip (return 0, NOT a finding) inputs whose v4/v5
    // header claims a logical/frame length that a <=262144-byte fuzz input
    // can never legitimately produce. Keeps the fuzzer focused on real paths
    // and prevents pathological allocations even outside ASan.
    if (declares_implausible_lengths(data, size)) return 0;
    const std::string path = traceon_fuzz::write_input_to_tmp(data, size);
    // Alternate between the two public entry points so both get coverage:
    // SmartStrategy::loadBinary (core) and Cache::restore (public API).
    for (bool via_cache : {false, true}) {
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
                    TracEon::SmartStrategy restored;
                    restored.loadBinary(out);
                    if (restored.getFileCacheSize() != s.getFileCacheSize() ||
                        restored.getDetectedFormat() != s.getDetectedFormat()) std::abort();
                    for (const auto& key : s.getAllKeys()) {
                        if (!restored.hasSequence(key) || restored.getView(key) != s.getView(key) ||
                            restored.getQuality(key) != s.getQuality(key)) std::abort();
                    }
                }
            }
        } catch (const std::exception&) {
            // Expected: malformed/implausible/truncated/corrupt caches must be
            // rejected with an exception, never a crash.
        }
    }
    return 0;
}
