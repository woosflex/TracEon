// K-mer experimental target tests (kmer_tests).
//
// Covers: encode_kmer guards (vuln-0001), C API happy path + error model
// (vuln-0003: no exception can cross the C boundary), null-handle
// validation, freeze/pointer-stability enforcement, caller-owned iterator
// correctness + reentrancy, and KmerReferenceIndex TRKI malformed-file
// rejection + valid round-trip (vuln-0006).

#include <catch2/catch_test_macros.hpp>

#include "KmerEncoding.h"
#include "KmerReferenceIndex.h"
#include "kmerindex_c_api.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <filesystem>
#include <set>
#include <string>
#include <string_view>
#include <thread>
#include <atomic>
#include <vector>
#include <stdexcept>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace fs = std::filesystem;

// ── helpers ─────────────────────────────────────────────────────────────────

// RAII restore of RLIMIT_AS (used by the exception-boundary test).
#ifndef _WIN32
struct RlimitAsGuard {
    struct rlimit original{};
    bool applied = false;
    explicit RlimitAsGuard(rlim_t limit_bytes) {
        if (getrlimit(RLIMIT_AS, &original) != 0) return;
        struct rlimit constrained = original;
        constrained.rlim_cur = limit_bytes;
        if (setrlimit(RLIMIT_AS, &constrained) == 0) applied = true;
    }
    ~RlimitAsGuard() {
        if (applied) setrlimit(RLIMIT_AS, &original);
    }
};
#endif

// Writes a TRKI file with explicit little-endian bytes, exactly like the
// format spec: magic(4) + version/u32 + k/u32 + reserved/u32 + n_slots/u64
// + n_keys/u64 + n_positions/u64, then slots then positions.
static void write_trki(const std::string& path,
                       uint32_t version, uint32_t k,
                       uint64_t n_slots, uint64_t n_keys, uint64_t n_positions,
                       const std::vector<std::pair<uint64_t, uint64_t>>& slots = {},
                       const std::vector<uint64_t>& positions = {},
                       bool bad_magic = false) {
    std::ofstream f(path, std::ios::binary);
    REQUIRE(f.good());
    auto put32 = [&](uint32_t v) {
        char b[4];
        for (int i = 0; i < 4; ++i) b[i] = (char)((v >> (8 * i)) & 0xFFu);
        f.write(b, 4);
    };
    auto put64 = [&](uint64_t v) {
        char b[8];
        for (int i = 0; i < 8; ++i) b[i] = (char)((v >> (8 * i)) & 0xFFu);
        f.write(b, 8);
    };
    f.write(bad_magic ? "XXXX" : "TRKI", 4);
    put32(version);
    put32(k);
    put32(0); // reserved
    put64(n_slots);
    put64(n_keys);
    put64(n_positions);
    for (auto& [key, value] : slots) { put64(key); put64(value); }
    for (uint64_t p : positions) put64(p);
    f.close();
}

static void append_bytes(const std::string& path, const char* bytes, size_t n) {
    std::ofstream f(path, std::ios::binary | std::ios::app);
    f.write(bytes, (std::streamsize)n);
}

// ── encode_kmer (vuln-0001) ─────────────────────────────────────────────────

TEST_CASE("encode_kmer guards: short view / bad k never read out of bounds",
          "[kmer][encode][vuln-0001]") {
    using TracEon::encode_kmer;

    // Short view with oversized k -> false, out untouched (the vuln-0001 PoC).
    const std::string seq8(8, 'A');
    std::string_view short_view(seq8.data(), 8); // non-NUL-terminated view
    uint64_t out = 0xDEADBEEFDEADBEEFull;
    REQUIRE_FALSE(encode_kmer(short_view, 12, out));
    REQUIRE(out == 0xDEADBEEFDEADBEEFull); // out must be untouched on failure

    // Empty / truncated views.
    out = 7;
    REQUIRE_FALSE(encode_kmer("", 1, out));
    REQUIRE(out == 7);
    REQUIRE_FALSE(encode_kmer("ACG", 4, out)); // size 3 < k 4
    REQUIRE(out == 7);

    // k out of the valid [1, 32] range.
    REQUIRE_FALSE(encode_kmer("ACGTACGT", 0, out));
    REQUIRE_FALSE(encode_kmer("ACGTACGT", -3, out));
    REQUIRE_FALSE(encode_kmer(std::string(40, 'A'), 33, out)); // k > 32 packing limit

    // Valid encodes: 2-bit packing A=0 C=1 G=2 T=3, little endian of window.
    REQUIRE(encode_kmer("ACGT", 4, out));
    REQUIRE(out == 27); // 0b00011011
    REQUIRE(encode_kmer("acgt", 4, out)); // case-insensitive
    REQUIRE(out == 27);
    REQUIRE(encode_kmer("AAAA", 4, out));
    REQUIRE(out == 0);
    REQUIRE(encode_kmer("TTTT", 4, out));
    REQUIRE(out == 255); // 0b11111111

    // Only the first k bases are consumed.
    REQUIRE(encode_kmer("ACGTACGTACGT", 4, out));
    REQUIRE(out == 27);

    // Non-ACGT base inside the window -> false; outside the window -> fine.
    REQUIRE_FALSE(encode_kmer("ACNT", 4, out));
    REQUIRE(encode_kmer("ACGTN", 4, out)); // N at index 4 is outside [0,4)
    REQUIRE(out == 27);
}

TEST_CASE("encode_kmer packs 32-mers (max packing limit)", "[kmer][encode]") {
    using TracEon::encode_kmer;
    const std::string seq32(32, 'C'); // C=1 -> all bits 01
    uint64_t out = 0;
    REQUIRE(encode_kmer(seq32, 32, out));
    REQUIRE(out == 0x5555555555555555ull);
}

// ── C API happy path + minimap2 semantics ───────────────────────────────────

TEST_CASE("C API happy path: create/reserve/insert/get/iter/destroy",
          "[kmer][capi]") {
    kmerindex_t* h = kmerindex_create();
    REQUIRE(h != nullptr);
    REQUIRE(kmerindex_size(h) == 0);
    REQUIRE(kmerindex_reserve(h, 8) == 1);

    // minimap2 khash semantics: bit0 of the key is the single-occurrence
    // flag; hash and equality both ignore it.
    REQUIRE(kmerindex_insert(h, (2ull << 1) | 1u, 100) == 1); // final_key=5, base 2
    REQUIRE(kmerindex_insert(h, (5ull << 1) | 0u, 200) == 1); // final_key=10, base 5
    REQUIRE(kmerindex_insert(h, (9ull << 1) | 1u, 300) == 1); // final_key=19, base 9
    REQUIRE(kmerindex_size(h) == 3);

    // Collision: same base key as the first insert (base 2) -> 0, no error.
    REQUIRE(kmerindex_insert(h, (2ull << 1) | 0u, 999) == 0);
    REQUIRE(std::string(kmerindex_last_error(h)).empty());

    // get with query_key bit0=0 finds entries stored with bit0=1.
    uint64_t matched = 0;
    const uint64_t* v = kmerindex_get(h, 2ull << 1, &matched);
    REQUIRE(v != nullptr);
    REQUIRE(*v == 100);
    REQUIRE(matched == 5); // actual stored key carries the flag bit

    v = kmerindex_get(h, 5ull << 1, &matched);
    REQUIRE(v != nullptr);
    REQUIRE(*v == 200);
    REQUIRE(matched == 10);

    // Miss -> nullptr, no error.
    REQUIRE(kmerindex_get(h, 123ull << 1, &matched) == nullptr);
    REQUIRE(std::string(kmerindex_last_error(h)).empty());

    // Iteration sees every entry exactly once.
    kmerindex_iter_t it;
    REQUIRE(kmerindex_iter_begin(h, &it) == 1);
    std::set<std::pair<uint64_t, uint64_t>> seen;
    uint64_t k = 0, val = 0;
    while (kmerindex_iter_next(&it, &k, &val) == 1) seen.emplace(k, val);
    REQUIRE(seen == std::set<std::pair<uint64_t, uint64_t>>(
                        {{5, 100}, {10, 200}, {19, 300}}));
    // Past-the-end keeps returning 0.
    REQUIRE(kmerindex_iter_next(&it, &k, &val) == 0);

    kmerindex_destroy(h);
    kmerindex_destroy(nullptr); // NULL-safe
}

// ── Exception boundary (vuln-0003) ──────────────────────────────────────────

// ASan + setrlimit(RLIMIT_AS) are incompatible: lowering the address-space
// cap breaks ASan's shadow-memory mmaps (the core suite's OOM test has the
// same limitation). The absurd-size reserve below exercises the exception
// boundary without any rlimit, so it runs under every sanitizer too.
#if defined(__SANITIZE_ADDRESS__)
constexpr bool kUnderAsan = true;
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
constexpr bool kUnderAsan = true;
#else
constexpr bool kUnderAsan = false;
#endif
#else
constexpr bool kUnderAsan = false;
#endif

#ifndef _WIN32
TEST_CASE("C API: absurd reserve size returns 0 (length_error caught), never aborts",
          "[kmer][capi][vuln-0003][oom]") {
    // reserve(2^62) exceeds any allocatable size -> ankerl throws (bad_alloc
    // or length_error) before any allocation; the C boundary must translate
    // it to return 0 + diagnostic, never an exception crossing into C.
    // Skipped under ASan: ASan's allocator aborts on failed huge
    // allocations before the C++ exception can be thrown, so the boundary
    // path cannot be exercised there (verified under Release + UBSan).
    if (kUnderAsan) {
        WARN("skipping absurd-size reserve under AddressSanitizer (allocator aborts on OOM before bad_alloc)");
        return;
    }
    kmerindex_t* h = kmerindex_create();
    REQUIRE(h != nullptr);
    REQUIRE(kmerindex_reserve(h, 1ull << 62) == 0);
    REQUIRE(std::string(kmerindex_last_error(h)).find("allocation") != std::string::npos);
    // Handle still usable after the failed reserve.
    REQUIRE(kmerindex_size(h) == 0);
    REQUIRE(kmerindex_insert(h, 42, 1) == 1);
    uint64_t matched = 0;
    const uint64_t* v = kmerindex_get(h, 42, &matched);
    REQUIRE(v != nullptr);
    REQUIRE(*v == 1);
    kmerindex_destroy(h);
}

TEST_CASE("C API: huge reserve under memory pressure returns 0, never aborts",
          "[kmer][capi][vuln-0003][oom]") {
    // Real memory-pressure path: lower RLIMIT_AS to 256 MiB, then reserve
    // 2^34 entries (~64 GiB of buckets) -> std::bad_alloc inside reserve().
    // If the exception leaked across the C boundary the process would abort
    // (SIGABRT) and this test would fail loudly; expected is a clean 0
    // return with a diagnostic. Skipped under ASan (rlimit breaks ASan's
    // shadow mmaps) -- the length_error path above covers the boundary.
    if (kUnderAsan) {
        WARN("skipping RLIMIT_AS trigger under AddressSanitizer (shadow mmap conflict)");
        return;
    }
    // Construct the handle BEFORE constraining the process.
    kmerindex_t* h = kmerindex_create();
    REQUIRE(h != nullptr);

    RlimitAsGuard guard(256ull * 1024 * 1024);
    if (guard.applied) {
        REQUIRE(kmerindex_reserve(h, 1ull << 34) == 0);
        REQUIRE(std::string(kmerindex_last_error(h)).find("allocation") != std::string::npos);
        // Handle still usable after the failed reserve.
        REQUIRE(kmerindex_size(h) == 0);
        REQUIRE(kmerindex_insert(h, 42, 1) == 1);
        uint64_t matched = 0;
        const uint64_t* v = kmerindex_get(h, 42, &matched);
        REQUIRE(v != nullptr);
        REQUIRE(*v == 1);
    } else {
        WARN("setrlimit(RLIMIT_AS) not permitted in this environment — skipping OOM trigger check");
    }
    kmerindex_destroy(h);
}
#endif

// ── Null validation ─────────────────────────────────────────────────────────

TEST_CASE("C API: NULL handle/output pointers never crash", "[kmer][capi][null]") {
    kmerindex_t* h = kmerindex_create();
    REQUIRE(h != nullptr);

    REQUIRE(kmerindex_size(nullptr) == 0);
    REQUIRE(std::string(kmerindex_last_error(h)).find("null") != std::string::npos);

    REQUIRE(kmerindex_reserve(nullptr, 16) == 0);
    REQUIRE(std::string(kmerindex_last_error(h)).find("null") != std::string::npos);

    REQUIRE(kmerindex_insert(nullptr, 1, 1) == 0); // 0 = failure value, no crash
    uint64_t matched = 0;
    REQUIRE(kmerindex_get(nullptr, 1, &matched) == nullptr);
    REQUIRE(kmerindex_get(h, 1, nullptr) == nullptr);
    REQUIRE(std::string(kmerindex_last_error(h)).find("null") != std::string::npos);

    REQUIRE(kmerindex_freeze(nullptr) == 0);
    REQUIRE(kmerindex_iter_begin(nullptr, nullptr) == 0);
    REQUIRE(kmerindex_iter_begin(h, nullptr) == 0);

    kmerindex_iter_t it{};
    REQUIRE(kmerindex_iter_next(nullptr, &matched, &matched) == 0);
    REQUIRE(kmerindex_iter_next(&it, nullptr, &matched) == 0);
    REQUIRE(kmerindex_iter_next(&it, &matched, nullptr) == 0);

    // last_error is never NULL even with a NULL handle.
    REQUIRE(kmerindex_last_error(nullptr) != nullptr);

    kmerindex_destroy(h);
}

// ── Freeze / pointer stability ──────────────────────────────────────────────

TEST_CASE("C API: freeze enforcement — get/iter/freeze block later inserts",
          "[kmer][capi][freeze]") {
    // (a) get() freezes: interior pointers stay stable.
    kmerindex_t* h = kmerindex_create();
    REQUIRE(h != nullptr);
    REQUIRE(kmerindex_insert(h, 10, 1) == 1);
    uint64_t matched = 0;
    const uint64_t* p = kmerindex_get(h, 10, &matched);
    REQUIRE(p != nullptr);
    REQUIRE(*p == 1);

    REQUIRE(kmerindex_insert(h, 11, 2) == 0); // frozen
    REQUIRE(std::string(kmerindex_last_error(h)).find("frozen") != std::string::npos);
    REQUIRE(kmerindex_reserve(h, 1000) == 0); // frozen
    REQUIRE(std::string(kmerindex_last_error(h)).find("frozen") != std::string::npos);
    REQUIRE(kmerindex_size(h) == 1); // no growth happened

    // Reads still work after freezing.
    REQUIRE(kmerindex_get(h, 10, &matched) == p);
    REQUIRE(*p == 1);
    kmerindex_destroy(h);

    // (b) iter_begin() freezes.
    kmerindex_t* h2 = kmerindex_create();
    REQUIRE(h2 != nullptr);
    REQUIRE(kmerindex_insert(h2, 10, 1) == 1);
    kmerindex_iter_t it{};
    REQUIRE(kmerindex_iter_begin(h2, &it) == 1);
    REQUIRE(kmerindex_insert(h2, 11, 2) == 0);
    REQUIRE(std::string(kmerindex_last_error(h2)).find("frozen") != std::string::npos);
    kmerindex_destroy(h2);

    // (c) explicit freeze().
    kmerindex_t* h3 = kmerindex_create();
    REQUIRE(h3 != nullptr);
    REQUIRE(kmerindex_freeze(h3) == 1);
    REQUIRE(kmerindex_insert(h3, 1, 1) == 0);
    REQUIRE(std::string(kmerindex_last_error(h3)).find("frozen") != std::string::npos);
    REQUIRE(kmerindex_reserve(h3, 8) == 0);
    kmerindex_destroy(h3);

    // (d) get() on a MISS also freezes (read phase has begun).
    kmerindex_t* h4 = kmerindex_create();
    REQUIRE(h4 != nullptr);
    REQUIRE(kmerindex_get(h4, 999, &matched) == nullptr);
    REQUIRE(kmerindex_insert(h4, 1, 1) == 0);
    REQUIRE(std::string(kmerindex_last_error(h4)).find("frozen") != std::string::npos);
    kmerindex_destroy(h4);
}

// ── Caller-owned iterator: correctness + reentrancy ─────────────────────────

TEST_CASE("C API: caller-owned iterators are independent and reentrant",
          "[kmer][capi][iter]") {
    const uint64_t n = 100;
    kmerindex_t* h = kmerindex_create();
    REQUIRE(h != nullptr);
    REQUIRE(kmerindex_reserve(h, n) == 1);
    std::set<uint64_t> keys;
    for (uint64_t i = 1; i <= n; ++i) {
        REQUIRE(kmerindex_insert(h, i * 2, i * 7) == 1);
        keys.insert(i * 2);
    }

    // Two iterators over the same index, interleaved step by step: each
    // must independently visit all n entries (no shared cursor).
    kmerindex_iter_t it1{};
    kmerindex_iter_t it2{};
    REQUIRE(kmerindex_iter_begin(h, &it1) == 1);
    REQUIRE(kmerindex_iter_begin(h, &it2) == 1);

    std::set<uint64_t> seen1, seen2;
    uint64_t k1 = 0, v1 = 0, k2 = 0, v2 = 0;
    for (;;) {
        int r1 = kmerindex_iter_next(&it1, &k1, &v1);
        int r2 = kmerindex_iter_next(&it2, &k2, &v2);
        if (!r1 && !r2) break;
        if (r1) { seen1.insert(k1); REQUIRE(v1 == k1 * 7 / 2); }
        if (r2) { seen2.insert(k2); REQUIRE(v2 == k2 * 7 / 2); }
    }
    REQUIRE(seen1 == keys);
    REQUIRE(seen2 == keys);

    // A third iterator started after the first two finished also sees all.
    kmerindex_iter_t it3{};
    REQUIRE(kmerindex_iter_begin(h, &it3) == 1);
    std::set<uint64_t> seen3;
    while (kmerindex_iter_next(&it3, &k1, &v1)) seen3.insert(k1);
    REQUIRE(seen3 == keys);

    // An uninitialized iterator fails cleanly instead of crashing.
    kmerindex_iter_t uninit{};
    REQUIRE(kmerindex_iter_next(&uninit, &k1, &v1) == 0);
    REQUIRE(std::string(kmerindex_last_error(h)).find("not initialized") != std::string::npos);

    kmerindex_destroy(h);
}

// ── Concurrent freeze/get (minimap2-style multithreaded lookups) ───────────

TEST_CASE("C API: concurrent get after freeze is race-free; inserts fail",
          "[kmer][capi][thread][freeze]") {
    // Minimap2 runs lookups from mapping worker threads while every get()
    // WRITES the shim's freeze flag and a straggler insert() READS it --
    // before the fix that flag was a plain bool, so concurrent read/write
    // was UB (the design-review defect). After an explicit freeze the gets
    // only re-store the (already-true) flag and every insert must fail
    // before touching the map, so the only shared state touched
    // concurrently is the atomic flag itself. This test pins that contract:
    // no crash, no hang, no data race (see the TSan harness run).
    constexpr uint64_t kNumKeys = 100'000;
    kmerindex_t* h = kmerindex_create();
    REQUIRE(h != nullptr);
    REQUIRE(kmerindex_reserve(h, kNumKeys) == 1);
    for (uint64_t i = 0; i < kNumKeys; ++i) {
        REQUIRE(kmerindex_insert(h, (i << 1) | 1u, i * 3) == 1);
    }
    REQUIRE(kmerindex_size(h) == kNumKeys);
    REQUIRE(kmerindex_freeze(h) == 1);

    constexpr int kNumReaders = 8;
    constexpr uint64_t kReadsPerThread = 25'000;
    constexpr uint64_t kInsertAttempts = 50'000;

    std::atomic<bool> read_violation{false};
    std::atomic<uint64_t> hits{0}, misses{0};

    std::vector<std::thread> readers;
    readers.reserve(kNumReaders);
    for (int t = 0; t < kNumReaders; ++t) {
        readers.emplace_back([&, t] {
            const uint64_t stride = 2'654'435'761ull; // odd => wide base spread
            for (uint64_t r = 0;
                 r < kReadsPerThread && !read_violation.load(std::memory_order_relaxed);
                 ++r) {
                const uint64_t base = (r * 7 + (uint64_t)t * stride) % (2 * kNumKeys);
                uint64_t matched = 0;
                const uint64_t* v = kmerindex_get(h, base << 1, &matched);
                if (base < kNumKeys) {
                    // Guaranteed hit: stored as (base<<1)|1, value base*3;
                    // lookups use bit0=0 (minimap2 ignore-bit0 equality).
                    if (v == nullptr || *v != base * 3 || matched != ((base << 1) | 1u)) {
                        read_violation.store(true, std::memory_order_relaxed);
                        return;
                    }
                    ++hits;
                } else {
                    if (v != nullptr) {
                        read_violation.store(true, std::memory_order_relaxed);
                        return;
                    }
                    ++misses;
                }
            }
        });
    }

    // Main thread: every insert attempt must fail (return 0) -- never
    // succeed (1) and never hit the exception path (-1). Keys are unique
    // and absent, so the only legal 0 is the frozen rejection; the
    // thread-local diagnostic must read "frozen" on THIS thread.
    for (uint64_t i = 0; i < kInsertAttempts; ++i) {
        REQUIRE(kmerindex_insert(h, (kNumKeys + i) << 1, i) == 0);
    }
    REQUIRE(std::string(kmerindex_last_error(h)).find("frozen") != std::string::npos);

    for (auto& th : readers) th.join();
    REQUIRE_FALSE(read_violation.load());
    REQUIRE(hits.load() + misses.load() == kNumReaders * kReadsPerThread);
    REQUIRE(hits.load() > 0);  // mixed workload: guaranteed hits happened
    REQUIRE(misses.load() > 0); // ...and guaranteed misses happened

    // Reads still valid after the storm: interior pointers stable.
    for (uint64_t i = 0; i < kNumKeys; i += 9973) {
        uint64_t matched = 0;
        const uint64_t* v = kmerindex_get(h, i << 1, &matched);
        REQUIRE(v != nullptr);
        REQUIRE(*v == i * 3);
        REQUIRE(matched == ((i << 1) | 1u));
    }
    kmerindex_destroy(h);
}

// ── KmerReferenceIndex: malformed TRKI rejection (vuln-0006) ────────────────

TEST_CASE("TRKI: crafted header with n_slots=2^60 is rejected (no overflow)",
          "[kmer][trki][vuln-0006][malformed]") {
    // vuln-0006 case (a): n_slots=2^60 wraps n_slots*16 to 0 in unchecked
    // arithmetic; the old code accepted a 56-byte file and then built a
    // wild pointer. The new loader must reject it.
    const std::string path = "trki_evil_a.trki";
    write_trki(path, 1, 4, 1ull << 60, 1, 2, {}, {0, 0});
    REQUIRE(fs::file_size(path) == 56);
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(path), std::runtime_error);
    fs::remove(path);
}

TEST_CASE("TRKI: slot with out-of-range offset/count is rejected at load",
          "[kmer][trki][vuln-0006][malformed]") {
    // vuln-0006 case (b): value=(1,000,000<<32)|4 with n_positions=1 used
    // to build an ~8 MB out-of-bounds span. Now rejected during the load
    // time slot walk; even if it slipped through, query() re-validates.
    const std::string path = "trki_evil_b.trki";
    write_trki(path, 1, 4, 1, 1, 1,
               {{42, (1'000'000ull << 32) | 4u}}, {0});
    REQUIRE(fs::file_size(path) == 64);
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(path), std::runtime_error);
    fs::remove(path);
}

TEST_CASE("TRKI: truncated and oversized files are rejected",
          "[kmer][trki][vuln-0006][malformed]") {
    // Declared tables larger than the actual file (truncated).
    const std::string trunc = "trki_truncated.trki";
    write_trki(trunc, 1, 4, 4, 1, 1, {}, {0}); // claims 4 slots but writes none
    REQUIRE(fs::file_size(trunc) == 40 + 8);   // header + 1 position only
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(trunc), std::runtime_error);
    fs::remove(trunc);

    // Declared tables smaller than the actual file (trailing garbage).
    const std::string garbage = "trki_trailing.trki";
    write_trki(garbage, 1, 4, 1, 1, 1, {{42, 0}}, {0});
    append_bytes(garbage, "JUNKJUNK", 8);
    REQUIRE(fs::file_size(garbage) == 40 + 16 + 8 + 8);
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(garbage), std::runtime_error);
    fs::remove(garbage);

    // File too small to hold even a header.
    const std::string tiny = "trki_tiny.trki";
    write_trki(tiny, 1, 4, 1, 1, 1, {{42, 0}}, {0});
    {
        std::ofstream f(tiny, std::ios::binary | std::ios::trunc);
        f.write("TRKI", 4); // 4 bytes only
    }
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(tiny), std::runtime_error);
    fs::remove(tiny);
}

TEST_CASE("TRKI: bad magic / version / k / slot-count / n_keys rejected",
          "[kmer][trki][vuln-0006][malformed]") {
    const std::string p = "trki_badfields.trki";

    write_trki(p, 1, 4, 1, 1, 1, {{42, 0}}, {0}, /*bad_magic=*/true);
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(p), std::runtime_error);

    write_trki(p, 2, 4, 1, 1, 1, {{42, 0}}, {0}); // version 2 unsupported
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(p), std::runtime_error);

    write_trki(p, 1, 0, 1, 1, 1, {{42, 0}}, {0}); // k = 0
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(p), std::runtime_error);
    write_trki(p, 1, 33, 1, 1, 1, {{42, 0}}, {0}); // k = 33 > packing limit
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(p), std::runtime_error);

    write_trki(p, 1, 4, 3, 3, 3, {{1, 0}, {2, 0}, {3, 0}}, {0, 1, 2}); // n_slots=3 not pow2
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(p), std::runtime_error);

    write_trki(p, 1, 4, 2, 5, 1, {{1, 0}, {2, 0}}, {0}); // n_keys > n_slots
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(p), std::runtime_error);

    // slot value 0 encodes count=0 -> "slot references out-of-range positions"
    write_trki(p, 1, 4, 1, 1, 1, {{42, 0}}, {0});
    REQUIRE_THROWS_AS(TracEon::KmerReferenceIndex::mmap_open(p), std::runtime_error);

    fs::remove(p);
}

TEST_CASE("TRKI: empty index (40-byte file) loads and queries empty",
          "[kmer][trki]") {
    const std::string path = "trki_empty.trki";
    write_trki(path, 1, 4, 0, 0, 0);
    REQUIRE(fs::file_size(path) == 40);
    auto idx = TracEon::KmerReferenceIndex::mmap_open(path);
    REQUIRE(idx.n_keys() == 0);
    REQUIRE(idx.query(27).empty());
    fs::remove(path);
}

// ── KmerReferenceIndex: valid round-trip ────────────────────────────────────

TEST_CASE("TRKI: build -> save -> mmap_open round-trip preserves queries",
          "[kmer][trki][roundtrip]") {
    const std::string seq = "ACGTACGTACGTACGTACGT"; // 20 bases, k=4
    const std::string path = "trki_roundtrip.trki";

    TracEon::KmerReferenceIndex idx;
    idx.build_from_sequence(seq, 4);
    REQUIRE(idx.n_keys() == 4); // ACGT, CGTA, GTAC, TACG
    REQUIRE(idx.k() == 4);

    // Pre-save query snapshot (ACGT=27 -> positions 0,4,8,12,16).
    const uint64_t acgt = 27;
    auto pre = idx.query(acgt);
    REQUIRE(pre.size() == 5);
    const std::vector<uint64_t> pre_pos(pre.begin(), pre.end());
    REQUIRE(pre_pos == std::vector<uint64_t>({0, 4, 8, 12, 16}));
    REQUIRE(idx.query(0).empty()); // AAAA never occurs

    idx.save(path);
    // 20 bases / k=4 -> 17 windows: ACGT x5, CGTA x4, GTAC x4, TACG x4.
    // File layout: 40-byte header + n_slots(16)*16 + n_positions(17)*8 = 432 bytes.
    REQUIRE(fs::file_size(path) == 40 + 16 * 16 + 17 * 8);

    auto loaded = TracEon::KmerReferenceIndex::mmap_open(path);
    REQUIRE(loaded.k() == 4);
    REQUIRE(loaded.n_keys() == 4);

    auto post = loaded.query(acgt);
    REQUIRE(post.size() == 5);
    const std::vector<uint64_t> post_pos(post.begin(), post.end());
    REQUIRE(post_pos == pre_pos);

    // Every distinct k-mer round-trips identically.
    for (int i = 0; i + 4 <= (int)seq.size(); ++i) {
        uint64_t code = 0;
        REQUIRE(TracEon::encode_kmer(std::string_view(seq).substr(i, 4), 4, code));
        auto a = idx.query(code);
        auto b = loaded.query(code);
        REQUIRE(std::vector<uint64_t>(a.begin(), a.end()) ==
                std::vector<uint64_t>(b.begin(), b.end()));
    }

    fs::remove(path);
}

TEST_CASE("TRKI: build_from_sequence with invalid k yields an empty index",
          "[kmer][trki]") {
    TracEon::KmerReferenceIndex idx;
    idx.build_from_sequence("ACGTACGT", 0);   // k=0
    REQUIRE(idx.n_keys() == 0);
    REQUIRE(idx.query(27).empty());
    idx.build_from_sequence("ACGTACGT", 33);  // k > 32
    REQUIRE(idx.n_keys() == 0);
    idx.build_from_sequence("ACGTACGT", 100); // k > seq size
    REQUIRE(idx.n_keys() == 0);
    REQUIRE(idx.query(0).empty());
}
