// K-mer experimental target tests (kmer_tests).
//
// Covers: encode_kmer guards (vuln-0001) and the C API happy path + error
// model (vuln-0003: no exception can cross the C boundary), null-handle
// validation, and freeze/pointer-stability enforcement.

#include <catch2/catch_test_macros.hpp>

#include "KmerEncoding.h"
#include "kmerindex_c_api.h"

#include <cstdint>
#include <thread>
#include <atomic>
#include <vector>
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

TEST_CASE("C API happy path: create/reserve/insert/get/destroy",
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
    // reserve(2^62) is far above the container's max_size() (~2^32 entries
    // for this uint32 value_idx map). ankerl::unordered_dense 4.4.0
    // reserve() CLAMPS to max_size() instead of throwing (bad_alloc or
    // length_error), so without a pre-check this would SILENTLY SUCCEED
    // (return 1) on machines where the huge allocation lands -- the bug
    // fixed in GitHub issue #9. The C boundary now pre-checks
    // n > map.max_size() and rejects such requests deterministically,
    // returning 0 with the "requested size too large" diagnostic BEFORE
    // any allocation; the result is platform-independent.
    // No ASan skip is needed for this test: the pre-check fires BEFORE any
    // allocation, so no huge allocation is ever attempted and the boundary
    // path works identically under AddressSanitizer.
    kmerindex_t* h = kmerindex_create();
    REQUIRE(h != nullptr);
    REQUIRE(kmerindex_reserve(h, 1ull << 62) == 0);
    // The C boundary's pre-check rejects the size BEFORE any allocation,
    // so the diagnostic is "requested size too large" -- asserting
    // "allocation" would encode the old container-throwing behavior.
    REQUIRE(std::string(kmerindex_last_error(h)).find("requested size too large") != std::string::npos);
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
    // 2^30 entries (~24-32 GB of buckets). The size stays BELOW the
    // container's max_size() (~2^32), so it passes the deterministic
    // pre-check and genuinely exercises reserve()'s allocation failure:
    // the ~24-32 GB it requests is far above the 256 MiB rlimit, forcing
    // std::bad_alloc. (A size ABOVE max_size() would be rejected by the
    // pre-check instead -- that is the other [oom] test's path.) If the
    // exception leaked across the C boundary the process would abort
    // (SIGABRT) and this test would fail loudly; expected is a clean 0
    // return with a diagnostic. Skipped under ASan (rlimit breaks ASan's
    // shadow mmaps) -- the absurd-size pre-check path above covers the
    // boundary.
    if (kUnderAsan) {
        WARN("skipping RLIMIT_AS trigger under AddressSanitizer (shadow mmap conflict)");
        return;
    }
    // Construct the handle BEFORE constraining the process.
    kmerindex_t* h = kmerindex_create();
    REQUIRE(h != nullptr);

    RlimitAsGuard guard(256ull * 1024 * 1024);
    if (guard.applied) {
        REQUIRE(kmerindex_reserve(h, 1ull << 30) == 0);
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

    // last_error is never NULL even with a NULL handle.
    REQUIRE(kmerindex_last_error(nullptr) != nullptr);

    kmerindex_destroy(h);
}

// ── Freeze / pointer stability ──────────────────────────────────────────────

TEST_CASE("C API: get() freezes the index and blocks later inserts",
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

    // (b) get() on a MISS also freezes (read phase has begun).
    kmerindex_t* h4 = kmerindex_create();
    REQUIRE(h4 != nullptr);
    REQUIRE(kmerindex_get(h4, 999, &matched) == nullptr);
    REQUIRE(kmerindex_insert(h4, 1, 1) == 0);
    REQUIRE(std::string(kmerindex_last_error(h4)).find("frozen") != std::string::npos);
    kmerindex_destroy(h4);
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
    uint64_t freeze_key = 0;
    REQUIRE(kmerindex_get(h, 0, &freeze_key) != nullptr); // freeze via get()

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
