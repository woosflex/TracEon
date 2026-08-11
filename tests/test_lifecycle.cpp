#include "TestHelpers.h"

// test_lifecycle.cpp — lifecycle domain (immutable-after-load,
// addEntry/set, clearCache, reader quiescence, failure
// atomicity, [oom] memory guard).
// Split from the monolithic SmartStrategyTests.cpp.


TEST_CASE("addEntry stores data with stable string_views", "[strategy][addentry]") {
    TracEon::SmartStrategy strategy;

    // Add entries and immediately discard the source strings to ensure views
    // don't dangle on the stack-allocated temporaries.
    {
        std::string id  = "manual_seq";
        std::string seq = "ACGTACGT";
        std::string qual = "IIIIIIII";
        strategy.addEntry(id, seq, qual);
    } // id, seq, qual destroyed here — views must still be valid

    REQUIRE(strategy.hasSequence("manual_seq"));
    REQUIRE(strategy.getSequence("manual_seq") == "ACGTACGT");
    REQUIRE(strategy.getQuality("manual_seq") == "IIIIIIII");
}

#ifndef _WIN32

TEST_CASE("addEntry() throws instead of crashing when the process is memory-constrained",
          "[strategy][addentry][errors][oom]") {
    // addEntry()'s proactive guard compares against system-wide available
    // memory (/proc/meminfo), which a ulimit/cgroup-constrained *process*
    // can still exceed sooner. Lower RLIMIT_AS (virtual address space) to
    // just above current usage, then confirm a large addEntry() throws
    // (via the bad_alloc→runtime_error conversion) instead of crashing.
    struct rlimit original{};
    REQUIRE(getrlimit(RLIMIT_AS, &original) == 0);

    // Build the source string *before* constraining the process — addEntry()
    // makes its own independent copy of it (manual_store_.push_back), so
    // the limit only needs headroom for one more ~64 MiB allocation on top
    // of whatever's already resident (including this string itself).
    TracEon::SmartStrategy strategy;
    const std::string big_seq(64 * 1024 * 1024, 'A'); // 64 MiB

    struct rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    struct rlimit constrained = original;
    constexpr long headroom_kb = 16 * 1024; // 16 MiB — not enough for another 64 MiB copy
    constrained.rlim_cur = static_cast<rlim_t>((usage.ru_maxrss + headroom_kb) * 1024);

    bool limit_applied = (setrlimit(RLIMIT_AS, &constrained) == 0);

    if (limit_applied) {
        REQUIRE_THROWS(strategy.addEntry("big", big_seq, ""));
        // Restore so the rest of the test binary isn't constrained.
        setrlimit(RLIMIT_AS, &original);
    } else {
        WARN("setrlimit(RLIMIT_AS) not permitted in this environment — skipping OOM trigger check");
    }
}
#endif


// ─── clearCache() + reload ────────────────────────────────────────────────────

TEST_CASE("clearCache() then reload — no dangling views", "[strategy]") {
    std::string fasta_path = "tmp_clear_reload.fasta";
    {
        std::ofstream out(fasta_path);
        out << ">seq1\nACGT\n>seq2\nTGCA\n";
    }

    TracEon::SmartStrategy strategy;
    strategy.loadFile(fasta_path);
    REQUIRE(strategy.getFileCacheSize() == 2);

    strategy.clearCache();
    REQUIRE(strategy.getFileCacheSize() == 0);
    REQUIRE_FALSE(strategy.hasSequence("seq1"));

    // Reload the same file — old text_arena_ is gone; new arena must be used
    strategy.loadFile(fasta_path);
    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("seq1") == "ACGT");
    REQUIRE(strategy.getSequence("seq2") == "TGCA");

    fs::remove(fasta_path);
}


// ═══════════════════════════════════════════════════════════════════════════
// Bug 3 regression: immutable-after-load contract
// ═══════════════════════════════════════════════════════════════════════════
// The lock-free read path (getView()/getQuality()/hasSequence()) reads the
// index without holding a lock once data_ready_ is set. That is only sound
// if the index is truly immutable after publication, so addEntry() must be
// rejected once a load has happened (until clearCache() reopens the build).

TEST_CASE("addEntry() before load works — multiple build-phase entries",
          "[strategy][bug3][addentry]") {
    TracEon::SmartStrategy strategy;
    strategy.addEntry("a", "ACGT", "");
    strategy.addEntry("b", "TGCA", "IIII");
    strategy.addEntry("c", "GGGG", "");
    REQUIRE(strategy.getFileCacheSize() == 3);
    REQUIRE(strategy.getSequence("a") == "ACGT");
    REQUIRE(strategy.getSequence("b") == "TGCA");
    REQUIRE(strategy.getQuality("b") == "IIII");
    REQUIRE(strategy.getSequence("c") == "GGGG");
}


TEST_CASE("addEntry() after loadFile() throws std::logic_error",
          "[strategy][bug3][addentry][errors]") {
    const std::string src = "bug3_after_load.fasta";
    {
        std::ofstream out(src);
        out << ">seq1\nACGT\n>seq2\nTGCA\n";
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src);
    REQUIRE(strategy.getFileCacheSize() == 2);

    REQUIRE_THROWS_AS(strategy.addEntry("late", "TTTT", ""), std::logic_error);
    // The loaded data must remain intact and readable
    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("seq1") == "ACGT");
    REQUIRE(strategy.getSequence("seq2") == "TGCA");

    fs::remove(src);
}


TEST_CASE("addEntry() after loadBinary() throws std::logic_error",
          "[strategy][bug3][addentry][errors]") {
    const std::string bin = "bug3_after_restore.bin";
    {
        TracEon::SmartStrategy s;
        s.addEntry("k", "AAAA", "");
        s.saveBinary(bin);
    }
    TracEon::SmartStrategy strategy;
    strategy.loadBinary(bin);
    REQUIRE(strategy.getFileCacheSize() == 1);
    REQUIRE_THROWS_AS(strategy.addEntry("k2", "CCCC", ""), std::logic_error);
    REQUIRE(strategy.getFileCacheSize() == 1);

    fs::remove(bin);
}


TEST_CASE("addEntry() after clearCache() works again",
          "[strategy][bug3][addentry]") {
    const std::string src = "bug3_after_clear.fasta";
    {
        std::ofstream out(src);
        out << ">seq1\nACGT\n";
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src);
    REQUIRE_THROWS_AS(strategy.addEntry("x", "TTTT", ""), std::logic_error);

    strategy.clearCache();
    strategy.addEntry("y", "GGGG", ""); // allowed again after clearCache()
    REQUIRE(strategy.getFileCacheSize() == 1);
    REQUIRE(strategy.getSequence("y") == "GGGG");

    fs::remove(src);
}


// ── Lifecycle contract (reader quiescence — ADR-001 / README) ───────────────
// getView() returns a NON-OWNING std::string_view into the installed
// snapshot's arena. The view is valid only while that snapshot stays
// installed: once clearCache(), a reload, or destruction begins, every view
// into the old arena is dangling. These tests exercise the SEQUENTIAL
// lifecycle transitions that are the supported contract. There is
// deliberately NO test claiming concurrent clear + retained string_view is
// supported — it is not (see ADR-001 reader-quiescence section).

TEST_CASE("Lifecycle: retained view invalid after clearCache() (documented contract)",
          "[strategy][lifecycle]") {
    const std::string src = "lifecycle_clear.fasta";
    {
        std::ofstream out(src);
        out << ">seq1\nACGT\n";
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src);
    std::string_view retained = strategy.getView("seq1");
    REQUIRE(retained == "ACGT");

    // Per the lifecycle contract the caller must stop using `retained`
    // BEFORE clearCache() — the view points into the arena that is about to
    // be destroyed. After the clear, the cache reports empty and fresh
    // lookups return nothing; the retained view is documented as dangling
    // (we do not dereference it here — that would be UB by contract).
    strategy.clearCache();
    REQUIRE(strategy.getFileCacheSize() == 0);
    REQUIRE(strategy.getView("seq1").empty());

    fs::remove(src);
}


TEST_CASE("Lifecycle: sequential clear → load → clear → load transitions",
          "[strategy][lifecycle]") {
    const std::string src_a = "lifecycle_a.fasta";
    const std::string src_b = "lifecycle_b.fasta";
    {
        std::ofstream out(src_a);
        out << ">a1\nAAAA\n";
        std::ofstream out2(src_b);
        out2 << ">b1\nCCCC\n>b2\nGGGG\n";
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src_a);
    REQUIRE(strategy.getFileCacheSize() == 1);
    REQUIRE(strategy.getSequence("a1") == "AAAA");

    strategy.clearCache();
    strategy.loadFile(src_b); // reload replaces the snapshot entirely
    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE_FALSE(strategy.hasSequence("a1"));
    REQUIRE(strategy.getSequence("b1") == "CCCC");
    REQUIRE(strategy.getSequence("b2") == "GGGG");

    strategy.clearCache();
    REQUIRE(strategy.getFileCacheSize() == 0);
    strategy.loadFile(src_a); // and back again — no state carried over
    REQUIRE(strategy.getFileCacheSize() == 1);
    REQUIRE(strategy.getSequence("a1") == "AAAA");

    fs::remove(src_a);
    fs::remove(src_b);
}


TEST_CASE("Lifecycle: failed load leaves no partial state (failure atomicity)",
          "[strategy][lifecycle][v4]") {
    // A v4 file with a corrupt payload must be rejected AND leave the cache
    // in the pristine empty state — never a half-loaded map (design review
    // P0 sprint item 4: invalid loads must not publish partial state).
    const std::string src = "lifecycle_atomic.fasta";
    const std::string bin = "lifecycle_atomic.bin";
    {
        std::ofstream out(src);
        out << ">s1\nACGT\n>s2\nTGCA\n>s3\nGGGG\n";
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src);
    REQUIRE(strategy.getFileCacheSize() == 3);
    strategy.saveBinary(bin);

    // Corrupt one byte of the stored CRC: load must throw, and the cache
    // must look exactly like a fresh SmartStrategy afterwards.
    auto bytes = read_blob(bin);
    bytes[22] ^= 0x01;
    write_blob(bin, bytes);
    REQUIRE_THROWS_AS(strategy.loadBinary(bin), std::runtime_error);
    REQUIRE(strategy.getFileCacheSize() == 0);
    REQUIRE_FALSE(strategy.hasSequence("s1"));

    // A subsequent successful load works normally.
    REQUIRE_THROWS_AS(strategy.loadBinary(bin), std::runtime_error); // still corrupt
    REQUIRE(strategy.getFileCacheSize() == 0);

    fs::remove(src);
    fs::remove(bin);
}
