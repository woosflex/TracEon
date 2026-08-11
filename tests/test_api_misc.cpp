#include "TestHelpers.h"

// test_api_misc.cpp — everything else (initial state, getView,
// size, keys, duplicate keys, SIMD, NGS mode, format
// detection, save/restore basics).
// Split from the monolithic SmartStrategyTests.cpp.


TEST_CASE("SmartStrategy initial state", "[SmartStrategy]") {
    TracEon::SmartStrategy strategy;
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::UNKNOWN);
    REQUIRE(strategy.getFileCacheSize() == 0);
}


TEST_CASE("Lock-Free Read Performance", "[strategy][performance]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "simple.fasta";
    {
        std::ofstream out(test_file_path);
        out << ">seq1\nACGT\n>seq2\nTGCA\n";
        out.close();
    }
    strategy.loadFile(test_file_path);
    
    SECTION("Can perform concurrent reads without blocking") {
        const int NUM_THREADS = 4;
        const int LOOKUPS_PER_THREAD = 10000;
        std::vector<std::thread> threads;
        std::atomic<size_t> total_hits{0};
        
        auto worker = [&]() {
            size_t hits = 0;
            for (int i = 0; i < LOOKUPS_PER_THREAD; ++i) {
                auto view = strategy.getView("seq1");
                if (!view.empty()) hits++;
            }
            total_hits += hits;
        };
        
        for (int i = 0; i < NUM_THREADS; ++i) threads.emplace_back(worker);
        for (auto& t : threads) t.join();
        REQUIRE(total_hits == NUM_THREADS * LOOKUPS_PER_THREAD);
    }
    fs::remove(test_file_path);
}


TEST_CASE("String Key Storage Validation", "[strategy][architecture]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "key_test.fasta";
    {
        std::ofstream out(test_file_path);
        out << ">seq1\nACGT\n";
        out.close();
    }
    strategy.loadFile(test_file_path);
    REQUIRE(strategy.getIndexMode() == TracEon::IndexMode::GENOME);
    auto keys = strategy.getAllKeys();
    REQUIRE(keys.size() == 1);
    REQUIRE(keys[0] == "seq1");
    fs::remove(test_file_path);
}


TEST_CASE("loadFile throws on missing file", "[strategy][errors]") {
    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadFile("nonexistent_file_xyz.fasta"), std::runtime_error);
}


// ── simd_find_char unit tests ─────────────────────────────────────────────────
// These tests verify correctness of the SIMD and scalar paths regardless of
// which ISA extension is active at compile time.  The cases are designed so
// that both the "full SIMD lane" path and the scalar tail path are exercised:
//   AVX2  → full loop: strings ≥32 bytes; tail: strings 1–31 bytes
//   NEON  → full loop: strings ≥16 bytes; tail: strings 1–15 bytes
//   memchr → all sizes covered by the libc implementation

TEST_CASE("simd_find_char — basic correctness", "[simd]") {
    using TracEon::simd_find_char;

    // 64 bytes: exercises at least one full SIMD lane on both AVX2 and NEON.
    // The target char '>' is placed at index 32 (just past one AVX2 lane).
    const std::string data(32, 'A');                     // 32 × 'A'
    const std::string rest = ">" + std::string(31, 'C'); // '>' then 31 × 'C'
    const std::string buf  = data + rest;                // 64 bytes total
    const char* b = buf.data();
    const char* e = b + buf.size();

    REQUIRE(simd_find_char(b, e, '>') == b + 32);  // found exactly at index 32
    REQUIRE(simd_find_char(b, e, 'A') == b);        // first byte matches
    REQUIRE(simd_find_char(b, e, 'C') == b + 33);  // first 'C' after the '>'
    REQUIRE(simd_find_char(b, e, 'Z') == e);        // not present → returns end
}


TEST_CASE("simd_find_char — scalar tail (< SIMD lane size)", "[simd]") {
    using TracEon::simd_find_char;

    // 7 bytes: always handled by the scalar tail on any SIMD width.
    const std::string buf = "ACGT\nAT";
    const char* b = buf.data();
    const char* e = b + buf.size();

    REQUIRE(simd_find_char(b, e, '\n') == b + 4);
    REQUIRE(simd_find_char(b, e, 'T')  == b + 3);  // first 'T' at index 3
    REQUIRE(simd_find_char(b, e, 'Z')  == e);
}


TEST_CASE("simd_find_char — empty range", "[simd]") {
    using TracEon::simd_find_char;
    const char* p = "X";
    REQUIRE(simd_find_char(p, p, 'X') == p);  // end == begin → never match
}


TEST_CASE("simd_find_char — first and last byte", "[simd]") {
    using TracEon::simd_find_char;

    // 33 bytes: straddles one AVX2 lane boundary.
    const std::string buf = ">" + std::string(31, 'N') + "<";
    const char* b = buf.data();
    const char* e = b + buf.size();

    REQUIRE(simd_find_char(b, e, '>') == b);       // first byte
    REQUIRE(simd_find_char(b, e, '<') == e - 1);   // last byte (scalar tail)
    REQUIRE(simd_find_char(b, e, 'N') == b + 1);   // second byte
}


// ── simd_find_char large-data cross-platform dispatch ──────────────────────────
// Verifies that simd_find_char dispatches correctly on the current platform
// by exercising large buffers that force full SIMD lane usage (not just the
// scalar tail).  On x86-64 this exercises the AVX2 path when available,
// on ARM64 the NEON path, and memchr on other platforms.

TEST_CASE("simd_find_char — large buffer dispatch (>=512 bytes)", "[simd][dispatch]") {
    using TracEon::simd_find_char;

    // 512 bytes: many SIMD lanes on any architecture (32 lanes of AVX2, 16 of NEON).
    std::string buf(512, 'A');
    buf[127] = '>';    // within first AVX2 lane group (bytes 96-127)
    buf[255] = '\n';   // at a lane boundary (byte 255 = 8*32 - 1 for AVX2)
    buf[383] = '>';    // well into the middle
    buf[511] = 'G';    // last byte

    const char* b = buf.data();
    const char* e = b + buf.size();

    REQUIRE(simd_find_char(b, e, '>') == b + 127);
    REQUIRE(simd_find_char(b, e, '\n') == b + 255);
    REQUIRE(simd_find_char(b, e, 'G') == b + 511);     // last byte
    REQUIRE(simd_find_char(b, e, 'Z') == e);            // not present
}


TEST_CASE("simd_find_char — multilane boundary scanning", "[simd][dispatch]") {
    using TracEon::simd_find_char;

    // 100 bytes: straddles AVX2 lane boundaries multiple times.
    std::string buf = std::string(31, 'A') + ">" + std::string(31, 'C') + "\n" + std::string(31, 'G');
    REQUIRE(buf.size() == 95);

    const char* b = buf.data();
    const char* e = b + buf.size();

    REQUIRE(simd_find_char(b, e, '>') == b + 31);   // at byte 31 (AVX2 lane boundary)
    REQUIRE(simd_find_char(b, e, '\n') == b + 63);  // at byte 63
    REQUIRE(simd_find_char(b, e, 'G') == b + 64);   // first G at byte 64
    REQUIRE(simd_find_char(b, e, 'Z') == e);
}


TEST_CASE("simd_find_char — search for '>' in FASTA-like content", "[simd][dispatch]") {
    using TracEon::simd_find_char;

    // Simulates what normalizeFastaArena/parseFastaInternal do:
    // scanning sequence data for '>' to find record boundaries.
    std::string seq_data(200, 'A');
    seq_data += ">next_record\n";
    seq_data += std::string(100, 'T');

    const char* b = seq_data.data();
    const char* e = b + seq_data.size();

    // '>' should be found exactly at offset 200
    REQUIRE(simd_find_char(b, e, '>') == b + 200);
    // '\n' should be at offset 212 (after ">next_record")
    REQUIRE(simd_find_char(b, e, '\n') == b + 212);
    // 'T' should be at offset 213
    REQUIRE(simd_find_char(b, e, 'T') == b + 213);
}


// ─── NGSIndex Mode ───────────────────────────────────────────────────────────

TEST_CASE("NGSIndex mode — wired up via constructor", "[strategy][ngs]") {
    TracEon::SmartStrategy strategy(TracEon::IndexMode::NGS);
    REQUIRE(strategy.getIndexMode() == TracEon::IndexMode::NGS);
}


TEST_CASE("NGSIndex mode — FASTA load and lookup", "[strategy][ngs]") {
    std::string fasta_path = "tmp_ngs_fasta.fasta";
    {
        std::ofstream out(fasta_path);
        out << ">read_A\nACGTACGT\n>read_B\nTTTTGGGG\n>read_C\nCCCCAAAA\n";
    }

    TracEon::SmartStrategy strategy(TracEon::IndexMode::NGS);
    strategy.loadFile(fasta_path);

    REQUIRE(strategy.getIndexMode() == TracEon::IndexMode::NGS);
    REQUIRE(strategy.getFileCacheSize() == 3);
    REQUIRE(strategy.hasSequence("read_A"));
    REQUIRE(strategy.getSequence("read_A") == "ACGTACGT");
    REQUIRE(strategy.hasSequence("read_B"));
    REQUIRE(strategy.getSequence("read_B") == "TTTTGGGG");
    REQUIRE(strategy.hasSequence("read_C"));
    REQUIRE(strategy.getSequence("read_C") == "CCCCAAAA");
    REQUIRE_FALSE(strategy.hasSequence("nonexistent"));

    std::vector<std::string> keys = strategy.getAllKeys();
    REQUIRE(keys.size() == 3);

    fs::remove(fasta_path);
}


TEST_CASE("NGSIndex mode — FASTQ load and lookup", "[strategy][ngs]") {
    std::string fastq_path = "tmp_ngs_fastq.fastq";
    {
        std::ofstream out(fastq_path);
        out << "@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nHHHH\n";
    }

    TracEon::SmartStrategy strategy(TracEon::IndexMode::NGS);
    strategy.loadFile(fastq_path);

    REQUIRE(strategy.getIndexMode() == TracEon::IndexMode::NGS);
    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("read1") == "ACGT");
    REQUIRE(strategy.getQuality("read1") == "IIII");
    REQUIRE(strategy.getSequence("read2") == "TGCA");
    REQUIRE(strategy.getQuality("read2") == "HHHH");

    fs::remove(fastq_path);
}


TEST_CASE("NGSIndex mode — save and restore round-trip", "[strategy][ngs]") {
    std::string fasta_path = "tmp_ngs_rtrip.fasta";
    std::string bin_path   = "tmp_ngs_rtrip.bin";
    {
        std::ofstream out(fasta_path);
        out << ">seq1\nACGT\n>seq2\nGGGG\n";
    }

    {
        TracEon::SmartStrategy s1(TracEon::IndexMode::NGS);
        s1.loadFile(fasta_path);
        s1.saveBinary(bin_path);
    }

    TracEon::SmartStrategy s2;
    s2.loadBinary(bin_path);

    REQUIRE(s2.getIndexMode() == TracEon::IndexMode::NGS);
    REQUIRE(s2.getFileCacheSize() == 2);
    REQUIRE(s2.getSequence("seq1") == "ACGT");
    REQUIRE(s2.getSequence("seq2") == "GGGG");

    fs::remove(fasta_path);
    fs::remove(bin_path);
}


TEST_CASE("NGSIndex mode — addEntry round-trip (immutable-after-load contract)", "[strategy][ngs]") {
    // Reworked for the immutable-after-load contract (Bug 3): addEntry() is
    // only legal before any load or after clearCache(), so the former
    // loadFile-then-addEntry combination is rejected by design. The NGS
    // load→save→restore path is covered by the "save and restore round-trip"
    // test above; this covers the NGS addEntry→save→restore path.
    std::string bin_path = "tmp_ngs_addentry.bin";

    {
        TracEon::SmartStrategy s1(TracEon::IndexMode::NGS);
        s1.addEntry("manual_seq", "TTTT", "");
        REQUIRE(s1.getFileCacheSize() == 1);
        REQUIRE(s1.getSequence("manual_seq") == "TTTT");
        s1.saveBinary(bin_path);
    }

    // Restore into a fresh default-mode strategy — loadBinary must set NGS mode
    TracEon::SmartStrategy s2;
    s2.loadBinary(bin_path);

    REQUIRE(s2.getIndexMode() == TracEon::IndexMode::NGS);
    REQUIRE(s2.getFileCacheSize() == 1);
    REQUIRE(s2.getSequence("manual_seq") == "TTTT");

    fs::remove(bin_path);
}


// ─── Format Detection ─────────────────────────────────────────────────────────

TEST_CASE("RNA FASTA format detection", "[strategy][format]") {
    std::string fasta_path = "tmp_rna.fasta";
    {
        std::ofstream out(fasta_path);
        out << ">transcript1\nAUGCUAUGCUAUGCUA\n>transcript2\nUUUUAAAACCCCGGGG\n";
    }

    TracEon::SmartStrategy strategy;
    strategy.loadFile(fasta_path);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::RNA_FASTA);

    fs::remove(fasta_path);
}


TEST_CASE("Protein FASTA format detection", "[strategy][format]") {
    std::string fasta_path = "tmp_protein.fasta";
    {
        std::ofstream out(fasta_path);
        out << ">prot1\nMKTLLLTLVVVTIVLAMGLSLSEEKE\n>prot2\nACDEFGHIKLMNPQRSTWYVACDE\n";
    }

    TracEon::SmartStrategy strategy;
    strategy.loadFile(fasta_path);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::PROTEIN_FASTA);

    fs::remove(fasta_path);
}


// ═══════════════════════════════════════════════════════════════════════════
// Bug 2 regression: duplicate-key set() must not corrupt save/restore
// ═══════════════════════════════════════════════════════════════════════════
// addEntry() used to call emplace() (a no-op on a duplicate key) while
// unconditionally incrementing serialized_size_estimate_, so the binary-cache
// header declared more bytes than the payload held → restore() threw
// "Binary cache v3 decompressed size mismatch". Semantics: FIRST value wins
// (emplace keeps the first insertion); duplicates are silent no-ops.

TEST_CASE("Duplicate key set() — first wins, save/restore consistent (GENOME)",
          "[strategy][bug2][addentry]") {
    const std::string bin = "bug2_dup.bin";
    {
        TracEon::SmartStrategy strategy;
        strategy.addEntry("dup", "AAAA", "IIII");
        strategy.addEntry("dup", "CCCC", "JJJJ"); // duplicate → no-op, first wins
        REQUIRE(strategy.getFileCacheSize() == 1);
        REQUIRE(strategy.getSequence("dup") == "AAAA");
        REQUIRE(strategy.getQuality("dup") == "IIII");

        strategy.saveBinary(bin); // header must not over-declare the payload
    }
    {
        TracEon::SmartStrategy restored;
        REQUIRE_NOTHROW(restored.loadBinary(bin)); // was: throws decompressed size mismatch
        REQUIRE(restored.getFileCacheSize() == 1);
        REQUIRE(restored.getSequence("dup") == "AAAA");
        REQUIRE(restored.getQuality("dup") == "IIII");
    }
    fs::remove(bin);
}


TEST_CASE("Duplicate key set() — first wins, save/restore consistent (NGS)",
          "[strategy][bug2][addentry][ngs]") {
    const std::string bin = "bug2_dup_ngs.bin";
    {
        TracEon::SmartStrategy strategy(TracEon::IndexMode::NGS);
        strategy.addEntry("dup", "AAAA", "");
        strategy.addEntry("dup", "CCCC", "");
        REQUIRE(strategy.getFileCacheSize() == 1);
        REQUIRE(strategy.getSequence("dup") == "AAAA");
        strategy.saveBinary(bin);
    }
    {
        TracEon::SmartStrategy restored;
        REQUIRE_NOTHROW(restored.loadBinary(bin));
        REQUIRE(restored.getFileCacheSize() == 1);
        REQUIRE(restored.getSequence("dup") == "AAAA");
    }
    fs::remove(bin);
}
