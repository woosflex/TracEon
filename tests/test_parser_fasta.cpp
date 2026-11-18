#include "TestHelpers.h"

// test_parser_fasta.cpp — FASTA parser domain (normalize,
// header-only, abutting header, terminator, CRLF, MT).
// Split from the monolithic SmartStrategyTests.cpp.


TEST_CASE("Multi-line FASTA normalization", "[strategy][fasta]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "multiline.fasta";
    {
        std::ofstream out(test_file_path);
        // Classic NCBI-style 60-char wrapped FASTA
        out << ">seq1 chromosome 1\n"
               "ACGTACGTACGTACGT\n"
               "TTTTGGGGCCCCAAAA\n"
               ">seq2 chromosome 2\n"
               "GGGGCCCC\n";
        out.close();
    }
    strategy.loadFile(test_file_path);

    REQUIRE(strategy.getFileCacheSize() == 2);
    // Sequence must be contiguous — no embedded newlines
    REQUIRE(strategy.getSequence("seq1") == "ACGTACGTACGTACGTTTTTGGGGCCCCAAAA");
    REQUIRE(strategy.getSequence("seq2") == "GGGGCCCC");

    fs::remove(test_file_path);
}


// ── CRLF (Windows-style) FASTA ─────────────────────────────────────────────────
// Verifies that normalizeFastaArena() and the downstream parsers correctly handle
// \r\n line endings that appear in files created on Windows.

TEST_CASE("FASTA with CRLF line endings", "[strategy][fasta][crlf]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "crlf_test.fasta";
    {
        std::ofstream out(test_file_path, std::ios::binary);
        out << ">seq1 chromosome 1\r\n"
               "ACGTACGTACGTACGT\r\n"
               "TTTTGGGGCCCCAAAA\r\n"
               ">seq2 chromosome 2\r\n"
               "GGGGCCCC\r\n";
        out.close();
    }
    strategy.loadFile(test_file_path);

    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("seq1") == "ACGTACGTACGTACGTTTTTGGGGCCCCAAAA");
    REQUIRE(strategy.getSequence("seq2") == "GGGGCCCC");
    REQUIRE(strategy.getQuality("seq1").empty());
    REQUIRE(strategy.getQuality("seq2").empty());

    fs::remove(test_file_path);
}


TEST_CASE("FASTA with mixed CRLF and LF line endings", "[strategy][fasta][crlf]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "mixed_eol_test.fasta";
    {
        std::ofstream out(test_file_path, std::ios::binary);
        out << ">seq1 first record\r\n"
               "ACGTACGT\r\n"
               ">seq2 second record\n"
               "GGGGCCCC\n";
        out.close();
    }
    strategy.loadFile(test_file_path);

    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("seq1") == "ACGTACGT");
    REQUIRE(strategy.getSequence("seq2") == "GGGGCCCC");

    fs::remove(test_file_path);
}


TEST_CASE("FASTA single record with CRLF", "[strategy][fasta][crlf]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "single_crlf.fasta";
    {
        std::ofstream out(test_file_path, std::ios::binary);
        out << ">chr1\r\nACGTACGTACGTACGTTTTT\r\n";
        out.close();
    }
    strategy.loadFile(test_file_path);

    REQUIRE(strategy.getFileCacheSize() == 1);
    REQUIRE(strategy.getSequence("chr1") == "ACGTACGTACGTACGTTTTT");

    fs::remove(test_file_path);
}


// ── Last record without trailing newline ──────────────────────────────────────
// Files produced by some tools or partial downloads may omit the final \n.
// The parser must still capture the last sequence correctly.

TEST_CASE("FASTA last record without trailing newline", "[strategy][fasta][notrailingnl]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "no_trail_nl.fasta";
    {
        std::ofstream out(test_file_path, std::ios::binary);
        out << ">seq1\nACGTACGTACGTACGT\n"
               ">seq2\nGGGGCCCC";  // <-- no trailing newline
        out.close();
    }
    strategy.loadFile(test_file_path);

    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("seq1") == "ACGTACGTACGTACGT");
    REQUIRE(strategy.getSequence("seq2") == "GGGGCCCC");

    fs::remove(test_file_path);
}


TEST_CASE("FASTA single record without trailing newline", "[strategy][fasta][notrailingnl]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "single_no_trail.fasta";
    {
        std::ofstream out(test_file_path, std::ios::binary);
        out << ">chr1\nACGTACGTACGT";  // no trailing newline at all
        out.close();
    }
    strategy.loadFile(test_file_path);

    REQUIRE(strategy.getFileCacheSize() == 1);
    REQUIRE(strategy.getSequence("chr1") == "ACGTACGTACGT");

    fs::remove(test_file_path);
}


TEST_CASE("FASTA single record with CRLF and no trailing newline", "[strategy][fasta][crlf][notrailingnl]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "crlf_no_trail.fasta";
    {
        std::ofstream out(test_file_path, std::ios::binary);
        out << ">chr1\r\nACGTACGTACGT";  // CRLF header, no trailing newline on sequence
        out.close();
    }
    strategy.loadFile(test_file_path);

    REQUIRE(strategy.getFileCacheSize() == 1);
    REQUIRE(strategy.getSequence("chr1") == "ACGTACGTACGT");

    fs::remove(test_file_path);
}


TEST_CASE("Multithreaded FASTA parsing with SIMD (>10 MB)", "[strategy][multithreaded][fasta][simd]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "mt_fasta_10mb.fasta";

    // Generate >10 MB of FASTA data
    const std::string content = make_large_fasta(11u * 1024u * 1024u); // ~11 MB

    {
        std::ofstream out(test_file_path, std::ios::binary);
        out << content;
    }

    strategy.loadFile(test_file_path);

    // Total record count should be > 0
    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // Verify first, middle, and last records exist by their exact names
    REQUIRE(strategy.hasSequence("record_1"));
    REQUIRE(strategy.hasSequence("record_2"));
    REQUIRE(strategy.hasSequence("record_3"));

    if (total > 3) {
        const std::string mid_key = "record_" + std::to_string(total / 2);
        REQUIRE(strategy.hasSequence(mid_key));
        const std::string last_key = "record_" + std::to_string(total);
        REQUIRE(strategy.hasSequence(last_key));
    }

    // Verify sequence content of first record
    std::string_view seq1 = strategy.getView("record_1");
    REQUIRE_FALSE(seq1.empty());
    // Must be a multiple of 4 (ACGT repeated)
    REQUIRE(seq1.size() % 4 == 0);
    // First 4 chars should be ACGT
    REQUIRE(seq1.substr(0, 4) == "ACGT");
    // Last 4 chars should be ACGT
    REQUIRE(seq1.substr(seq1.size() - 4) == "ACGT");

    // No quality data for FASTA
    REQUIRE(strategy.getQuality("record_1").empty());

    // Verify format detection
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTA);

    fs::remove(test_file_path);
}


TEST_CASE("Multithreaded FASTA with CRLF line endings (>10 MB)", "[strategy][multithreaded][fasta][crlf][simd]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "mt_crlf_fasta.fasta";

    const size_t TARGET = 11u * 1024u * 1024u;
    const size_t SEQ_LEN = 4000;
    const size_t REPEATS = SEQ_LEN / 4;
    const std::string seq = std::string(REPEATS, 'A') + std::string(REPEATS, 'C') +
                            std::string(REPEATS, 'G') + std::string(REPEATS, 'T');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << ">record_" << n << "\r\n";
            out << seq << "\r\n";
            written += 12 + (n > 9 ? (n > 99 ? 4 : 3) : 2) +  // ">record_N\r\n"
                       seq.size() + 2;                           // seq + "\r\n"
            ++n;
        }
    }

    strategy.loadFile(test_file_path);

    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // Verify first record
    REQUIRE(strategy.hasSequence("record_1"));
    std::string_view seq1 = strategy.getView("record_1");
    REQUIRE_FALSE(seq1.empty());
    // seq is built as string(REPEATS,'A') + string(REPEATS,'C') + ...
    REQUIRE(seq1.substr(0, 4) == "AAAA");

    // Verify a middle record and last record
    if (total > 3) {
        const std::string mid_key = "record_" + std::to_string(total / 2);
        REQUIRE(strategy.hasSequence(mid_key));
    }
    const std::string last_key = "record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    REQUIRE_FALSE(strategy.getView(last_key).empty());

    fs::remove(test_file_path);
}


TEST_CASE("Multithreaded FASTA — last record without trailing newline (>10 MB)", "[strategy][multithreaded][fasta][notrailingnl][simd]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "mt_notrail_fasta.fasta";

    const size_t TARGET = 11u * 1024u * 1024u;
    const size_t SEQ_LEN = 4000;
    const size_t REPEATS = SEQ_LEN / 4;
    const std::string seq = std::string(REPEATS, 'A') + std::string(REPEATS, 'C') +
                            std::string(REPEATS, 'G') + std::string(REPEATS, 'T');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            const bool is_last = (written + seq.size() + 20 >= TARGET);
            out << ">record_" << n << "\n";
            if (is_last) {
                out << seq;           // no trailing newline on last record
            } else {
                out << seq << "\n";
            }
            written += 10 + (n > 9 ? (n > 99 ? 4 : 3) : 2) +
                       seq.size() + (is_last ? 0 : 1);
            ++n;
        }
    }

    strategy.loadFile(test_file_path);

    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // Verify first record
    REQUIRE(strategy.hasSequence("record_1"));
    REQUIRE_FALSE(strategy.getView("record_1").empty());

    // Verify last record — the one without trailing newline
    const std::string last_key = "record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    std::string_view last_seq = strategy.getView(last_key);
    REQUIRE_FALSE(last_seq.empty());
    REQUIRE(last_seq.size() == SEQ_LEN);
    // seq is built as string(REPEATS,'A') + string(REPEATS,'C') + ...
    REQUIRE(last_seq.substr(0, 4) == "AAAA");

    fs::remove(test_file_path);
}


// ═══════════════════════════════════════════════════════════════════════════
// Pre-reserved Thread-Local Map Tests
// ═══════════════════════════════════════════════════════════════════════════
//
// These tests verify that the pre-reservation of thread-local maps in
// parseFastaMultithreadedTemplate() and parseFastqMultithreadedTemplate()
// works correctly.  Pre-reservation avoids mid-parse rehashing by estimating
// the record count per thread chunk and calling reserve() upfront.
//
// Tests cover correctness of record capture, no data loss, concurrent read
// safety after MT parsing, and various record-size heuristics.
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Pre-reserved MT FASTA — all records present after threaded parse",
          "[strategy][multithreaded][prereserve][fasta]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "prereserve_mt_fasta.fasta";

    // Generate FASTA large enough to trigger multithreading (>10 MB)
    // Use a mix of short and long sequence sizes to exercise the heuristic
    const size_t TARGET = 11u * 1024u * 1024u;
    const std::string short_seq = std::string(200, 'A') + std::string(200, 'C');
    const std::string long_seq = std::string(8000, 'G') + std::string(8000, 'T');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            const std::string& seq = (n % 5 == 0) ? long_seq : short_seq;
            const std::string line = ">record_" + std::to_string(n) + "\n" + seq + "\n";
            out << line;
            written += line.size();
            ++n;
        }
    }

    strategy.loadFile(test_file_path);
    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // Verify first, middle, and last records by getView
    REQUIRE_FALSE(strategy.getView("record_1").empty());
    REQUIRE(strategy.getSequence("record_1") == short_seq);

    if (total > 3) {
        const std::string mid_key = "record_" + std::to_string(total / 2);
        REQUIRE(strategy.hasSequence(mid_key));
        REQUIRE_FALSE(strategy.getView(mid_key).empty());
    }

    const std::string last_key = "record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    REQUIRE_FALSE(strategy.getView(last_key).empty());

    // Long records (n % 5 == 0) should have size 16000
    REQUIRE(strategy.getView("record_5").size() == 16000);

    // Short records should have size 400
    REQUIRE(strategy.getView("record_1").size() == 400);
    REQUIRE(strategy.getView("record_2").size() == 400);

    // Verify format detection
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTA);

    fs::remove(test_file_path);
}


TEST_CASE("Pre-reserved MT FASTA — concurrent read safety after threaded parse",
          "[strategy][multithreaded][prereserve][concurrent]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "prereserve_concurrent.fasta";

    const size_t TARGET = 11u * 1024u * 1024u;
    const std::string seq(4000, 'A');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << ">record_" << n << "\n" << seq << "\n";
            written += 12 + (n > 9 ? (n > 99 ? 4 : 3) : 2) + seq.size() + 1;
            ++n;
        }
    }

    strategy.loadFile(test_file_path);
    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // Launch concurrent readers while main thread also reads
    // This validates the lock-free data_ready_ signal works after
    // multithreaded parsing with pre-reserved thread-local maps.
    constexpr int NUM_READERS = 8;
    constexpr int LOOKUPS_PER_READER = 5000;
    std::vector<std::thread> readers;
    std::atomic<size_t> total_hits{0};
    std::atomic<bool> stop{false};

    auto reader_worker = [&]() {
        size_t hits = 0;
        for (int i = 0; i < LOOKUPS_PER_READER && !stop.load(); ++i) {
            // Vary the lookup key across first, middle, and last records
            int idx = 0;
            if (i % 3 == 0) idx = 1;
            else if (i % 3 == 1) idx = static_cast<int>(total / 2);
            else idx = static_cast<int>(total);
            std::string key = "record_" + std::to_string(idx);
            auto view = strategy.getView(key);
            if (!view.empty() && view.size() == 4000) ++hits;
        }
        total_hits += hits;
    };

    for (int i = 0; i < NUM_READERS; ++i) readers.emplace_back(reader_worker);

    // Main thread also reads while readers are active
    size_t main_hits = 0;
    for (int i = 0; i < LOOKUPS_PER_READER; ++i) {
        auto view = strategy.getView("record_1");
        if (!view.empty() && view.size() == 4000) ++main_hits;
    }

    for (auto& t : readers) t.join();
    REQUIRE(main_hits == LOOKUPS_PER_READER);
    REQUIRE(total_hits.load() == NUM_READERS * LOOKUPS_PER_READER);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTA);

    fs::remove(test_file_path);
}


TEST_CASE("Pre-reserved MT FASTA — getAllKeys enumerates all records",
          "[strategy][multithreaded][prereserve][keys]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "prereserve_keys.fasta";

    const size_t TARGET = 11u * 1024u * 1024u;
    const std::string seq(500, 'G');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << ">record_" << n << "\n" << seq << "\n";
            written += 12 + (n > 9 ? (n > 99 ? 4 : 3) : 2) + seq.size() + 1;
            ++n;
        }
    }

    strategy.loadFile(test_file_path);
    const size_t total = strategy.getFileCacheSize();

    // getAllKeys must return ALL record names (count matches)
    auto keys = strategy.getAllKeys();
    REQUIRE(keys.size() == total);

    // Verify first and last records are findable by direct lookup
    REQUIRE(strategy.hasSequence("record_1"));
    REQUIRE(strategy.getView("record_1").size() == 500);

    const std::string last_key = "record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    REQUIRE(strategy.getView(last_key).size() == 500);

    // Verify a sample of keys from getAllKeys can be looked up
    for (size_t i = 0; i < std::min(keys.size(), size_t(10)); ++i) {
        REQUIRE(strategy.hasSequence(keys[i]));
        REQUIRE(strategy.getView(keys[i]).size() == 500);
    }

    fs::remove(test_file_path);
}


TEST_CASE("Pre-reserved MT FASTA — many small records (heuristic stress)",
          "[strategy][multithreaded][prereserve][stress]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "prereserve_small_records.fasta";

    // Many very small records: the thread-local reserve heuristic uses
    // chunk_size / 100 which over-estimates for tiny records.
    // This test verifies no data loss even when the heuristic over-reserves.
    const size_t TARGET = 11u * 1024u * 1024u;
    const std::string seq(50, 'T');  // tiny sequences

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << ">r" << n << "\n" << seq << "\n";
            written += 4 + (n > 9 ? (n > 99 ? 3 : 2) : 1) + seq.size() + 1;
            ++n;
        }
    }

    strategy.loadFile(test_file_path);
    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // All records should be accessible
    REQUIRE(strategy.hasSequence("r1"));
    REQUIRE(strategy.getView("r1").size() == 50);

    if (total > 1) {
        const std::string last_key = "r" + std::to_string(total);
        REQUIRE(strategy.hasSequence(last_key));
        REQUIRE(strategy.getView(last_key).size() == 50);
    }

    // getAllKeys should report correct count
    auto keys = strategy.getAllKeys();
    REQUIRE(keys.size() == total);

    fs::remove(test_file_path);
}


TEST_CASE("Pre-reserved MT FASTA — very large records (heuristic under-estimate)",
          "[strategy][multithreaded][prereserve][stress]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "prereserve_large_records.fasta";

    // Few very large records: the heuristic chunk_size / 100 will
    // under-estimate since each record is ~1 MB. Reserve() sets a minimum
    // capacity so even an under-estimate is safe (just triggers rehash).
    const size_t TARGET = 12u * 1024u * 1024u;
    const std::string seq(500000, 'A');  // 500 KB sequences

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << ">large_record_" << n << "\n" << seq << "\n";
            written += 16 + (n > 9 ? (n > 99 ? 3 : 2) : 1) + seq.size() + 1;
            ++n;
        }
    }

    strategy.loadFile(test_file_path);
    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // Verify first and last records
    REQUIRE(strategy.hasSequence("large_record_1"));
    REQUIRE(strategy.getView("large_record_1").size() == 500000);
    REQUIRE(strategy.getView("large_record_1").substr(0, 4) == "AAAA");

    const std::string last_key = "large_record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    REQUIRE(strategy.getView(last_key).size() == 500000);

    // Verify at least a few records
    REQUIRE(strategy.getFileCacheSize() >= 10);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTA);

    fs::remove(test_file_path);
}


TEST_CASE("Pre-reserved MT FASTA — CRLF line endings with pre-reserved maps",
          "[strategy][multithreaded][prereserve][crlf]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "prereserve_crlf.fasta";

    const size_t TARGET = 11u * 1024u * 1024u;
    const std::string seq(4000, 'A');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << ">record_" << n << "\r\n" << seq << "\r\n";
            written += 13 + (n > 9 ? (n > 99 ? 4 : 3) : 2) + seq.size() + 2;
            ++n;
        }
    }

    strategy.loadFile(test_file_path);
    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // Records accessible
    REQUIRE(strategy.hasSequence("record_1"));
    REQUIRE(strategy.getView("record_1").size() == 4000);
    REQUIRE(strategy.getQuality("record_1").empty());

    // Last record
    const std::string last_key = "record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    REQUIRE(strategy.getView(last_key).size() == 4000);

    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTA);
    fs::remove(test_file_path);
}


TEST_CASE("Pre-reserved MT FASTA — last record without trailing newline",
          "[strategy][multithreaded][prereserve][notrailingnl]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "prereserve_notrail.fasta";

    const size_t TARGET = 11u * 1024u * 1024u;
    const std::string seq(4000, 'A');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            const bool is_last = (written + seq.size() + 20 >= TARGET);
            out << ">record_" << n << "\n" << seq;
            if (!is_last) out << "\n";
            written += 12 + (n > 9 ? (n > 99 ? 4 : 3) : 2) +
                       seq.size() + (is_last ? 0 : 1);
            ++n;
        }
    }

    strategy.loadFile(test_file_path);
    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    REQUIRE(strategy.hasSequence("record_1"));
    REQUIRE(strategy.getView("record_1").size() == 4000);

    const std::string last_key = "record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    REQUIRE(strategy.getView(last_key).size() == 4000);

    fs::remove(test_file_path);
}


// ═══════════════════════════════════════════════════════════════════════════
// Bug A regression: validation must reject malformed boundaries before
// normalizeFastaArena can truncate or reinterpret the file.
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("FASTA — a header marker must begin a line",
          "[strategy][fasta][integrity]") {
    TracEon::SmartStrategy strategy;
    const std::string path = "bugA_abut.fasta";
    { std::ofstream out(path, std::ios::binary); out << ">a\nACGT>b\nTTTT\n"; }
    REQUIRE_THROWS_AS(strategy.loadFile(path), std::runtime_error);
    REQUIRE(strategy.getFileCacheSize() == 0);
    REQUIRE(strategy.getArenaBytes() == 0);
    fs::remove(path);
}


TEST_CASE("FASTA — header-only record rejects the entire input",
          "[strategy][fasta][integrity]") {
    TracEon::SmartStrategy strategy;
    std::string path = "bugA_headonly.fasta";
    {
        std::ofstream out(path, std::ios::binary);
        out << ">a\n>b\nACGT\n";
    }
    REQUIRE_THROWS_AS(strategy.loadFile(path), std::runtime_error);
    REQUIRE(strategy.getFileCacheSize() == 0);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::UNKNOWN);
    fs::remove(path);
}


TEST_CASE("FASTA — header-only first record rejects later prefix",
          "[strategy][fasta][integrity]") {
    TracEon::SmartStrategy strategy;
    std::string path = "bugA_headonly_start.fasta";
    {
        std::ofstream out(path, std::ios::binary);
        out << ">empty\n>full\nACGT\n>more\nTTTT\n";
    }
    REQUIRE_THROWS_AS(strategy.loadFile(path), std::runtime_error);
    REQUIRE(strategy.getFileCacheSize() == 0);
    fs::remove(path);
}


TEST_CASE("FASTA — multi-line sequence normalization (control)",
          "[strategy][fasta][bugA]") {
    // Control: ordinary multi-line FASTA must normalize to a single line.
    TracEon::SmartStrategy strategy;
    std::string path = "bugA_multiline.fasta";
    {
        std::ofstream out(path, std::ios::binary);
        out << ">a\nACGT\nTTTT\n>b\nGG\n";
    }
    strategy.loadFile(path);
    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("a") == "ACGTTTTT"); // "ACGT" + "TTTT"
    REQUIRE(strategy.getSequence("b") == "GG");
    fs::remove(path);
}


TEST_CASE("FASTA — EOF without trailing newline (control)",
          "[strategy][fasta][bugA]") {
    TracEon::SmartStrategy strategy;
    std::string path = "bugA_eof.fasta";
    {
        std::ofstream out(path, std::ios::binary);
        // Last record has no trailing newline: the +1 resize slot must hold
        // the generated terminator.
        out << ">a\nACGT";
    }
    strategy.loadFile(path);
    REQUIRE(strategy.getFileCacheSize() == 1);
    REQUIRE(strategy.getSequence("a") == "ACGT");
    fs::remove(path);
}


TEST_CASE("FASTA — CRLF records and header-only rejection",
          "[strategy][fasta][bugA][crlf]") {
    TracEon::SmartStrategy strategy;
    std::string path = "bugA_crlf.fasta";
    {
        std::ofstream out(path, std::ios::binary);
        out << ">a\r\nACGT\r\n>b\r\nTTTT\r\n>c\r\nGGGG\r\n";
    }
    strategy.loadFile(path);
    REQUIRE(strategy.getFileCacheSize() == 3);
    REQUIRE(strategy.getSequence("a") == "ACGT");
    REQUIRE(strategy.getSequence("b") == "TTTT");
    REQUIRE(strategy.getSequence("c") == "GGGG");
    fs::remove(path);

    // CRLF header-only record in the middle.
    std::string path2 = "bugA_crlf_headonly.fasta";
    {
        std::ofstream out2(path2, std::ios::binary);
        out2 << ">x\r\n>y\r\nACGT\r\n";
    }
    REQUIRE_THROWS_AS(strategy.loadFile(path2), std::runtime_error);
    REQUIRE(strategy.getFileCacheSize() == 0);
    fs::remove(path2);
}


TEST_CASE("FASTA — header-only record at EOF rejects valid prefix",
          "[strategy][fasta][integrity]") {
    TracEon::SmartStrategy strategy;
    std::string path = "bugA_headonly_eof.fasta";
    {
        std::ofstream out(path, std::ios::binary);
        out << ">a\nACGT\n>b\n"; // 'b' is header-only at EOF
    }
    REQUIRE_THROWS_AS(strategy.loadFile(path), std::runtime_error);
    REQUIRE(strategy.getFileCacheSize() == 0);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::UNKNOWN);
    fs::remove(path);
}

TEST_CASE("Strict FASTA header grammar is identical for serial/MT and plain/GZIP",
          "[strategy][fasta][integrity][grammar]") {
    const std::string seq(1000, 'A');
    auto expand = [&](const std::string& first) {
        std::string data = first;
        for (size_t i = 0; data.size() <= 11u * 1024u * 1024u; ++i)
            data += ">r" + std::to_string(i) + "\n" + seq + "\n";
        return data;
    };
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"leading-space", " \t>bad\nACGT\n"},
        {"header-space", "> bad\nACGT\n"},
    };
    for (const auto& [kind, small] : cases) {
        for (const auto& [size_name, content] :
             std::vector<std::pair<std::string, std::string>>{{"serial", small}, {"mt", expand(small)}}) {
            for (bool gzip : {false, true}) {
                INFO(kind << ' ' << size_name << (gzip ? " gzip" : " plain"));
                const std::string path = "strict_fasta_" + kind + '_' + size_name +
                    (gzip ? ".fa.gz" : ".fa");
                if (gzip) (void)compress_to_gzip_bytes(content, path);
                else { std::ofstream out(path, std::ios::binary); out << content; }
                TracEon::SmartStrategy strategy;
                REQUIRE_THROWS_AS(strategy.loadFile(path), std::runtime_error);
                REQUIRE(strategy.getFileCacheSize() == 0);
                REQUIRE(strategy.getArenaBytes() == 0);
                REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::UNKNOWN);
                const std::string valid = "strict_fasta_reuse.fa";
                { std::ofstream out(valid, std::ios::binary); out << ">ok\nACGT\n"; }
                REQUIRE_NOTHROW(strategy.loadFile(valid));
                REQUIRE(strategy.hasSequence("ok"));
                fs::remove(valid);
                fs::remove(path);
            }
        }
    }
}

TEST_CASE("Trailing FASTA header rejects the entire plain or GZIP input",
          "[strategy][fasta][integrity]") {
    const std::string small = ">ok\nACGT\n>bad\n";
    std::string large;
    large.reserve(12u * 1024u * 1024u);
    const std::string seq(1000, 'A');
    for (size_t i = 0; large.size() <= 11u * 1024u * 1024u; ++i)
        large += ">r" + std::to_string(i) + "\n" + seq + "\n";
    large += ">bad\n";

    for (const auto& [label, content] :
         std::vector<std::pair<std::string, std::string>>{{"serial", small}, {"multithreaded", large}}) {
        for (bool gzip : {false, true}) {
            INFO(label << (gzip ? " gzip" : " plain"));
            const std::string path = "trailing_header_" + label + (gzip ? ".fasta.gz" : ".fasta");
            if (gzip)
                (void)compress_to_gzip_bytes(content, path);
            else {
                std::ofstream out(path, std::ios::binary);
                out << content;
            }
            TracEon::SmartStrategy strategy;
            REQUIRE_THROWS_AS(strategy.loadFile(path), std::runtime_error);
            REQUIRE(strategy.getFileCacheSize() == 0);
            REQUIRE_FALSE(strategy.hasSequence("ok"));
            REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::UNKNOWN);
            fs::remove(path);
        }
    }
}

TEST_CASE("FASTA rejects empty CR-delimited identifiers and preserves payload CR bytes", "[audit][fasta]") {
    const std::string path = "audit_cr.fasta";
    {
        std::ofstream out(path, std::ios::binary);
        out << ">\rdescription\nACGT\n";
    }
    TracEon::SmartStrategy cache;
    REQUIRE_THROWS(cache.loadFile(path));
    REQUIRE(cache.getFileCacheSize() == 0);
    {
        std::ofstream out(path, std::ios::binary);
        out << ">r\r\n\rACGT\r\r\n";
    }
    cache.loadFile(path);
    REQUIRE(cache.getView("r") == "\rACGT\r");
    fs::remove(path);
}
