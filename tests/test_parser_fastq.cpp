#include "TestHelpers.h"

// test_parser_fastq.cpp — FASTQ parser domain (framing, empty
// quality/seq, MT chunk classifier, bug1/bugB repro cases).
// Split from the monolithic SmartStrategyTests.cpp.


TEST_CASE("Multithreaded FASTQ parsing with SIMD (>10 MB)", "[strategy][multithreaded][fastq][simd]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "mt_fastq_10mb.fastq";

    const std::string content = make_large_fastq(11u * 1024u * 1024u); // ~11 MB

    {
        std::ofstream out(test_file_path, std::ios::binary);
        out << content;
    }

    strategy.loadFile(test_file_path);

    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // Verify first, middle, and last records
    REQUIRE(strategy.hasSequence("record_1"));
    REQUIRE(strategy.hasSequence("record_2"));
    REQUIRE(strategy.hasSequence("record_3"));

    if (total > 3) {
        const std::string mid_key = "record_" + std::to_string(total / 2);
        REQUIRE(strategy.hasSequence(mid_key));
        const std::string last_key = "record_" + std::to_string(total);
        REQUIRE(strategy.hasSequence(last_key));
    }

    // Verify sequence + quality for first record
    std::string_view seq1 = strategy.getView("record_1");
    REQUIRE_FALSE(seq1.empty());
    REQUIRE(seq1.size() % 4 == 0);
    REQUIRE(strategy.hasSequence("record_1"));

    // Quality should be present for FASTQ
    std::string qual1 = strategy.getQuality("record_1");
    REQUIRE_FALSE(qual1.empty());
    REQUIRE(qual1.size() == seq1.size());
    // Quality chars should be 'I'
    REQUIRE(qual1[0] == 'I');

    // Verify format detection
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTQ);

    fs::remove(test_file_path);
}


TEST_CASE("Multithreaded FASTQ with CRLF line endings (>10 MB)", "[strategy][multithreaded][fastq][crlf][simd]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "mt_crlf_fastq.fastq";

    const size_t TARGET = 11u * 1024u * 1024u;
    const size_t SEQ_LEN = 4000;
    const size_t REPEATS = SEQ_LEN / 4;
    const std::string seq = std::string(REPEATS, 'A') + std::string(REPEATS, 'C') +
                            std::string(REPEATS, 'G') + std::string(REPEATS, 'T');
    const std::string qual(SEQ_LEN, 'I');
    const std::string plus = "+";

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << "@record_" << n << "\r\n";
            out << seq << "\r\n";
            out << "+\r\n";
            out << qual << "\r\n";
            written += 11 + (n > 9 ? (n > 99 ? 4 : 3) : 2) +
                       SEQ_LEN + 2 +  // seq + "\r\n"
                       3 +            // "+\r\n"
                       SEQ_LEN + 2;   // qual + "\r\n"
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
    REQUIRE(seq1.size() == SEQ_LEN);

    // Quality string should match
    const std::string q1 = strategy.getQuality("record_1");
    REQUIRE_FALSE(q1.empty());
    REQUIRE(q1.size() == SEQ_LEN);
    REQUIRE(q1[0] == 'I');

    // Verify last record
    const std::string last_key = "record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    REQUIRE_FALSE(strategy.getView(last_key).empty());

    fs::remove(test_file_path);
}


TEST_CASE("Multithreaded FASTQ — last record without trailing newline (>10 MB)", "[strategy][multithreaded][fastq][notrailingnl][simd]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "mt_notrail_fastq.fastq";

    const size_t TARGET = 11u * 1024u * 1024u;
    const size_t SEQ_LEN = 4000;
    const size_t REPEATS = SEQ_LEN / 4;
    const std::string seq = std::string(REPEATS, 'A') + std::string(REPEATS, 'C') +
                            std::string(REPEATS, 'G') + std::string(REPEATS, 'T');
    const std::string qual(SEQ_LEN, 'I');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << "@record_" << n << "\n";
            out << seq << "\n";
            out << "+\n";
            const bool is_last = (written + 30 + seq.size() + qual.size() >= TARGET);
            if (is_last) {
                out << qual;          // no trailing newline after quality on last record
            } else {
                out << qual << "\n";
            }
            written += 11 + (n > 9 ? (n > 99 ? 4 : 3) : 2) +
                       SEQ_LEN + 1 +   // seq + \n
                       2 +             // +\n
                       SEQ_LEN + (is_last ? 0 : 1);
            ++n;
        }
    }

    strategy.loadFile(test_file_path);

    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // Verify first record
    REQUIRE(strategy.hasSequence("record_1"));
    REQUIRE_FALSE(strategy.getView("record_1").empty());
    REQUIRE(strategy.getQuality("record_1").size() == SEQ_LEN);

    // Verify last record (no trailing newline)
    const std::string last_key = "record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    std::string_view last_seq = strategy.getView(last_key);
    REQUIRE_FALSE(last_seq.empty());
    REQUIRE(last_seq.size() == SEQ_LEN);

    std::string last_qual = strategy.getQuality(last_key);
    REQUIRE_FALSE(last_qual.empty());
    REQUIRE(last_qual.size() == SEQ_LEN);
    REQUIRE(last_qual[0] == 'I');

    fs::remove(test_file_path);
}


// ═══════════════════════════════════════════════════════════════════════════
// Bug 3 regression: '+'-leading quality lines must not drop MT chunks
// ═══════════════════════════════════════════════════════════════════════════
// Phred+33 quality Q10 is ASCII '+' (0x2B), so real FASTQ files routinely
// have quality lines that START with '+'. The multithreaded chunk-boundary
// scan used to classify a candidate '@' header by walking BACKWARD to the
// previous line and testing '+'; a '+'-leading quality line is the previous
// line of the first header at/after a chunk boundary, so that header was
// misclassified as a quality line and silently dropped — every worker after
// chunk 0 lost its entire chunk (live repro on 20483dc: 80k reads → 4,468
// loaded, 94.4% loss). The classifier now looks FORWARD: a genuine header
// has seq at line+1 and '+' at line+2; a quality '@' has the next record's
// header at line+1 and its seq at line+2, which never starts with '+'.

TEST_CASE("FASTQ '+'-leading quality — multithreaded (>10MB) record count preserved",
          "[strategy][fastq][bug-plus-quality][multithreaded]") {
    // Mirrors the live repro: 80k × 100 bp, every quality line starts '+'
    // (Q10). 80k × ~205 B ≈ 16.4 MB → multithreaded path.
    const size_t SEQ_LEN = 100;
    const size_t N = 80000;
    std::string path = "bug3_mt_plus_quality.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        for (size_t i = 0; i < N; ++i) {
            out << "@read" << i << "\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << "ACGT"[j % 4];
            out << "\n+\n";
            out << '+'; // quality starts with '+' (Phred+33 Q10)
            for (size_t j = 1; j < SEQ_LEN; ++j) out << 'I';
            out << "\n";
        }
    }
    REQUIRE(fs::file_size(path) > 10u * 1024u * 1024u); // must take the MT path

    TracEon::SmartStrategy strategy;
    strategy.loadFile(path);

    // Every record must be present — no loss, no garbage.
    REQUIRE(strategy.getFileCacheSize() == N);
    REQUIRE(strategy.hasSequence("read0"));
    REQUIRE(strategy.getSequence("read0").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("read0").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("read0")[0] == '+');
    REQUIRE(strategy.getQuality("read0")[1] == 'I');
    REQUIRE(strategy.hasSequence("read19999"));
    REQUIRE(strategy.hasSequence("read39999"));
    REQUIRE(strategy.hasSequence("read59999"));
    REQUIRE(strategy.hasSequence("read79999"));
    REQUIRE(strategy.getSequence("read79999").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("read79999")[0] == '+');

    fs::remove(path);
}


TEST_CASE("FASTQ '+'-leading quality with CRLF — multithreaded (>10MB)",
          "[strategy][fastq][bug-plus-quality][crlf][multithreaded]") {
    // '+'-leading quality + CRLF line endings: the forward classifier's
    // two-line lookahead must skip '\r' correctly at chunk boundaries.
    const size_t SEQ_LEN = 300;
    const size_t N = 40000; // 40000 × ~605 B ≈ 24 MB → multithreaded path
    std::string path = "bug3_mt_plus_quality_crlf.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        for (size_t i = 0; i < N; ++i) {
            out << "@read" << i << "\r\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << "ACGT"[j % 4];
            out << "\r\n+\r\n";
            out << '+';
            for (size_t j = 1; j < SEQ_LEN; ++j) out << 'I';
            out << "\r\n";
        }
    }
    REQUIRE(fs::file_size(path) > 10u * 1024u * 1024u);

    TracEon::SmartStrategy strategy;
    strategy.loadFile(path);

    REQUIRE(strategy.getFileCacheSize() == N);
    REQUIRE(strategy.getSequence("read0").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("read0")[0] == '+');
    REQUIRE(strategy.hasSequence("read19999"));
    REQUIRE(strategy.hasSequence("read39999"));
    REQUIRE(strategy.getQuality("read39999")[0] == '+');

    fs::remove(path);
}


TEST_CASE("FASTQ with CRLF line endings — single-threaded (<10MB)",
          "[strategy][fastq][crlf]") {
    // ST path has no dedicated CRLF FASTQ coverage until now.
    const size_t SEQ_LEN = 100;
    const size_t N = 5000; // ~1 MB → single-threaded path
    std::string path = "st_crlf_fastq.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        for (size_t i = 0; i < N; ++i) {
            out << "@read" << i << " desc\r\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << "ACGT"[j % 4];
            out << "\r\n+\r\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << 'I';
            out << "\r\n";
        }
    }
    REQUIRE(fs::file_size(path) < 10u * 1024u * 1024u); // ST path

    TracEon::SmartStrategy strategy;
    strategy.loadFile(path);

    REQUIRE(strategy.getFileCacheSize() == N);
    REQUIRE(strategy.getSequence("read0").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("read0") == std::string(SEQ_LEN, 'I'));
    REQUIRE(strategy.hasSequence("read4999"));
    REQUIRE(strategy.getSequence("read4999").size() == SEQ_LEN);

    fs::remove(path);
}


TEST_CASE("FASTQ without trailing newline — single-threaded (<10MB)",
          "[strategy][fastq][notrailingnl]") {
    // ST path: last record's quality line has no trailing newline.
    const size_t SEQ_LEN = 100;
    const size_t N = 5000; // ~1 MB → single-threaded path
    std::string path = "st_notrail_fastq.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        for (size_t i = 0; i < N - 1; ++i) {
            out << "@read" << i << "\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << "ACGT"[j % 4];
            out << "\n+\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << 'I';
            out << "\n";
        }
        // Last record: no trailing newline after the quality line.
        out << "@read" << (N - 1) << "\n";
        for (size_t j = 0; j < SEQ_LEN; ++j) out << "ACGT"[j % 4];
        out << "\n+\n";
        for (size_t j = 0; j < SEQ_LEN; ++j) out << 'I';
    }

    TracEon::SmartStrategy strategy;
    strategy.loadFile(path);

    REQUIRE(strategy.getFileCacheSize() == N);
    REQUIRE(strategy.getSequence("read0").size() == SEQ_LEN);
    const std::string last_key = "read" + std::to_string(N - 1);
    REQUIRE(strategy.hasSequence(last_key));
    REQUIRE(strategy.getSequence(last_key).size() == SEQ_LEN);
    REQUIRE(strategy.getQuality(last_key) == std::string(SEQ_LEN, 'I'));

    fs::remove(path);
}


TEST_CASE("Pre-reserved MT FASTQ — all records present after threaded parse",
          "[strategy][multithreaded][prereserve][fastq]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "prereserve_mt_fastq.fastq";

    const size_t TARGET = 11u * 1024u * 1024u;
    const size_t SEQ_LEN = 4000;
    const size_t REPEATS = SEQ_LEN / 4;
    const std::string seq = std::string(REPEATS, 'A') + std::string(REPEATS, 'C') +
                            std::string(REPEATS, 'G') + std::string(REPEATS, 'T');
    const std::string qual(SEQ_LEN, 'I');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << "@record_" << n << "\n";
            out << seq << "\n";
            out << "+\n";
            out << qual << "\n";
            written += 11 + (n > 9 ? (n > 99 ? 4 : 3) : 2) +
                       SEQ_LEN + 1 + 2 + SEQ_LEN + 1;
            ++n;
        }
    }

    strategy.loadFile(test_file_path);
    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // First record
    REQUIRE(strategy.hasSequence("record_1"));
    std::string_view seq1 = strategy.getView("record_1");
    REQUIRE_FALSE(seq1.empty());
    REQUIRE(seq1.size() == SEQ_LEN);
    REQUIRE(seq1.substr(0, 4) == "AAAA");

    // Quality
    std::string q1 = strategy.getQuality("record_1");
    REQUIRE_FALSE(q1.empty());
    REQUIRE(q1.size() == SEQ_LEN);
    REQUIRE(q1[0] == 'I');

    // Middle record
    if (total > 3) {
        const std::string mid_key = "record_" + std::to_string(total / 2);
        REQUIRE(strategy.hasSequence(mid_key));
        REQUIRE_FALSE(strategy.getView(mid_key).empty());
        REQUIRE(strategy.getQuality(mid_key).size() == SEQ_LEN);
    }

    // Last record
    const std::string last_key = "record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    REQUIRE_FALSE(strategy.getView(last_key).empty());
    REQUIRE(strategy.getQuality(last_key).size() == SEQ_LEN);

    // Format detection
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTQ);

    fs::remove(test_file_path);
}


TEST_CASE("Pre-reserved MT FASTQ — concurrent read safety after threaded parse",
          "[strategy][multithreaded][prereserve][concurrent][fastq]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "prereserve_concurrent_fastq.fastq";

    const size_t TARGET = 11u * 1024u * 1024u;
    const size_t SEQ_LEN = 400;
    const std::string seq = std::string(SEQ_LEN / 2, 'A') + std::string(SEQ_LEN / 2, 'C');
    const std::string qual(SEQ_LEN, 'J');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << "@record_" << n << "\n" << seq << "\n+\n" << qual << "\n";
            written += 11 + (n > 9 ? (n > 99 ? 4 : 3) : 2) +
                       SEQ_LEN + 1 + 2 + SEQ_LEN + 1;
            ++n;
        }
    }

    strategy.loadFile(test_file_path);
    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // Verify quality data accessible after MT parse
    REQUIRE(strategy.getQuality("record_1") == qual);
    REQUIRE(strategy.getQuality("record_1").size() == SEQ_LEN);

    // Concurrent readers
    constexpr int NUM_READERS = 4;
    constexpr int LOOKUPS = 2000;
    std::vector<std::thread> readers;
    std::atomic<size_t> total_seq_hits{0};
    std::atomic<size_t> total_qual_hits{0};

    for (int r = 0; r < NUM_READERS; ++r) {
        readers.emplace_back([&]() {
            size_t seq_hits = 0, qual_hits = 0;
            for (int i = 0; i < LOOKUPS; ++i) {
                auto view = strategy.getView("record_1");
                if (view.size() == SEQ_LEN) ++seq_hits;
                std::string q = strategy.getQuality("record_1");
                if (q.size() == SEQ_LEN && q[0] == 'J') ++qual_hits;
            }
            total_seq_hits += seq_hits;
            total_qual_hits += qual_hits;
        });
    }
    for (auto& t : readers) t.join();

    REQUIRE(total_seq_hits.load() == NUM_READERS * LOOKUPS);
    REQUIRE(total_qual_hits.load() == NUM_READERS * LOOKUPS);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTQ);

    fs::remove(test_file_path);
}


TEST_CASE("Pre-reserved MT FASTQ — CRLF line endings with pre-reserved maps",
          "[strategy][multithreaded][prereserve][crlf][fastq]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "prereserve_crlf_fastq.fastq";

    const size_t TARGET = 11u * 1024u * 1024u;
    const size_t SEQ_LEN = 400;
    const std::string seq = std::string(SEQ_LEN, 'A');
    const std::string qual(SEQ_LEN, 'I');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << "@record_" << n << "\r\n";
            out << seq << "\r\n";
            out << "+\r\n";
            out << qual << "\r\n";
            written += 12 + (n > 9 ? (n > 99 ? 4 : 3) : 2) +
                       SEQ_LEN + 2 + 3 + SEQ_LEN + 2;
            ++n;
        }
    }

    strategy.loadFile(test_file_path);
    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    // Verify record
    REQUIRE(strategy.hasSequence("record_1"));
    REQUIRE(strategy.getView("record_1").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("record_1").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("record_1")[0] == 'I');

    // Last record
    const std::string last_key = "record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    REQUIRE_FALSE(strategy.getView(last_key).empty());
    REQUIRE(strategy.getQuality(last_key).size() == SEQ_LEN);

    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTQ);
    fs::remove(test_file_path);
}


TEST_CASE("Pre-reserved MT FASTQ — last record without trailing newline",
          "[strategy][multithreaded][prereserve][notrailingnl][fastq]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "prereserve_notrail_fastq.fastq";

    const size_t TARGET = 11u * 1024u * 1024u;
    const size_t SEQ_LEN = 400;
    const std::string seq(SEQ_LEN, 'A');
    const std::string qual(SEQ_LEN, 'I');

    {
        std::ofstream out(test_file_path, std::ios::binary);
        int n = 1;
        size_t written = 0;
        while (written < TARGET) {
            out << "@record_" << n << "\n" << seq << "\n+\n";
            const bool is_last = (written + 20 + seq.size() + qual.size() >= TARGET);
            if (is_last) {
                out << qual;   // no trailing newline after quality
            } else {
                out << qual << "\n";
            }
            written += 11 + (n > 9 ? (n > 99 ? 4 : 3) : 2) +
                       SEQ_LEN + 1 + 2 + SEQ_LEN + (is_last ? 0 : 1);
            ++n;
        }
    }

    strategy.loadFile(test_file_path);
    const size_t total = strategy.getFileCacheSize();
    REQUIRE(total > 0);

    REQUIRE(strategy.hasSequence("record_1"));
    REQUIRE(strategy.getView("record_1").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("record_1").size() == SEQ_LEN);

    const std::string last_key = "record_" + std::to_string(total);
    REQUIRE(strategy.hasSequence(last_key));
    REQUIRE(strategy.getView(last_key).size() == SEQ_LEN);
    REQUIRE(strategy.getQuality(last_key).size() == SEQ_LEN);
    REQUIRE(strategy.getQuality(last_key)[0] == 'I');

    fs::remove(test_file_path);
}


// ═══════════════════════════════════════════════════════════════════════════
// Bug 1 regression: FASTQ parsers must not misparse '@' in quality lines
// ═══════════════════════════════════════════════════════════════════════════
// Phred+33 quality scores ≥ 64 (Q31+) are ASCII '@' or higher, so real FASTQ
// files routinely contain '@' in quality lines. The old multithreaded
// chunk-boundary scan tested *every* '@' character (not just line-start
// headers), and its lookbehind landed on the wrong line for mid-line '@'s,
// so a quality '@' could be mistaken for a new record header — shifting
// records and silently dropping most of the file (repro: 100k reads → 16,703).
// The parsers now consume strict 4-line cycles (header, seq, '+', qual) and
// the boundary scan only ever inspects line-start '@' characters.

TEST_CASE("FASTQ '@' in quality — multithreaded (>10MB) record count preserved",
          "[strategy][fastq][bug1][multithreaded]") {
    const size_t SEQ_LEN = 300;
    const size_t N = 40000; // 40000 × ~305 B ≈ 12.2 MB → multithreaded path
    std::string path = "bug1_mt_at_quality.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        for (size_t i = 0; i < N; ++i) {
            out << "@read" << i << " desc\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << "ACGT"[j % 4];
            out << "\n+\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << '@'; // Phred Q31+ — every qual char is '@'
            out << "\n";
        }
    }
    REQUIRE(fs::file_size(path) > 10u * 1024u * 1024u); // must take the multithreaded path

    TracEon::SmartStrategy strategy;
    strategy.loadFile(path);

    // Every record must be present with intact content — no loss, no garbage.
    REQUIRE(strategy.getFileCacheSize() == N);
    REQUIRE(strategy.getSequence("read0").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("read0") == std::string(SEQ_LEN, '@'));
    REQUIRE(strategy.getSequence("read19999").size() == SEQ_LEN);
    REQUIRE(strategy.hasSequence("read39999"));
    REQUIRE(strategy.getQuality("read39999") == std::string(SEQ_LEN, '@'));

    fs::remove(path);
}


TEST_CASE("FASTQ '@' in quality — multithreaded, '@' mid-line",
          "[strategy][fastq][bug1][multithreaded]") {
    // '@' NOT at the start of the quality line — the worst case for the old
    // character-scanning boundary logic.
    const size_t SEQ_LEN = 300;
    const size_t N = 40000;
    std::string path = "bug1_mt_midline_at.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        for (size_t i = 0; i < N; ++i) {
            out << "@read" << i << "\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << "ACGT"[j % 4];
            out << "\n+\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << ((j % 3 == 0) ? '@' : 'I');
            out << "\n";
        }
    }
    REQUIRE(fs::file_size(path) > 10u * 1024u * 1024u);

    TracEon::SmartStrategy strategy;
    strategy.loadFile(path);
    REQUIRE(strategy.getFileCacheSize() == N);
    REQUIRE(strategy.hasSequence("read39999"));
    REQUIRE(strategy.getQuality("read1").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("read1")[0] == '@');
    REQUIRE(strategy.getQuality("read1")[1] == 'I');

    fs::remove(path);
}


TEST_CASE("FASTQ '@' in quality — single-threaded (<10MB)",
          "[strategy][fastq][bug1]") {
    const size_t SEQ_LEN = 100;
    const size_t N = 5000; // ~500 KB → single-threaded path
    std::string path = "bug1_st_at_quality.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        for (size_t i = 0; i < N; ++i) {
            out << "@read" << i << " desc\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << "ACGT"[j % 4];
            out << "\n+\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << ((j % 2) ? '@' : '~');
            out << "\n";
        }
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(path);
    REQUIRE(strategy.getFileCacheSize() == N);
    REQUIRE(strategy.getSequence("read0").size() == SEQ_LEN);
    std::string expected_qual;
    expected_qual.reserve(SEQ_LEN);
    for (size_t j = 0; j < SEQ_LEN; ++j) expected_qual += (j % 2) ? '@' : '~';
    REQUIRE(strategy.getQuality("read0") == expected_qual);
    REQUIRE(strategy.hasSequence("read4999"));

    fs::remove(path);
}


// ═══════════════════════════════════════════════════════════════════════════
// Bug B regression: FASTQ newline-run collapse breaks the 4-line cycle
// ═══════════════════════════════════════════════════════════════════════════
// The parse loops used to collapse RUNS of newlines after each line. That is
// wrong for positional FASTQ when a line is EMPTY: an empty quality line
// ('@e0\nAAA\n+\n\n@e1\n...') made qual_start land on '@e1' (record e0
// stored with the wrong quality, then the cycle broke and every later record
// was dropped), and an empty sequence line ('@e1\n\n+\n...') made seq_start
// land on '+' (mis-parse, cycle break, whole file lost). Zero-length reads
// are valid FASTQ (fastp/cutadapt with --length-required 0 emit them). Every
// boundary now consumes exactly ONE line terminator ('\n', with the optional
// preceding '\r' already part of the previous line).

TEST_CASE("FASTQ — empty quality line — single-threaded (LF)",
          "[strategy][fastq][bugB]") {
    TracEon::SmartStrategy strategy;
    std::string path = "bugB_st_emptyqual.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        out << "@e0\nAAA\n+\n\n@e1\nTTT\n+\nIII\n";
    }
    strategy.loadFile(path);
    REQUIRE(strategy.getFileCacheSize() == 2); // was: 1 (e1 dropped)
    REQUIRE(strategy.getSequence("e0") == "AAA");
    REQUIRE(strategy.getQuality("e0").empty()); // was: quality = "@e1"
    REQUIRE(strategy.getSequence("e1") == "TTT");
    REQUIRE(strategy.getQuality("e1") == "III");
    fs::remove(path);
}


TEST_CASE("FASTQ — empty quality line — single-threaded (CRLF)",
          "[strategy][fastq][bugB][crlf]") {
    TracEon::SmartStrategy strategy;
    std::string path = "bugB_st_emptyqual_crlf.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        out << "@e0\r\nAAA\r\n+\r\n\r\n@e1\r\nTTT\r\n+\r\nIII\r\n";
    }
    strategy.loadFile(path);
    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("e0") == "AAA");
    REQUIRE(strategy.getQuality("e0").empty());
    REQUIRE(strategy.getSequence("e1") == "TTT");
    REQUIRE(strategy.getQuality("e1") == "III");
    fs::remove(path);
}


TEST_CASE("FASTQ — empty sequence line (zero-length read) — single-threaded",
          "[strategy][fastq][bugB]") {
    // Zero-length read e1 ('@e1\n\n+\n\n') sits between two normal records.
    // Empty-seq records are not stored (existing invariant), but the cycle
    // must continue so e2 survives (was: cycle broke at the '+' check and
    // the whole file was lost).
    TracEon::SmartStrategy strategy;
    std::string path = "bugB_st_emptyseq.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        out << "@e0\nAAA\n+\nIII\n@e1\n\n+\n\n@e2\nTTT\n+\nJJJ\n";
    }
    strategy.loadFile(path);
    REQUIRE(strategy.getFileCacheSize() == 2); // e0 + e2 (e1 dropped)
    REQUIRE(strategy.getSequence("e0") == "AAA");
    REQUIRE(strategy.getQuality("e0") == "III");
    REQUIRE_FALSE(strategy.hasSequence("e1"));
    REQUIRE(strategy.getSequence("e2") == "TTT");
    REQUIRE(strategy.getQuality("e2") == "JJJ");
    fs::remove(path);
}


TEST_CASE("FASTQ — empty quality line — multithreaded (>10MB, LF)",
          "[strategy][fastq][bugB][multithreaded]") {
    const size_t SEQ_LEN = 200;
    const size_t N = 40000; // ~16.6 MB → multithreaded path
    std::string path = "bugB_mt_emptyqual.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        for (size_t i = 0; i < N; ++i) {
            out << "@read" << i << "\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << "ACGT"[j % 4];
            out << "\n+\n";
            if (i % 997 == 0) {
                out << "\n"; // empty quality line (record still stored with qual='')
            } else {
                for (size_t j = 0; j < SEQ_LEN; ++j) out << 'I';
                out << "\n";
            }
        }
    }
    REQUIRE(fs::file_size(path) > 10u * 1024u * 1024u);

    TracEon::SmartStrategy strategy;
    strategy.loadFile(path);

    REQUIRE(strategy.getFileCacheSize() == N); // was: ~4.4k (loss after first empty qual)
    // First record has empty quality (i=0): must be stored with qual=''.
    REQUIRE(strategy.getSequence("read0").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("read0").empty());
    // Record after an empty-quality record must carry its OWN quality.
    REQUIRE(strategy.getSequence("read1").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("read1") == std::string(SEQ_LEN, 'I'));
    REQUIRE(strategy.getQuality("read997").empty());
    REQUIRE(strategy.getQuality("read998") == std::string(SEQ_LEN, 'I'));
    REQUIRE(strategy.getQuality("read1994").empty());
    REQUIRE(strategy.getQuality("read1995") == std::string(SEQ_LEN, 'I'));
    // Last record intact.
    const std::string last = "read" + std::to_string(N - 1);
    REQUIRE(strategy.hasSequence(last));
    REQUIRE(strategy.getSequence(last).size() == SEQ_LEN);
    REQUIRE(strategy.getQuality(last) == std::string(SEQ_LEN, 'I'));

    fs::remove(path);
}


TEST_CASE("FASTQ — empty quality line — multithreaded (>10MB, CRLF)",
          "[strategy][fastq][bugB][multithreaded][crlf]") {
    const size_t SEQ_LEN = 300;
    const size_t N = 40000; // ~25 MB → multithreaded path
    std::string path = "bugB_mt_emptyqual_crlf.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        for (size_t i = 0; i < N; ++i) {
            out << "@read" << i << "\r\n";
            for (size_t j = 0; j < SEQ_LEN; ++j) out << "ACGT"[j % 4];
            out << "\r\n+\r\n";
            if (i % 997 == 0) {
                out << "\r\n"; // empty quality line
            } else {
                for (size_t j = 0; j < SEQ_LEN; ++j) out << 'I';
                out << "\r\n";
            }
        }
    }
    REQUIRE(fs::file_size(path) > 10u * 1024u * 1024u);

    TracEon::SmartStrategy strategy;
    strategy.loadFile(path);

    REQUIRE(strategy.getFileCacheSize() == N);
    REQUIRE(strategy.getQuality("read0").empty());
    REQUIRE(strategy.getQuality("read1") == std::string(SEQ_LEN, 'I'));
    REQUIRE(strategy.getQuality("read997").empty());
    REQUIRE(strategy.getQuality("read998") == std::string(SEQ_LEN, 'I'));
    REQUIRE(strategy.getSequence("read0").size() == SEQ_LEN);
    const std::string last = "read" + std::to_string(N - 1);
    REQUIRE(strategy.hasSequence(last));
    REQUIRE(strategy.getQuality(last) == std::string(SEQ_LEN, 'I'));

    fs::remove(path);
}


TEST_CASE("FASTQ — empty sequence line (zero-length reads) — multithreaded (>10MB)",
          "[strategy][fastq][bugB][multithreaded]") {
    const size_t SEQ_LEN = 200;
    const size_t N = 40000;
    const size_t EMPTY_EVERY = 1001; // every 1001st read is zero-length
    std::string path = "bugB_mt_emptyseq.fastq";
    {
        std::ofstream out(path, std::ios::binary);
        for (size_t i = 0; i < N; ++i) {
            out << "@read" << i << "\n";
            if (i % EMPTY_EVERY == 0) {
                out << "\n+\n\n"; // zero-length read
            } else {
                for (size_t j = 0; j < SEQ_LEN; ++j) out << "ACGT"[j % 4];
                out << "\n+\n";
                for (size_t j = 0; j < SEQ_LEN; ++j) out << 'I';
                out << "\n";
            }
        }
    }
    REQUIRE(fs::file_size(path) > 10u * 1024u * 1024u);

    TracEon::SmartStrategy strategy;
    strategy.loadFile(path);

    // Zero-length reads are not stored (existing invariant), but the cycle
    // must NOT break: all non-empty records survive (was: whole file lost
    // from the first zero-length read onward).
    size_t empty_count = 0;
    for (size_t i = 0; i < N; ++i) if (i % EMPTY_EVERY == 0) ++empty_count;
    REQUIRE(strategy.getFileCacheSize() == N - empty_count);
    REQUIRE_FALSE(strategy.hasSequence("read0"));      // zero-length
    REQUIRE(strategy.getSequence("read1").size() == SEQ_LEN);
    REQUIRE(strategy.getQuality("read1") == std::string(SEQ_LEN, 'I'));
    REQUIRE_FALSE(strategy.hasSequence("read1001"));   // zero-length
    REQUIRE(strategy.getSequence("read1002").size() == SEQ_LEN);
    REQUIRE(strategy.getSequence("read1000").size() == SEQ_LEN);
    const std::string last = "read" + std::to_string(N - 1);
    REQUIRE(strategy.hasSequence(last));
    REQUIRE(strategy.getSequence(last).size() == SEQ_LEN);
    REQUIRE(strategy.getQuality(last) == std::string(SEQ_LEN, 'I'));

    fs::remove(path);
}
