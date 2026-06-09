#include <catch2/catch_test_macros.hpp>
#include "SmartStrategy.h"
#include "SimdUtils.h"
#include <fstream>
#include <filesystem>
#include <thread>
#include <atomic>
#include <zlib.h>

namespace fs = std::filesystem;

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

TEST_CASE("GZIP Decompression Correctness", "[strategy][gzip]") {
    std::string original = ">seq1\nACGT\n>seq2\nTGCA\n";
    std::string gz_filename = "temp_test.fasta.gz";

    gzFile file = gzopen(gz_filename.c_str(), "wb");
    REQUIRE(file != nullptr);
    int written = gzwrite(file, original.data(), (unsigned int)original.size());
    REQUIRE(written == (int)original.size());
    gzclose(file);

    TracEon::SmartStrategy strategy;
    
    SECTION("Can explicit load .gz file") {
        strategy.loadGzipFile(gz_filename);
        REQUIRE(strategy.getFileCacheSize() == 2);
        REQUIRE(strategy.getSequence("seq1") == "ACGT");
    }

    SECTION("Can auto-detect .gz extension") {
        strategy.loadFile(gz_filename);
        REQUIRE(strategy.getFileCacheSize() == 2);
        REQUIRE(strategy.getSequence("seq1") == "ACGT");
    }
    
    if (fs::exists(gz_filename)) fs::remove(gz_filename);
}

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

TEST_CASE("loadFile throws on missing file", "[strategy][errors]") {
    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadFile("nonexistent_file_xyz.fasta"), std::runtime_error);
}

TEST_CASE("loadBinary throws on missing file", "[strategy][errors]") {
    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadBinary("nonexistent_cache_xyz.bin"), std::runtime_error);
}

TEST_CASE("loadBinary rejects corrupt/truncated file", "[strategy][errors]") {
    std::string bad_file = "corrupt.bin";
    {
        std::ofstream out(bad_file, std::ios::binary);
        // Write valid magic but then truncate immediately
        out.write("TRO\x01", 4);
        out.write("\x00", 1); // mode byte
        // Missing count and records — deliberately truncated
    }
    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadBinary(bad_file), std::runtime_error);
    fs::remove(bad_file);
}

TEST_CASE("loadBinary rejects old MMAP format", "[strategy][errors]") {
    std::string old_file = "old_format.bin";
    {
        std::ofstream out(old_file, std::ios::binary);
        out.write("MMAP", 4); // old v1.0 magic
        uint8_t mode = 0; out.write(reinterpret_cast<const char*>(&mode), 1);
        uint64_t count = 0; out.write(reinterpret_cast<const char*>(&count), 8);
    }
    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadBinary(old_file), std::runtime_error);
    fs::remove(old_file);
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
               "        >seq2\nGGGGCCCC";  // <-- no trailing newline
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

// ═══════════════════════════════════════════════════════════════════════════
// Multithreaded SIMD Parser Tests
// ═══════════════════════════════════════════════════════════════════════════
//
// These tests exercise parseFastaMultithreadedTemplate() and
// parseFastqMultithreadedTemplate() via the public loadFile() API by creating
// content larger than MULTITHREAD_THRESHOLD (10 MB).
//
// The threshold check lives in parseArena():
//   if (content.size() > 10*1024*1024) {
//       if (isFastq) parseFastqMultithreadedTemplate(...);
//       else         parseFastaMultithreadedTemplate(...);
//   }
//
// ── Helpers ───────────────────────────────────────────────────────────────────

namespace {
    // Build a FASTA string larger than 10 MB with record_count records.
    // Each record: ">record_N\n{pattern}\n"
    // pattern is repeated_chunk repeated to fill approx seq_len bytes per record.
    std::string make_large_fasta(size_t target_bytes, size_t seq_len = 4000) {
        const size_t REPEATS = seq_len / 4;
        const std::string chunk = "ACGT";
        std::string seq;
        seq.reserve(seq_len);
        for (size_t i = 0; i < REPEATS; ++i) seq += chunk;

        std::string result;
        result.reserve(target_bytes + 4096);
        int n = 1;
        while (result.size() < target_bytes) {
            result += ">record_" + std::to_string(n) + "\n";
            result += seq;
            result += '\n';
            ++n;
        }
        return result;
    }

    // Build a FASTQ string larger than 10 MB with record_count records.
    std::string make_large_fastq(size_t target_bytes, size_t seq_len = 4000) {
        const size_t REPEATS = seq_len / 4;
        const std::string seq = std::string(REPEATS, 'A') + std::string(REPEATS, 'C') +
                                std::string(REPEATS, 'G') + std::string(REPEATS, 'T');
        const std::string qual(seq_len, 'I');

        std::string result;
        result.reserve(target_bytes + 4096);
        int n = 1;
        while (result.size() < target_bytes) {
            result += "@record_" + std::to_string(n) + "\n";
            result += seq;
            result += '\n';
            result += "+\n";
            result += qual;
            result += '\n';
            ++n;
        }
        return result;
    }
} // anonymous namespace

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
