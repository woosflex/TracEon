#include <catch2/catch_test_macros.hpp>
#include "SmartStrategy.h"
#include "SimdUtils.h"
#include <fstream>
#include <filesystem>
#include <thread>
#include <atomic>
#include <zlib.h>
#include <lz4.h>
#include <cstring>
#include <cstdint>
#ifndef _WIN32
#include <sys/resource.h>
#endif

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

TEST_CASE("addEntry() round-trips through saveBinary/loadBinary",
          "[strategy][addentry][smart_compression]") {
    // Under the immutable-after-load contract (ADR-001, Bug 3 fix), addEntry()
    // is only legal before any load or after clearCache(), so a manual
    // build can no longer be mixed with a loadFile() into one object. This
    // still exercises addEntry()'s incremental serialized_size_estimate_
    // updates on top of the sizeof(uint64_t) record-count base, and the
    // loadBinary→saveBinary path (refreshPayloadEstimate seeding) is covered
    // by the Smart Compression round-trip tests below.
    const std::string bin = "addentry_roundtrip.bin";

    TracEon::SmartStrategy strategy;
    strategy.addEntry("manual1", "GATTACA", "IIIIIII");
    strategy.addEntry("manual2", "CATCATCAT", "");
    REQUIRE(strategy.getFileCacheSize() == 2);

    strategy.saveBinary(bin);

    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getFileCacheSize() == 2);
    REQUIRE(restored.getSequence("manual1") == "GATTACA");
    REQUIRE(restored.getQuality("manual1") == "IIIIIII");
    REQUIRE(restored.getSequence("manual2") == "CATCATCAT");

    fs::remove(bin);
}

TEST_CASE("loadFile throws on missing file", "[strategy][errors]") {
    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadFile("nonexistent_file_xyz.fasta"), std::runtime_error);
}

// ── Smart Compression tests ───────────────────────────────────────────────────

TEST_CASE("Smart Compression — large DNA payload uses LZ4_HC (v3, better ratio)",
          "[strategy][smart_compression][lz4hc]") {
    const std::string src = "sc_large_dna.fasta";
    const std::string bin = "sc_large_dna.bin";
    {
        std::ofstream out(src);
        // 3000 × 4000-char sequences ≈ 12 MiB payload — above the 10 MiB HC threshold
        const std::string seq(4000, 'A');
        for (int i = 1; i <= 3000; ++i)
            out << ">seq" << i << '\n' << seq << '\n';
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTA);
    strategy.saveBinary(bin);

    // Verify v3 magic (LZ4_HC and LZ4_default both write the streaming v3 format)
    {
        std::ifstream f(bin, std::ios::binary);
        char magic[4];
        f.read(magic, 4);
        REQUIRE(magic[3] == '\x03');
    }

    // LZ4_HC on homopolymer runs achieves extreme ratios — expect very small file
    REQUIRE(fs::file_size(bin) < 500'000); // < 500 KB for ~12 MiB of AAAA...

    // Full round-trip correctness
    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getFileCacheSize() == 3000);
    REQUIRE(restored.getSequence("seq1")    == std::string(4000, 'A'));
    REQUIRE(restored.getSequence("seq3000") == std::string(4000, 'A'));

    fs::remove(src);
    fs::remove(bin);
}

TEST_CASE("Smart Compression — large protein payload uses LZ4_default (not HC)",
          "[strategy][smart_compression][lz4default_protein]") {
    const std::string src = "sc_large_protein.fasta";
    const std::string bin = "sc_large_protein.bin";
    {
        std::ofstream out(src);
        // 3000 × 4000-char protein sequences ≈ 12 MiB payload
        // 'E' (glutamic acid) is not in the nucleotide set — forces PROTEIN_FASTA
        const std::string seq(4000, 'E');
        for (int i = 1; i <= 3000; ++i)
            out << ">prot" << i << '\n' << seq << '\n';
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::PROTEIN_FASTA);
    strategy.saveBinary(bin);

    // v3 magic expected (payload is above 64 KiB, so not uncompressed)
    {
        std::ifstream f(bin, std::ios::binary);
        char magic[4];
        f.read(magic, 4);
        REQUIRE(magic[3] == '\x03');
    }

    // Round-trip correctness
    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getFileCacheSize() == 3000);
    REQUIRE(restored.getSequence("prot1")    == std::string(4000, 'E'));
    REQUIRE(restored.getSequence("prot3000") == std::string(4000, 'E'));

    fs::remove(src);
    fs::remove(bin);
}

TEST_CASE("Smart Compression — small DNA payload uses LZ4_default (below HC threshold)",
          "[strategy][smart_compression][lz4default_small]") {
    const std::string src = "sc_small_dna.fasta";
    const std::string bin = "sc_small_dna.bin";
    {
        std::ofstream out(src);
        out << ">chr1\nACGT\n>chr2\nTGCA\n";
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTA);
    strategy.saveBinary(bin);

    // v3 magic: LZ4_default path still writes the streaming v3 format
    {
        std::ifstream f(bin, std::ios::binary);
        char magic[4];
        f.read(magic, 4);
        REQUIRE(magic[3] == '\x03');
    }

    // Round-trip correctness
    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getFileCacheSize() == 2);
    REQUIRE(restored.getSequence("chr1") == "ACGT");
    REQUIRE(restored.getSequence("chr2") == "TGCA");

    fs::remove(src);
    fs::remove(bin);
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

// ── Backward compatibility: v1/v2 binary cache readers ──────────────────────
// saveBinary() now always writes v3 (streaming LZ4 Frame), but loadBinary()
// must still read caches written by older versions of TracEon. These tests
// hand-construct valid v1 (uncompressed) and v2 (single-block LZ4) binaries
// using the same layout the old writers produced, to verify the reader
// branches for those formats were not broken by the v3 changes.

static void write_genome_payload(std::vector<char>& buf, const std::vector<std::tuple<std::string, std::string, std::string>>& records) {
    uint64_t count = records.size();
    buf.insert(buf.end(), reinterpret_cast<const char*>(&count), reinterpret_cast<const char*>(&count) + sizeof(count));
    for (const auto& [id, seq, qual] : records) {
        auto append_field = [&](const std::string& field) {
            uint32_t len = static_cast<uint32_t>(field.size());
            buf.insert(buf.end(), reinterpret_cast<const char*>(&len), reinterpret_cast<const char*>(&len) + sizeof(len));
            buf.insert(buf.end(), field.data(), field.data() + len);
        };
        append_field(id);
        append_field(seq);
        append_field(qual);
    }
}

TEST_CASE("loadBinary reads legacy v1 (uncompressed) format", "[strategy][compat][v1]") {
    std::string bin = "legacy_v1.bin";
    std::vector<std::tuple<std::string, std::string, std::string>> records = {
        {"seq1", "ACGT", ""}, {"seq2", "TGCA", "IIII"}
    };
    {
        std::vector<char> payload;
        write_genome_payload(payload, records);

        std::ofstream out(bin, std::ios::binary);
        out.write("TRO\x01", 4);
        uint8_t mode = 0;
        out.write(reinterpret_cast<const char*>(&mode), 1);
        out.write(payload.data(), static_cast<std::streamsize>(payload.size()));
    }

    TracEon::SmartStrategy strategy;
    strategy.loadBinary(bin);
    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("seq1") == "ACGT");
    REQUIRE(strategy.getSequence("seq2") == "TGCA");
    REQUIRE(strategy.getQuality("seq2") == "IIII");

    fs::remove(bin);
}

TEST_CASE("loadBinary reads legacy v2 (single-block LZ4) format", "[strategy][compat][v2]") {
    std::string bin = "legacy_v2.bin";
    std::vector<std::tuple<std::string, std::string, std::string>> records = {
        {"seq1", "ACGT", ""}, {"seq2", "TGCA", "IIII"}
    };
    std::vector<char> payload;
    write_genome_payload(payload, records);

    const size_t max_compressed = static_cast<size_t>(LZ4_compressBound(static_cast<int>(payload.size())));
    std::vector<char> compressed(max_compressed);
    int compressed_size = LZ4_compress_default(payload.data(), compressed.data(),
                                                static_cast<int>(payload.size()),
                                                static_cast<int>(max_compressed));
    REQUIRE(compressed_size > 0);
    compressed.resize(static_cast<size_t>(compressed_size));

    {
        std::ofstream out(bin, std::ios::binary);
        out.write("TRO\x02", 4);
        uint8_t mode = 0;
        out.write(reinterpret_cast<const char*>(&mode), 1);
        uint64_t original_size = static_cast<uint64_t>(payload.size());
        uint64_t compressed_len = static_cast<uint64_t>(compressed_size);
        out.write(reinterpret_cast<const char*>(&original_size), sizeof(original_size));
        out.write(reinterpret_cast<const char*>(&compressed_len), sizeof(compressed_len));
        out.write(compressed.data(), compressed_size);
    }

    TracEon::SmartStrategy strategy;
    strategy.loadBinary(bin);
    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("seq1") == "ACGT");
    REQUIRE(strategy.getSequence("seq2") == "TGCA");
    REQUIRE(strategy.getQuality("seq2") == "IIII");

    fs::remove(bin);
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

// Helper: compress bytes with zlib-ng gzwrite into a temp file, read back raw bytes
static std::vector<char> compress_to_gzip_bytes(const std::string& data, const std::string& tmp_path) {
    gzFile f = gzopen(tmp_path.c_str(), "wb");
    gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
    gzclose(f);

    std::ifstream in(tmp_path, std::ios::binary | std::ios::ate);
    size_t sz = in.tellg();
    in.seekg(0);
    std::vector<char> bytes(sz);
    in.read(bytes.data(), sz);
    return bytes;
}

TEST_CASE("Parallel GZIP decompression — concatenated streams", "[strategy][gzip][parallel]") {
    std::string part1 = ">seq1\nACGTACGT\n>seq2\nTGCATGCA\n";
    std::string part2 = ">seq3\nAAAAAAAA\n>seq4\nCCCCCCCC\n";

    std::string tmp1 = "tmp_stream1.gz";
    std::string tmp2 = "tmp_stream2.gz";
    std::string concat_path = "tmp_concat.fasta.gz";

    auto bytes1 = compress_to_gzip_bytes(part1, tmp1);
    auto bytes2 = compress_to_gzip_bytes(part2, tmp2);

    {
        std::ofstream out(concat_path, std::ios::binary);
        out.write(bytes1.data(), bytes1.size());
        out.write(bytes2.data(), bytes2.size());
    }

    SECTION("All sequences present after concatenated load") {
        TracEon::SmartStrategy strategy;
        strategy.loadFile(concat_path);

        REQUIRE(strategy.getFileCacheSize() == 4);
        REQUIRE(strategy.getSequence("seq1") == "ACGTACGT");
        REQUIRE(strategy.getSequence("seq2") == "TGCATGCA");
        REQUIRE(strategy.getSequence("seq3") == "AAAAAAAA");
        REQUIRE(strategy.getSequence("seq4") == "CCCCCCCC");
    }

    fs::remove(tmp1);
    fs::remove(tmp2);
    fs::remove(concat_path);
}

TEST_CASE("Parallel GZIP — single-stream fallback unchanged", "[strategy][gzip][parallel]") {
    std::string data = ">chr1\nACGTACGTACGT\n>chr2\nTTTTGGGGCCCC\n";
    std::string gz_path = "tmp_single_stream.fasta.gz";

    gzFile f = gzopen(gz_path.c_str(), "wb");
    gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
    gzclose(f);

    TracEon::SmartStrategy strategy;
    strategy.loadFile(gz_path);

    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("chr1") == "ACGTACGTACGT");
    REQUIRE(strategy.getSequence("chr2") == "TTTTGGGGCCCC");

    fs::remove(gz_path);
}

TEST_CASE("Parallel GZIP — two large streams actually exercise loadGzipParallel", "[strategy][gzip][parallel]") {
    // Each compressed stream must exceed PARALLEL_GZIP_THRESHOLD/2 (512 KB) so
    // the concatenated file is > 1 MB, triggering loadGzipParallel.
    // FASTA nucleotide sequences compress very well with zlib (50-200:1 for
    // periodic patterns). Solution: use gzip level 0 (store mode), which wraps
    // data in a valid GZIP container without compression. This keeps the file
    // size close to the uncompressed size, guaranteeing the threshold is crossed.

    const int SEQ_LEN = 800;
    const char BASES[] = "ACGT";

    auto make_fasta = [&](const std::string& prefix, int count) {
        std::string s;
        for (int i = 0; i < count; ++i) {
            s += '>';
            s += prefix;
            s += std::to_string(i);
            s += '\n';
            for (int j = 0; j < SEQ_LEN; ++j)
                s += BASES[(i + j) % 4];
            s += '\n';
        }
        return s;
    };

    // Each stream: 1500 records × ~800 bytes ≈ 1.2 MB.
    // At gzip level 0 (no compression), file size ≈ uncompressed ≈ 1.2 MB per stream.
    // Total: ~2.4 MB → well above PARALLEL_GZIP_THRESHOLD (1 MB).
    std::string part1 = make_fasta("stream1_seq", 1500);
    std::string part2 = make_fasta("stream2_seq", 1500);

    std::string tmp1 = "tmp_large_stream1.gz";
    std::string tmp2 = "tmp_large_stream2.gz";
    std::string concat_path = "tmp_large_concat.fasta.gz";

    // "wb0" = gzip store mode (level 0): valid GZIP format, no compression.
    auto compress_stream = [](const std::string& data, const std::string& path) {
        gzFile f = gzopen(path.c_str(), "wb0");
        REQUIRE(f != nullptr);
        gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
        gzclose(f);
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        size_t sz = static_cast<size_t>(in.tellg());
        in.seekg(0);
        std::vector<char> bytes(sz);
        in.read(bytes.data(), sz);
        return bytes;
    };

    auto bytes1 = compress_stream(part1, tmp1);
    auto bytes2 = compress_stream(part2, tmp2);

    {
        std::ofstream out(concat_path, std::ios::binary);
        out.write(bytes1.data(), bytes1.size());
        out.write(bytes2.data(), bytes2.size());
    }

    // Guard: concatenated file must exceed threshold so this test is meaningful.
    size_t concat_size = static_cast<size_t>(
        std::filesystem::file_size(concat_path));
    REQUIRE(concat_size > 1024 * 1024); // must trigger loadGzipParallel

    TracEon::SmartStrategy strategy;
    strategy.loadFile(concat_path);

    REQUIRE(strategy.getFileCacheSize() == 3000);
    REQUIRE(strategy.hasSequence("stream1_seq0"));
    REQUIRE(static_cast<int>(strategy.getSequence("stream1_seq0").size()) == SEQ_LEN);
    REQUIRE(strategy.hasSequence("stream1_seq1499"));
    REQUIRE(strategy.hasSequence("stream2_seq0"));
    REQUIRE(static_cast<int>(strategy.getSequence("stream2_seq0").size()) == SEQ_LEN);
    REQUIRE(strategy.hasSequence("stream2_seq1499"));

    fs::remove(tmp1);
    fs::remove(tmp2);
    fs::remove(concat_path);
}

TEST_CASE("Parallel GZIP — coincidental magic bytes in a single-member stream don't cause a false split",
          "[strategy][gzip][parallel]") {
    // Regression test: scanGzipStreams() used to look for the raw byte
    // sequence 0x1f 0x8b 0x08 anywhere in the file to find concatenated
    // stream boundaries. That sequence can appear by chance *inside* a
    // single member's compressed payload (observed on real NCBI FASTA.gz
    // files). Splitting there feeds inflate() a misaligned/truncated
    // deflate stream, which can explode into gigabytes of garbage via LZ77
    // back-references — tripping the OOM guard on a perfectly valid file.
    //
    // Level 0 (store mode) gzip embeds input bytes close to verbatim, so
    // planting 0x1f 0x8b 0x08 in the plaintext reliably reproduces a
    // coincidental match in the compressed output.
    const int SEQ_LEN = 800;
    const char BASES[] = "ACGT";

    auto make_fasta = [&](const std::string& prefix, int count) {
        std::string s;
        for (int i = 0; i < count; ++i) {
            s += '>';
            s += prefix;
            s += std::to_string(i);
            s += '\n';
            for (int j = 0; j < SEQ_LEN; ++j)
                s += BASES[(i + j) % 4];
            s += '\n';
        }
        return s;
    };

    std::string content = make_fasta("bogus_magic_seq", 1500); // ~1.2 MB, single member

    // Plant a coincidental gzip-header-looking byte triple in the middle of
    // a sequence line (not on a newline boundary, so it can't be mistaken
    // for record structure).
    size_t plant_at = content.size() / 2;
    content[plant_at]     = static_cast<char>(0x1f);
    content[plant_at + 1] = static_cast<char>(0x8b);
    content[plant_at + 2] = static_cast<char>(0x08);

    std::string gz_path = "tmp_false_positive_magic.fasta.gz";
    gzFile f = gzopen(gz_path.c_str(), "wb0"); // store mode: preserves planted bytes verbatim
    REQUIRE(f != nullptr);
    gzwrite(f, content.data(), static_cast<unsigned>(content.size()));
    gzclose(f);

    size_t gz_size = static_cast<size_t>(std::filesystem::file_size(gz_path));
    REQUIRE(gz_size > 1024 * 1024); // must trigger loadGzipParallel's dispatch path

    TracEon::SmartStrategy strategy;
    REQUIRE_NOTHROW(strategy.loadFile(gz_path));

    REQUIRE(strategy.getFileCacheSize() == 1500);
    REQUIRE(strategy.hasSequence("bogus_magic_seq0"));
    REQUIRE(static_cast<int>(strategy.getSequence("bogus_magic_seq0").size()) == SEQ_LEN);
    REQUIRE(strategy.hasSequence("bogus_magic_seq1499"));
    REQUIRE(static_cast<int>(strategy.getSequence("bogus_magic_seq1499").size()) == SEQ_LEN);

    fs::remove(gz_path);
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

// ─── v1 Binary Format ────────────────────────────────────────────────────────

TEST_CASE("v1 binary format (TRO\\x01) load and retrieve", "[strategy][binary]") {
    // Construct a valid v1 blob in-memory and write to disk.
    // Format: magic(4) | mode(1) | count(8) | [id_len(4) id seq_len(4) seq qual_len(4)] ...

    auto append_u32 = [](std::vector<char>& buf, uint32_t v) {
        buf.insert(buf.end(), reinterpret_cast<const char*>(&v),
                               reinterpret_cast<const char*>(&v) + 4);
    };
    auto append_u64 = [](std::vector<char>& buf, uint64_t v) {
        buf.insert(buf.end(), reinterpret_cast<const char*>(&v),
                               reinterpret_cast<const char*>(&v) + 8);
    };
    auto append_str = [](std::vector<char>& buf, std::string_view s) {
        buf.insert(buf.end(), s.begin(), s.end());
    };

    std::vector<char> blob;
    // Magic + version
    const char magic[] = {'T','R','O','\x01'};
    blob.insert(blob.end(), magic, magic + 4);
    // mode = 0 (GenomeIndex)
    blob.push_back(0);
    // count = 2
    append_u64(blob, 2);
    // record 1: "chr1" / "ACGTACGT" / ""
    append_u32(blob, 4); append_str(blob, "chr1");
    append_u32(blob, 8); append_str(blob, "ACGTACGT");
    append_u32(blob, 0);
    // record 2: "chr2" / "TTTTGGGG" / ""
    append_u32(blob, 4); append_str(blob, "chr2");
    append_u32(blob, 8); append_str(blob, "TTTTGGGG");
    append_u32(blob, 0);

    std::string bin_path = "tmp_v1_format.bin";
    {
        std::ofstream out(bin_path, std::ios::binary);
        out.write(blob.data(), blob.size());
    }

    TracEon::SmartStrategy strategy;
    strategy.loadBinary(bin_path);

    REQUIRE(strategy.getIndexMode() == TracEon::IndexMode::GENOME);
    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("chr1") == "ACGTACGT");
    REQUIRE(strategy.getSequence("chr2") == "TTTTGGGG");
    REQUIRE_FALSE(strategy.hasSequence("chr3"));

    fs::remove(bin_path);
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
