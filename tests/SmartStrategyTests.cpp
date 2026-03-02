#include <catch2/catch_test_macros.hpp>
#include "SmartStrategy.h"
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