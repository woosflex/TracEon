#include <catch2/catch_test_macros.hpp>
#include "SmartStrategy.h"
#include <fstream>
#include <filesystem>
#include <thread>
#include <atomic>
#include <zlib.h> // FIX: Required for gzFile, gzopen

namespace fs = std::filesystem;

TEST_CASE("SmartStrategy determines format correctly", "[SmartStrategy]") {
    TracEon::SmartStrategy strategy;
    // ... (Assume tests here) ...
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
    REQUIRE(written == original.size());
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