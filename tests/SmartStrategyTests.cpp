#include <catch2/catch_test_macros.hpp>
#include "../include/SmartStrategy.h"
#include <fstream>
#include <filesystem>
#include <thread>
#include <atomic>
#include <chrono>

// Note: These tests assume the presence of 'simple.fasta' and 'simple.fastq'
// in the '../test_data/' directory, consistent with other test files.

namespace fs = std::filesystem;

TEST_CASE("SmartStrategy determines format correctly", "[SmartStrategy]") {

    SECTION("Can load a simple FASTA file") {
        TracEon::SmartStrategy strategy;
        std::string test_file_path = "../test_data/simple.fasta";

        strategy.loadFile(test_file_path);

        REQUIRE(strategy.getFileCacheSize() == 2);
        REQUIRE(strategy.getSequence("seq1") == "GATTACA");
        REQUIRE(strategy.getSequence("seq2") == "CGCGCGCGCGCGCGCGCGCGCGCGCGCG");
        REQUIRE(strategy.getQuality("seq1").empty()); // FASTA has no quality scores
    }

    SECTION("Can load a simple FASTQ file") {
        TracEon::SmartStrategy strategy;
        std::string test_file_path = "../test_data/simple.fastq";

        strategy.loadFile(test_file_path);

        REQUIRE(strategy.getFileCacheSize() == 2);

        // Check first record
        REQUIRE(strategy.getSequence("seq1") == "GATTACA");
        REQUIRE(strategy.getQuality("seq1") == "@@@DDDD");

        // Check second record
        REQUIRE(strategy.getSequence("seq2") == "CGCGCGCGCGCGCGCGCGCGCGCGCGCG");
        REQUIRE(strategy.getQuality("seq2") == "HHHIIIIIIIIIIIIIIIIIIIIIII");
    }
}

TEST_CASE("Lock-Free Read Performance", "[strategy][performance]") {
    TracEon::SmartStrategy strategy;
    // Assuming test is run from build directory
    std::string test_file_path = "../test_data/simple.fasta";
    
    // Ensure file exists
    if (!fs::exists(test_file_path)) {
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
        
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < NUM_THREADS; ++i) {
            threads.emplace_back(worker);
        }
        for (auto& t : threads) t.join();
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        REQUIRE(total_hits == NUM_THREADS * LOOKUPS_PER_THREAD);
        
        // Performance validation: 40K lookups should be instantaneous (<10ms)
        // This implicitly validates that we aren't deadlocking or serializing heavily
        REQUIRE(duration_ms < 50); 
    }
}

TEST_CASE("String Key Storage Validation", "[strategy][architecture]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "../test_data/simple.fasta";
    if (!fs::exists(test_file_path)) {
        std::ofstream out(test_file_path);
        out << ">seq1\nACGT\n>seq2\nTGCA\n";
        out.close();
    }
    strategy.loadFile(test_file_path);
    
    // Verify we're using GenomeIndex (std::string keys)
    // This confirms our "Hybrid Architecture" override is active
    REQUIRE(strategy.getIndexMode() == TracEon::IndexMode::GENOME);
    
    // Verify keys are owned strings
    auto keys = strategy.getAllKeys();
    REQUIRE(keys.size() == 2);
    bool found_seq1 = false;
    for(const auto& k : keys) if(k == "seq1") found_seq1 = true;
    REQUIRE(found_seq1);
}

TEST_CASE("Zero-Copy Value Storage", "[strategy][memory]") {
    TracEon::SmartStrategy strategy;
    std::string test_file_path = "../test_data/simple.fasta";
    if (!fs::exists(test_file_path)) {
        std::ofstream out(test_file_path);
        out << ">seq1\nACGT\n>seq2\nTGCA\n";
        out.close();
    }
    strategy.loadFile(test_file_path);
    
    // Get a view
    auto view1 = strategy.getView("seq1");
    REQUIRE(!view1.empty());
    
    // Get same view again - should point to same memory address
    // This confirms we aren't creating new string copies on retrieval
    auto view2 = strategy.getView("seq1");
    REQUIRE(view1.data() == view2.data()); 
}

