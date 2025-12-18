#include <catch2/catch_all.hpp>
#include <fstream>
#include <filesystem>
#include "Cache.h"

TEST_CASE("Cache Functionality", "[cache]") {

    SECTION("Can set and get a simple key-value pair") {
        TracEon::Cache cache;
        std::string key = "seq1";
        std::string value = "GATTACA";

        cache.set(key, value);
        std::string retrieved_value = cache.get(key);

        REQUIRE(retrieved_value == value);
        REQUIRE(cache.size() == 1);
    }

    SECTION("Can load a simple FASTA file") {
        TracEon::Cache cache;
        std::string test_file_path = "simple.fasta"; // Create locally for test
        
        {
             std::ofstream out(test_file_path);
             out << ">seq1\nGATTACA\n>seq2\nCGCGCGCGCGCGCGCGCGCGCGCGCGCG\n";
             out.close();
        }

        cache.loadFile(test_file_path);

        REQUIRE(cache.size() == 2);
        REQUIRE(cache.get("seq1") == "GATTACA");
        REQUIRE(cache.get("seq2") == "CGCGCGCGCGCGCGCGCGCGCGCGCGCG");
        
        std::filesystem::remove(test_file_path);
    }

    SECTION("Can save and restore the cache state") {
        std::string src_filename = "temp_persist.fasta";
        std::string bin_filename = "test_cache.bin";

        // 1. Create a source file (Architecturally correct source of truth)
        {
            std::ofstream out(src_filename);
            out << ">seq1\nGATTACA\n>seq2\nCGCGCGCG\n";
            out.close();
        }

        // 2. Load into Cache 1 (Populates SmartStrategy/Arena)
        TracEon::Cache cache1;
        cache1.loadFile(src_filename);
        REQUIRE(cache1.size() == 2);

        // 3. Save to Binary (Serializes SmartStrategy)
        cache1.save(bin_filename);

        // 4. Restore into Cache 2 (Maps Binary)
        TracEon::Cache cache2;
        cache2.restore(bin_filename);

        // 5. Verify Persistence
        REQUIRE(cache2.size() == cache1.size());
        REQUIRE(cache2.get("seq1") == "GATTACA");
        REQUIRE(cache2.get("seq2") == "CGCGCGCG");
        
        // Cleanup
        std::filesystem::remove(src_filename);
        std::filesystem::remove(bin_filename);
    }
}