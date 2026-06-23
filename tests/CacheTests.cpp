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
        std::string test_file_path = "simple.fasta";
        
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

    SECTION("Can load a multi-line FASTA file") {
        TracEon::Cache cache;
        std::string test_file_path = "multiline_cache.fasta";

        {
            std::ofstream out(test_file_path);
            // Two-line wrapped sequence (classic NCBI style)
            out << ">seq1\nACGT\nACGT\n>seq2\nTTTT\n";
            out.close();
        }

        cache.loadFile(test_file_path);

        REQUIRE(cache.size() == 2);
        REQUIRE(cache.get("seq1") == "ACGTACGT");
        REQUIRE(cache.get("seq2") == "TTTT");

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

    SECTION("LZ4 binary cache compression and round-trip integrity") {
        std::string src_filename = "temp_lz4_test.fasta";
        std::string bin_filename = "test_cache_lz4.bin";

        // 1. Create a source file with repetitive data (compresses well)
        {
            std::ofstream out(src_filename);
            out << ">seq1\n";
            for (int i = 0; i < 100; ++i) out << "ACGTACGTACGTACGT";
            out << "\n>seq2\n";
            for (int i = 0; i < 100; ++i) out << "GCTAGCTAGCTAGCTA";
            out << "\n";
            out.close();
        }

        // 2. Load and save to v2 binary format (LZ4-compressed)
        TracEon::Cache cache1;
        cache1.loadFile(src_filename);
        REQUIRE(cache1.size() == 2);

        cache1.save(bin_filename);
        size_t compressed_size = std::filesystem::file_size(bin_filename);

        // 3. Restore from v2 binary
        TracEon::Cache cache2;
        cache2.restore(bin_filename);

        // 4. Verify all sequences match exactly
        REQUIRE(cache2.size() == 2);
        REQUIRE(cache2.get("seq1").length() == 1600);
        REQUIRE(cache2.get("seq2").length() == 1600);

        // Verify content
        std::string seq1_retrieved = cache2.get("seq1");
        std::string seq1_expected(1600, '\0');
        for (int i = 0; i < 100; ++i) {
            seq1_expected.replace(i * 16, 16, "ACGTACGTACGTACGT");
        }
        REQUIRE(seq1_retrieved == seq1_expected);

        // 5. Verify file size is reasonable (compressed < uncompressed)
        std::ifstream src_check(src_filename);
        src_check.seekg(0, std::ios::end);
        size_t src_size = src_check.tellg();
        src_check.close();

        // Compressed should be significantly smaller than source
        REQUIRE(compressed_size < src_size);

        // Cleanup
        std::filesystem::remove(src_filename);
        std::filesystem::remove(bin_filename);
    }
}