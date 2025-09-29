#include <catch2/catch_all.hpp>
#include <fstream>
#include "../Cache.h"

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

    SECTION("Can get real compression with SmartStrategy") {
        TracEon::Cache cache;
        std::string key = "seq1";
        // Use a longer sequence to see compression benefits
        std::string value = "GATTACAGATTACAGATTACAGATTACA"; // 28 bytes

        cache.set(key, value);

        size_t stored_size = cache.getStoredSize(key);

        // SmartStrategy automatically compresses DNA
        REQUIRE(stored_size < value.size());
    }

    SECTION("Can load a simple FASTA file") {
        TracEon::Cache cache;
        std::string test_file_path = "../test_data/simple.fasta";
        cache.loadFile(test_file_path);

        REQUIRE(cache.size() == 2);
        REQUIRE(cache.get("seq1") == "GATTACA");
        REQUIRE(cache.get("seq2") == "CGCGCGCGCGCGCGCGCGCGCGCGCGCG");
    }

    SECTION("Can save and restore the cache state") {
        std::string test_filename = "test_cache.bin";

        TracEon::Cache cache1;
        // No setStrategy() call needed
        cache1.set("seq1", "GATTACA");
        cache1.set("seq2", "CGCGCGCG");

        cache1.save(test_filename);

        TracEon::Cache cache2;
        cache2.restore(test_filename);

        REQUIRE(cache2.size() == cache1.size());
        REQUIRE(cache2.get("seq1") == "GATTACA");
        REQUIRE(cache2.get("seq2") == "CGCGCGCG");
    }
}