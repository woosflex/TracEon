#include <catch2/catch_all.hpp>
#include "../Cache.h"

TEST_CASE("Cache Functionality", "[cache]") {
    SECTION("A newly created cache is empty") {
        TracEon::Cache cache;
        REQUIRE(cache.size() == 0);
    }

    SECTION("Can set and get a simple key-value pair") {
        TracEon::Cache cache;
        std::string key = "seq1";
        std::string value = "GATTACA";

        cache.set(key, value);
        std::string retrieved_value = cache.get(key);

        REQUIRE(retrieved_value == value);
        REQUIRE(cache.size() == 1);
    }
    SECTION("Can switch strategy and get real compression") {
        TracEon::Cache cache;
        std::string key = "seq1";
        std::string value = "GATTACA"; // 7 bytes

        cache.setStrategy(TracEon::StrategyType::TwoBit);
        cache.set(key, value);

        // We need a way to check the size of the stored data.
        size_t stored_size = cache.getStoredSize(key);

        // 4 bytes for length header + 2 bytes for data = 6 bytes
        REQUIRE(stored_size == 6);
        REQUIRE(stored_size < value.size());
    }

    SECTION("Can load a simple FASTA file") {
        TracEon::Cache cache;
        // The path is relative to the build directory, so we go up one level.
        std::string test_file_path = "../test_data/simple.fasta";

        // This method is not implemented yet, so this will cause a linker error.
        cache.loadFasta(test_file_path);

        REQUIRE(cache.size() == 2);
        REQUIRE(cache.get("seq1") == "GATTACA");
        REQUIRE(cache.get("seq2") == "CGCGCGCGCGCGCGCGCGCGCGCGCGCG");
    }
}