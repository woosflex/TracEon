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
        std::string value = "GATTACA"; // 7 bytes long

        cache.setStrategy(TracEon::StrategyType::TwoBit);
        cache.set(key, value);

        size_t stored_size = cache.getStoredSize(key);

        // Expect 4 bytes for length header + 2 bytes for data = 6 bytes
        REQUIRE(stored_size == 6);
        REQUIRE(stored_size < value.size());
    }

    SECTION("Can load a simple FASTA file") {
        TracEon::Cache cache;
        // This path is relative to the build directory in cmake-build-debug-wsl/
        std::string test_file_path = "../test_data/simple.fasta";

        cache.loadFasta(test_file_path);

        REQUIRE(cache.size() == 2);
        REQUIRE(cache.get("seq1") == "GATTACA");
        REQUIRE(cache.get("seq2") == "CGCGCGCGCGCGCGCGCGCGCGCGCGCG");
    }
    SECTION("Can save and restore the cache state") {
        // Use a temporary filename for our test.
        std::string test_filename = "test_cache.bin";

        // 1. Create and populate an original cache.
        TracEon::Cache cache1;
        cache1.setStrategy(TracEon::StrategyType::TwoBit);
        cache1.set("seq1", "GATTACA");
        cache1.set("seq2", "CGCGCGCG");

        // This method is not implemented yet.
        cache1.save(test_filename);

        // 2. Create a new, empty cache.
        TracEon::Cache cache2;

        // 3. Restore the state from the file.
        cache2.restore(test_filename);

        // 4. Verify that the restored cache is identical to the original.
        REQUIRE(cache2.size() == cache1.size());
        REQUIRE(cache2.get("seq1") == "GATTACA");
        REQUIRE(cache2.get("seq2") == "CGCGCGCG");
        // Also check that the restored cache knows which strategy was used.
        REQUIRE(cache2.getStoredSize("seq1") == 6);
    }
}