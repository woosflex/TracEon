#include <catch2/catch_all.hpp>
#include "../Cache.h"

TEST_CASE("Cache Initialization", "[cache]") {
    SECTION("A newly created cache is empty") {
        TracEon::Cache cache;
        REQUIRE(cache.size() == 0);
    }

    SECTION("Can set and get a simple key-value pair") {
        TracEon::Cache cache;
        std::string key = "seq1";
        std::string value = "GATTACA";

        // These methods don't exist yet, so this will fail to compile.
        cache.set(key, value);
        std::string retrieved_value = cache.get(key);

        REQUIRE(retrieved_value == value);
        REQUIRE(cache.size() == 1);
    }
}

//
// Created by Adnan Raza on 9/26/2025.
//