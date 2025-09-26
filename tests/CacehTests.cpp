#include <catch2/catch_all.hpp>
#include "../Cache.h"

TEST_CASE("Cache Initialization", "[cache]") {
    SECTION("A newly created cache is empty") {
        TracEon::Cache cache;
        REQUIRE(cache.size() == 0);
    }
}

//
// Created by Adnan Raza on 9/26/2025.
//