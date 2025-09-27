#include <catch2/catch_all.hpp>
#include "../Cache.h"
#include "../DecodedRecordTypes.h" // This file doesn't exist yet

TEST_CASE("FASTQ File Loading", "[fastq]") {

    SECTION("Can load a simple FASTQ file") {
        TracEon::Cache cache;
        std::string test_file_path = "../test_data/simple.fastq";

        // This method is not implemented yet.
        cache.loadFastq(test_file_path);

        REQUIRE(cache.size() == 2);

        // We'll need a new way to get a full FASTQ record.
        auto record1 = cache.getFastqRecord("seq1");
        REQUIRE(record1.has_value());
        REQUIRE(record1->sequence == "GATTACA");
        REQUIRE(record1->quality == "!''*.~~");

        auto record2 = cache.getFastqRecord("seq2");
        REQUIRE(record2.has_value());
        REQUIRE(record2->sequence == "TTAACCGG");
        REQUIRE(record2->quality == "!''*+,-.");
    }
}
//
// Created by Adnan Raza on 9/27/2025.
//