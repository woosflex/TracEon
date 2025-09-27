#include <catch2/catch_all.hpp>
#include "../Cache.h"
#include "../DecodedRecordTypes.h"

TEST_CASE("FASTQ File Loading", "[fastq]") {

    SECTION("Can load a simple FASTQ file") {
        TracEon::Cache cache;
        std::string test_file_path = "../test_data/simple.fastq";

        cache.loadFastq(test_file_path);

        REQUIRE(cache.size() == 2);

        // Check the first record
        auto record1 = cache.getFastqRecord("seq1");
        REQUIRE(record1.has_value());
        REQUIRE(record1->sequence == "GATTACA");
        REQUIRE(record1->quality == "!''*.~~");

        // Check the second record
        auto record2 = cache.getFastqRecord("seq2");
        REQUIRE(record2.has_value());
        REQUIRE(record2->sequence == "TTAACCGG");
        REQUIRE(record2->quality == "!''*+,-.");
    }
}