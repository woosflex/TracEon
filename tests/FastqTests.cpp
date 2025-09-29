#include <catch2/catch_all.hpp>
#include <fstream>
#include "../Cache.h"
#include "../DecodedRecordTypes.h"

TEST_CASE("FASTQ File Loading", "[fastq]") {
    SECTION("Can load a simple FASTQ file") {
        TracEon::Cache cache;
        std::string test_file_path = "../test_data/simple.fastq";

        cache.loadFile(test_file_path);

        REQUIRE(cache.size() == 2);

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