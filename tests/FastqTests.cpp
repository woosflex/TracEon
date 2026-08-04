#include <catch2/catch_all.hpp>
#include <fstream>
#include <filesystem>
#include "Cache.h"
#include "DecodedRecordTypes.h"

TEST_CASE("FASTQ File Loading", "[fastq]") {
    SECTION("Can load a simple FASTQ file") {
        std::string test_file_path = "tmp_simple_test.fastq";

        {
            std::ofstream out(test_file_path);
            out << "@seq1\nGATTACA\n+\n@@@DDDD\n"
                << "@seq2\nCGCGCGCGCGCGCGCGCGCGCGCGCGCG\n+\nHHHIIIIIIIIIIIIIIIIIIIIIII\n";
        }

        TracEon::Cache cache;
        cache.loadFile(test_file_path);

        REQUIRE(cache.size() == 2);

        auto record1 = cache.getFastqRecord("seq1");
        REQUIRE(record1.has_value());
        REQUIRE(record1->sequence == "GATTACA");
        REQUIRE(record1->quality == "@@@DDDD");

        auto record2 = cache.getFastqRecord("seq2");
        REQUIRE(record2.has_value());
        REQUIRE(record2->sequence == "CGCGCGCGCGCGCGCGCGCGCGCGCGCG");
        REQUIRE(record2->quality == "HHHIIIIIIIIIIIIIIIIIIIIIII");

        std::filesystem::remove(test_file_path);
    }
}
