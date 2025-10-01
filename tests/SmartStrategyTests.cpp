#include <catch2/catch_all.hpp>
#include <vector>
#include <string>
#include "../SmartStrategy.h"

// Note: These tests assume the presence of 'simple.fasta' and 'simple.fastq'
// in the '../test_data/' directory, consistent with other test files.

TEST_CASE("SmartStrategy File Loading", "[strategy]") {

    SECTION("Can load a simple FASTA file") {
        TracEon::SmartStrategy strategy;
        std::string test_file_path = "../test_data/simple.fasta";

        strategy.loadFile(test_file_path);

        REQUIRE(strategy.getFileCacheSize() == 2);
        REQUIRE(strategy.getSequence("seq1") == "GATTACA");
        REQUIRE(strategy.getSequence("seq2") == "CGCGCGCGCGCGCGCGCGCGCGCGCGCG");
        REQUIRE(strategy.getQuality("seq1").empty()); // FASTA has no quality scores
    }

    SECTION("Can load a simple FASTQ file") {
        TracEon::SmartStrategy strategy;
        std::string test_file_path = "../test_data/simple.fastq";

        strategy.loadFile(test_file_path);

        REQUIRE(strategy.getFileCacheSize() == 2);

        // Check first record
        REQUIRE(strategy.getSequence("seq1") == "GATTACA");
        REQUIRE(strategy.getQuality("seq1") == "@@@DDDD");

        // Check second record
        REQUIRE(strategy.getSequence("seq2") == "CGCGCGCGCGCGCGCGCGCGCGCGCGCG");
        REQUIRE(strategy.getQuality("seq2") == "HHHIIIIIIIIIIIIIIIIIIIIIII");
    }
}

