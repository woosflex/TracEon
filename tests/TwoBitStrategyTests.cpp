#include <catch2/catch_all.hpp>
#include "../TwoBitDnaStrategy.h"

TEST_CASE("TwoBitDnaStrategy encoding and decoding", "[strategy]") {
    TwoBitDnaStrategy strategy;
    std::string original_dna = "GATTACA";

    SECTION("Encoding reduces size") {
        auto encoded_data = strategy.encode(original_dna);
        // Check that the encoded size is smaller than the original.
        REQUIRE(encoded_data.size() < original_dna.size());
        // Specifically, it should be 4 bytes for the length + 2 bytes for 7 bases of data.
        REQUIRE(encoded_data.size() == 6);
    }

    SECTION("Encode and Decode returns original string") {
        auto encoded_data = strategy.encode(original_dna);
        std::string decoded_dna = strategy.decode(encoded_data);
        REQUIRE(decoded_dna == original_dna);
    }
}