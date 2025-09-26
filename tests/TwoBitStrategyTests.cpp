//
// Created by Adnan Raza on 9/26/2025.
//

#include <catch2/catch_all.hpp>
#include "../TwoBitDnaStrategy.h" // This file doesn't exist yet

TEST_CASE("TwoBitDnaStrategy encoding and decoding", "[strategy]") {

    TwoBitDnaStrategy strategy;
    std::string original_dna = "GATTACA";

    SECTION("Encoding reduces size") {
        auto encoded_data = strategy.encode(original_dna);
        // 7 bases * 2 bits/base = 14 bits. This should fit in 2 bytes.
        REQUIRE(encoded_data.size() < original_dna.size());
        REQUIRE(encoded_data.size() == 2);
    }

    SECTION("Encode and Decode returns original string") {
        auto encoded_data = strategy.encode(original_dna);
        std::string decoded_dna = strategy.decode(encoded_data);
        REQUIRE(decoded_dna == original_dna);
    }
}