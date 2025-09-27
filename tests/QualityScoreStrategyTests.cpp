#include <catch2/catch_all.hpp>
#include "../QualityScoreRleStrategy.h" // This file doesn't exist yet

TEST_CASE("QualityScoreRleStrategy", "[strategy][quality]") {

    QualityScoreRleStrategy strategy;
    // A typical quality string with lots of repetition.
    std::string original_quality = "FFFFHHHHIIIIJJJJ"; // 16 bytes long

    SECTION("Encoding reduces size") {
        auto encoded_data = strategy.encode(original_quality);
        // We expect RLE to make this much smaller.
        REQUIRE(encoded_data.size() < original_quality.size());
    }

    SECTION("Encode and Decode returns original string") {
        auto encoded_data = strategy.encode(original_quality);
        std::string decoded_quality = strategy.decode(encoded_data);
        REQUIRE(decoded_quality == original_quality);
    }
}