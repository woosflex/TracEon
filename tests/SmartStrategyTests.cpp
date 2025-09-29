#include <catch2/catch_all.hpp>
#include "../SmartStrategy.h" // This file doesn't exist yet

TEST_CASE("SmartStrategy automatic encoding", "[strategy][smart]") {
    SmartStrategy strategy;

    SECTION("Correctly encodes and decodes DNA") {
        std::string dna = "GATTACAGATTACAGATTACAGATTACA"; // 28 bytes - longer sequence
        auto encoded = strategy.encode(dna);
        // Now we'll see actual compression benefits
        REQUIRE(encoded.size() < dna.size());
        REQUIRE(strategy.decode(encoded) == dna);
    }

    SECTION("Correctly encodes and decodes RNA") {
        std::string rna = "GAUUACAGAUUACAGAUUACAGAUUACA"; // 28 bytes - longer sequence
        auto encoded = strategy.encode(rna);
        REQUIRE(encoded.size() < rna.size());
        // Check that U is treated like T for compression
        REQUIRE(strategy.decode(encoded) == "GATTACAGATTACAGATTACAGATTACA");
    }

    SECTION("Correctly encodes and decodes Quality Scores with RLE") {
        std::string quality = "FFFFHHHHIIIIJJJJ"; // 16 bytes
        // Give a hint to the encoder that this is a quality string
        auto encoded = strategy.encode(quality, DataTypeHint::QualityScore);
        REQUIRE(encoded.size() < quality.size());
        REQUIRE(encoded.size() == 9); // 1 byte header + 8 bytes of RLE data
        REQUIRE(strategy.decode(encoded) == quality);
    }

    SECTION("Correctly handles Protein sequences (no compression)") {
        std::string protein = "LVFP"; // A non-compressible string
        auto encoded = strategy.encode(protein);
        // Expect 1-byte header + 4 bytes of plain text data
        REQUIRE(encoded.size() == protein.size() + 1);
        REQUIRE(strategy.decode(encoded) == protein);
    }
}