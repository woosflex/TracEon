#include <catch2/catch_all.hpp>
#include "../SmartStrategy.h" // This file doesn't exist yet

TEST_CASE("SmartStrategy automatic encoding", "[strategy][smart]") {
    SmartStrategy strategy;

    SECTION("Correctly encodes and decodes DNA") {
        std::string dna = "GATTACA";
        auto encoded = strategy.encode(dna);
        // We expect a 1-byte header + compressed data
        REQUIRE(encoded.size() < dna.size());
        REQUIRE(strategy.decode(encoded) == dna);
    }

    SECTION("Correctly encodes and decodes RNA") {
        std::string rna = "GAUUACA";
        auto encoded = strategy.encode(rna);
        REQUIRE(encoded.size() < rna.size());
        // Check that U is treated like T for compression
        REQUIRE(strategy.decode(encoded) == "GATTACA");
    }

    SECTION("Correctly encodes and decodes Quality Scores with RLE") {
        std::string quality = "FFFFHHHHIIIIJJJJ"; // 16 bytes
        auto encoded = strategy.encode(quality, true); // The 'true' flag hints that this is a quality string
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

    SECTION("Correctly handles DNA with ambiguous 'N' bases") {
        std::string dna_with_n = "GATNACAN";
        auto encoded = strategy.encode(dna_with_n);
        REQUIRE(encoded.size() < dna_with_n.size());
        REQUIRE(strategy.decode(encoded) == dna_with_n);
    }
}