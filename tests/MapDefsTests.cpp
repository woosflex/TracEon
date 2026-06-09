#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_all.hpp>
#include "MapDefs.h"

#include <string>
#include <string_view>
#include <vector>
#include <set>
#include <limits>

using namespace TracEon;

// ── HashMap Basic Operations ────────────────────────────────────────────────

TEST_CASE("HashMap<int, int> basic operations", "[MapDefs][HashMap]") {
    HashMap<int, int> map;

    SECTION("starts empty") {
        REQUIRE(map.empty());
        REQUIRE(map.size() == 0);
    }

    SECTION("insert and find") {
        map.insert({1, 100});
        map.insert({2, 200});
        map.insert({3, 300});

        REQUIRE(map.size() == 3);
        REQUIRE_FALSE(map.empty());

        auto it = map.find(1);
        REQUIRE(it != map.end());
        REQUIRE(it->second == 100);

        it = map.find(2);
        REQUIRE(it != map.end());
        REQUIRE(it->second == 200);

        it = map.find(4);
        REQUIRE(it == map.end());
    }

    SECTION("operator[] access") {
        map[1] = 10;
        map[2] = 20;
        map[3] = 30;

        REQUIRE(map[1] == 10);
        REQUIRE(map[2] == 20);
        REQUIRE(map[3] == 30);
        REQUIRE(map.size() == 3);
    }

    SECTION("erase removes element") {
        map[1] = 10;
        map[2] = 20;
        REQUIRE(map.size() == 2);

        map.erase(1);
        REQUIRE(map.size() == 1);
        REQUIRE(map.find(1) == map.end());
        REQUIRE(map.find(2) != map.end());
    }

    SECTION("clear removes all elements") {
        map[1] = 10;
        map[2] = 20;
        map[3] = 30;
        REQUIRE_FALSE(map.empty());

        map.clear();
        REQUIRE(map.empty());
        REQUIRE(map.size() == 0);
    }

    SECTION("contains returns correct results") {
        map[42] = 1;
        REQUIRE(map.contains(42));
        REQUIRE_FALSE(map.contains(99));
    }

    SECTION("count is 0 or 1 for unique keys") {
        map[7] = 1;
        REQUIRE(map.count(7) == 1);
        REQUIRE(map.count(8) == 0);
    }
}

TEST_CASE("HashMap<std::string, std::string> operations", "[MapDefs][HashMap]") {
    HashMap<std::string, std::string> map;

    SECTION("insert and retrieve strings") {
        map.insert({"key1", "value1"});
        map.insert({"key2", "value2"});
        map.insert({"key3", "value3"});

        REQUIRE(map.find("key1")->second == "value1");
        REQUIRE(map.find("key2")->second == "value2");
        REQUIRE(map.find("key3")->second == "value3");
    }

    SECTION("update existing key") {
        map.insert({"key1", "original"});
        REQUIRE(map.find("key1")->second == "original");

        map["key1"] = "updated";
        REQUIRE(map["key1"] == "updated");
        REQUIRE(map.size() == 1);
    }

    SECTION("erase nonexistent key is safe") {
        map.insert({"key1", "value1"});
        REQUIRE_NOTHROW(map.erase("nonexistent"));
        REQUIRE(map.size() == 1);
    }

    SECTION("iterate over elements") {
        map.insert({"a", "1"});
        map.insert({"b", "2"});
        map.insert({"c", "3"});

        int count = 0;
        for (const auto& [key, value] : map) {
            REQUIRE_FALSE(key.empty());
            REQUIRE_FALSE(value.empty());
            ++count;
        }
        REQUIRE(count == 3);
    }
}

// ── StringHashMap Transparent Lookup ────────────────────────────────────────

TEST_CASE("StringHashMap transparent lookup with string_view", "[MapDefs][StringHashMap][transparent]") {
    StringHashMap<int> map;

    // Insert with std::string keys
    map[std::string("sequence_chr1")]  = 1;
    map[std::string("sequence_chr2")]  = 2;
    map[std::string("sequence_chr3")]  = 3;
    map[std::string("sequence_chr10")] = 10;
    map[std::string("sequence_chrX")]  = 24;
    map[std::string("sequence_chrY")]  = 25;
    map[std::string("")]               = 0;  // empty key

    SECTION("find with string_view finds the key without conversion") {
        std::string_view sv1 = "sequence_chr1";
        auto it1 = map.find(sv1);
        REQUIRE(it1 != map.end());
        REQUIRE(it1->second == 1);

        std::string_view sv3 = "sequence_chr3";
        auto it3 = map.find(sv3);
        REQUIRE(it3 != map.end());
        REQUIRE(it3->second == 3);
    }

    SECTION("find with string_view returns end for missing key") {
        std::string_view sv = "nonexistent_sequence";
        REQUIRE(map.find(sv) == map.end());
    }

    SECTION("count with string_view") {
        std::string_view sv1 = "sequence_chr1";
        REQUIRE(map.count(sv1) == 1);

        std::string_view svMissing = "missing";
        REQUIRE(map.count(svMissing) == 0);
    }

    SECTION("contains with string_view") {
        std::string_view sv = "sequence_chrX";
        REQUIRE(map.contains(sv));

        std::string_view svBad = "chrM";
        REQUIRE_FALSE(map.contains(svBad));
    }

    SECTION("at() with string_view") {
        std::string_view sv = "sequence_chrY";
        REQUIRE(map.at(sv) == 25);
    }

    SECTION("transparent lookup finds empty key") {
        std::string_view empty = "";
        auto it = map.find(empty);
        REQUIRE(it != map.end());
        REQUIRE(it->second == 0);
    }

    SECTION("try_emplace with string_view key") {
        std::string_view sv = "sequence_chrM";
        auto [it, inserted] = map.try_emplace(sv, 26);
        REQUIRE(inserted);
        REQUIRE(it->second == 26);

        // Second try_emplace should not insert
        auto [it2, inserted2] = map.try_emplace(sv, 99);
        REQUIRE_FALSE(inserted2);
        REQUIRE(it2->second == 26);
    }

    SECTION("insert_or_assign with string_view key") {
        std::string_view sv = "sequence_chr1";  // already present
        auto [it, inserted] = map.insert_or_assign(sv, 100);
        REQUIRE_FALSE(inserted);  // was already present
        REQUIRE(it->second == 100);  // value was assigned

        // Verify via std::string lookup too
        REQUIRE(map[std::string("sequence_chr1")] == 100);
    }
}

// ── Heterogeneous Lookup with std::equal_to<> ──────────────────────────────

TEST_CASE("StringHashMap heterogeneous comparison", "[MapDefs][StringHashMap][heterogeneous]") {
    // Verify that std::equal_to<> enables comparison between string and string_view
    // without constructing a temporary std::string from the view.

    StringHashMap<int> map;
    map[std::string("chr1")] = 1;
    map[std::string("chr2")] = 2;

    SECTION("compare via const char*") {
        const char* key = "chr1";
        REQUIRE(map.contains(std::string_view(key)));
        REQUIRE(map.count(std::string_view(key)) == 1);
        REQUIRE(map.find(std::string_view(key))->second == 1);
    }

    SECTION("compare via mixed types") {
        std::string str_key = "chr2";
        std::string_view sv_key = "chr2";

        // Both should find the same entry
        REQUIRE(map.find(str_key)->second == 2);
        REQUIRE(map.find(sv_key)->second == 2);
    }
}

// ── TransparentStringHash correctness ───────────────────────────────────────

TEST_CASE("TransparentStringHash produces consistent hashes", "[MapDefs][hash]") {
    TransparentStringHash hasher;

    SECTION("same string produces same hash") {
        std::string s = "sequence_chr1";
        std::string_view sv = s;

        size_t h1 = hasher(s);
        size_t h2 = hasher(sv);
        REQUIRE(h1 == h2);
    }

    SECTION("different strings produce different hashes (with high probability)") {
        size_t h1 = hasher(std::string_view("chr1"));
        size_t h2 = hasher(std::string_view("chr2"));
        REQUIRE(h1 != h2);
    }

    SECTION("empty string hash is stable") {
        size_t h1 = hasher(std::string_view(""));
        size_t h2 = hasher(std::string(""));
        REQUIRE(h1 == h2);
    }

    SECTION("hash of identical content via string and string_view matches") {
        std::string s = "ACGTACGTACGTACGT";
        std::string_view sv = s;

        size_t h_s = hasher(s);
        size_t h_sv = hasher(sv);
        REQUIRE(h_s == h_sv);
    }

    SECTION("hash of substrings differs") {
        std::string s = "chr1_chr2";
        std::string_view sv1 = std::string_view(s.data(), 4);  // "chr1"
        std::string_view sv2 = std::string_view(s.data() + 5, 4);  // "chr2"

        REQUIRE(hasher(sv1) != hasher(sv2));
    }
}

TEST_CASE("is_transparent trait is correctly defined", "[MapDefs][hash]") {
    // Verify the transparent trait exists so that std::unordered_map-like APIs
    // can detect heterogeneous lookup support at compile time.
    REQUIRE(std::is_same_v<TransparentStringHash::is_transparent, void>);
}

// ── Edge Cases and Large Data ───────────────────────────────────────────────

TEST_CASE("HashMap edge cases", "[MapDefs][HashMap][edge]") {
    SECTION("HashMap with uint64_t keys") {
        HashMap<uint64_t, std::string> map;
        map[0] = "zero";
        map[std::numeric_limits<uint64_t>::max()] = "max";

        REQUIRE(map[0] == "zero");
        REQUIRE(map[std::numeric_limits<uint64_t>::max()] == "max");
        REQUIRE(map.size() == 2);
    }

    SECTION("StringHashMap with empty key") {
        StringHashMap<int> map;
        map[std::string("")] = -1;
        REQUIRE(map.count(std::string_view("")) == 1);
        REQUIRE(map.at(std::string_view("")) == -1);
    }

    SECTION("StringHashMap with long keys") {
        StringHashMap<int> map;
        std::string long_key(1000, 'A');
        map[long_key] = 42;

        std::string_view sv(long_key);
        REQUIRE(map.contains(sv));
        REQUIRE(map.at(sv) == 42);
    }

    SECTION("StringHashMap with special characters in keys") {
        StringHashMap<int> map;
        map[std::string("chr1_v2!@#$%")] = 1;
        map[std::string("chr2 (copy)")]  = 2;
        map[std::string("chr3\twith\ttabs")] = 3;

        REQUIRE(map.contains(std::string_view("chr1_v2!@#$%")));
        REQUIRE(map.contains(std::string_view("chr2 (copy)")));
        REQUIRE(map.contains(std::string_view("chr3\twith\ttabs")));
    }

    SECTION("HashMap with many entries (stress insertion)") {
        HashMap<int, int> map;
        constexpr int N = 10000;
        for (int i = 0; i < N; ++i) {
            map[i] = i * 2;
        }
        REQUIRE(map.size() == N);

        // Verify random access
        for (int i = 0; i < N; i += 100) {
            REQUIRE(map.contains(i));
            REQUIRE(map[i] == i * 2);
        }
    }

    SECTION("StringHashMap with many entries (stress insertion)") {
        StringHashMap<int> map;
        constexpr int N = 10000;
        for (int i = 0; i < N; ++i) {
            map[std::to_string(i)] = i;
        }
        REQUIRE(map.size() == N);

        // Verify transparent access
        for (int i = 0; i < N; i += 100) {
            std::string key = std::to_string(i);
            std::string_view sv(key);
            REQUIRE(map.contains(sv));
            REQUIRE(map.at(sv) == i);
        }
    }
}

// ── Rehash and Memory Stability ─────────────────────────────────────────────

TEST_CASE("StringHashMap survives rehash (views remain valid)", "[MapDefs][StringHashMap][rehash]") {
    // Verify that transparent lookup works correctly after the map rehashes
    // (grows beyond its initial bucket count).

    StringHashMap<int> map;

    // Insert enough elements to trigger at least one rehash
    constexpr int N = 1000;
    for (int i = 0; i < N; ++i) {
        map[std::to_string(i)] = i;
    }

    // All entries must still be findable via string_view
    for (int i = 0; i < N; ++i) {
        std::string key = std::to_string(i);
        std::string_view sv(key);
        REQUIRE(map.contains(sv));
        REQUIRE(map.at(sv) == i);
    }
}

// ── NGSIndex (HashMap<uint64_t, SequenceView>) usage pattern ────────────────

TEST_CASE("NGSIndex (uint64_t key) operations", "[MapDefs][NGSIndex]") {
    // NGSIndex is defined as HashMap<uint64_t, SequenceView> in SmartStrategy.h
    // This verifies the ankerl map works with small struct values.

    struct TestView {
        uint64_t offset;
        uint64_t length;
    };

    HashMap<uint64_t, TestView> index;

    SECTION("insert and find by numeric key") {
        index[100] = {0, 100};
        index[200] = {100, 200};
        index[300] = {300, 300};

        REQUIRE(index.contains(100));
        REQUIRE(index.at(100).offset == 0);
        REQUIRE(index.at(100).length == 100);

        REQUIRE(index.contains(200));
        REQUIRE(index.at(200).length == 200);
    }

    SECTION("lookup for non-sequential keys") {
        for (uint64_t i = 0; i < 500; ++i) {
            index[i * 7 + 3] = {i, i};
        }

        for (uint64_t i = 0; i < 500; ++i) {
            uint64_t key = i * 7 + 3;
            REQUIRE(index.contains(key));
            REQUIRE(index.at(key).length == i);
        }
    }
}

// ── API Compatibility (duck-typing) tests ───────────────────────────────────

TEST_CASE("HashMap API compatibility with std::unordered_map interface", "[MapDefs][HashMap][compat]") {
    // These tests verify that the ankerl map exposes the subset of
    // std::unordered_map's API that TracEon actually depends on.

    HashMap<std::string, int> map;

    // insert
    map.insert({"key", 1});
    REQUIRE(map.size() == 1);

    // emplace
    map.emplace(std::string("key2"), 2);
    REQUIRE(map[std::string("key2")] == 2);

    // try_emplace
    map.try_emplace(std::string("key3"), 3);
    REQUIRE(map[std::string("key3")] == 3);

    // extract (ankerl returns std::optional<value_type>)
    auto nh = map.extract(std::string("key"));
    REQUIRE(nh.has_value());
    REQUIRE(nh->first == "key");
    REQUIRE(nh->second == 1);
    REQUIRE_FALSE(map.contains(std::string("key")));
    REQUIRE(map.size() == 2);

    // Re-insert the extracted value
    map.insert(std::move(*nh));
    REQUIRE(map.contains(std::string("key")));
    REQUIRE(map.size() == 3);

    // equal_range
    auto [lo, hi] = map.equal_range(std::string("key2"));
    REQUIRE(lo != hi);
    REQUIRE(lo->second == 2);

    // reserve / load_factor / max_load_factor
    map.reserve(100);
    REQUIRE(map.load_factor() > 0.0f);
    REQUIRE(map.load_factor() <= map.max_load_factor());

    // bucket API (if available)
    REQUIRE(map.bucket_count() > 0);
}
