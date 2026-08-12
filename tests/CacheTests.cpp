#include <catch2/catch_all.hpp>
#include <fstream>
#include <filesystem>
#include "Cache.h"
#include "SmartStrategy.h"

TEST_CASE("Cache Functionality", "[cache]") {

    SECTION("Can set and get a simple key-value pair") {
        TracEon::Cache cache;
        std::string key = "seq1";
        std::string value = "GATTACA";

        cache.set(key, value);
        std::string retrieved_value = cache.get(key);

        REQUIRE(retrieved_value == value);
        REQUIRE(cache.size() == 1);
    }

    SECTION("Can load a simple FASTA file") {
        TracEon::Cache cache;
        std::string test_file_path = "simple.fasta";
        
        {
             std::ofstream out(test_file_path);
             out << ">seq1\nGATTACA\n>seq2\nCGCGCGCGCGCGCGCGCGCGCGCGCGCG\n";
             out.close();
        }

        cache.loadFile(test_file_path);

        REQUIRE(cache.size() == 2);
        REQUIRE(cache.get("seq1") == "GATTACA");
        REQUIRE(cache.get("seq2") == "CGCGCGCGCGCGCGCGCGCGCGCGCGCG");
        
        std::filesystem::remove(test_file_path);
    }

    SECTION("Can load a multi-line FASTA file") {
        TracEon::Cache cache;
        std::string test_file_path = "multiline_cache.fasta";

        {
            std::ofstream out(test_file_path);
            // Two-line wrapped sequence (classic NCBI style)
            out << ">seq1\nACGT\nACGT\n>seq2\nTTTT\n";
            out.close();
        }

        cache.loadFile(test_file_path);

        REQUIRE(cache.size() == 2);
        REQUIRE(cache.get("seq1") == "ACGTACGT");
        REQUIRE(cache.get("seq2") == "TTTT");

        std::filesystem::remove(test_file_path);
    }

    SECTION("Can save and restore the cache state") {
        std::string src_filename = "temp_persist.fasta";
        std::string bin_filename = "test_cache.bin";

        // 1. Create a source file (Architecturally correct source of truth)
        {
            std::ofstream out(src_filename);
            out << ">seq1\nGATTACA\n>seq2\nCGCGCGCG\n";
            out.close();
        }

        // 2. Load into Cache 1 (Populates SmartStrategy/Arena)
        TracEon::Cache cache1;
        cache1.loadFile(src_filename);
        REQUIRE(cache1.size() == 2);

        // 3. Save to Binary (Serializes SmartStrategy)
        cache1.save(bin_filename);

        // 4. Restore into Cache 2 (Maps Binary)
        TracEon::Cache cache2;
        cache2.restore(bin_filename);

        // 5. Verify Persistence
        REQUIRE(cache2.size() == cache1.size());
        REQUIRE(cache2.get("seq1") == "GATTACA");
        REQUIRE(cache2.get("seq2") == "CGCGCGCG");

        // Cleanup
        std::filesystem::remove(src_filename);
        std::filesystem::remove(bin_filename);
    }

    SECTION("LZ4 binary cache compression and round-trip integrity") {
        std::string src_filename = "temp_lz4_test.fasta";
        std::string bin_filename = "test_cache_lz4.bin";

        // 1. Create a source file with repetitive data (compresses well)
        {
            std::ofstream out(src_filename);
            out << ">seq1\n";
            for (int i = 0; i < 100; ++i) out << "ACGTACGTACGTACGT";
            out << "\n>seq2\n";
            for (int i = 0; i < 100; ++i) out << "GCTAGCTAGCTAGCTA";
            out << "\n";
            out.close();
        }

        // 2. Load and save to v2 binary format (LZ4-compressed)
        TracEon::Cache cache1;
        cache1.loadFile(src_filename);
        REQUIRE(cache1.size() == 2);

        cache1.save(bin_filename);
        size_t compressed_size = std::filesystem::file_size(bin_filename);

        // 3. Restore from v2 binary
        TracEon::Cache cache2;
        cache2.restore(bin_filename);

        // 4. Verify all sequences match exactly
        REQUIRE(cache2.size() == 2);
        REQUIRE(cache2.get("seq1").length() == 1600);
        REQUIRE(cache2.get("seq2").length() == 1600);

        // Verify content
        std::string seq1_retrieved = cache2.get("seq1");
        std::string seq1_expected(1600, '\0');
        for (int i = 0; i < 100; ++i) {
            seq1_expected.replace(i * 16, 16, "ACGTACGTACGTACGT");
        }
        REQUIRE(seq1_retrieved == seq1_expected);

        // 5. Verify file size is reasonable (compressed < uncompressed)
        std::ifstream src_check(src_filename);
        src_check.seekg(0, std::ios::end);
        size_t src_size = src_check.tellg();
        src_check.close();

        // Compressed should be significantly smaller than source
        REQUIRE(compressed_size < src_size);

        // Cleanup
        std::filesystem::remove(src_filename);
        std::filesystem::remove(bin_filename);
    }
}

TEST_CASE("Cache::getView returns zero-copy view", "[cache]") {
    std::string fasta_path = "tmp_cache_getview.fasta";
    {
        std::ofstream out(fasta_path);
        out << ">alpha\nACGTTGCA\n>beta\nTTTTAAAA\n";
    }

    TracEon::Cache cache;
    cache.loadFile(fasta_path);

    std::string_view v = cache.getView("alpha");
    REQUIRE(v == "ACGTTGCA");
    REQUIRE(cache.getView("beta") == "TTTTAAAA");
    REQUIRE(cache.getView("nonexistent").empty());

    std::filesystem::remove(fasta_path);
}

TEST_CASE("Cache::set() entries are persisted by save()", "[cache]") {
    std::string bin_path = "tmp_set_persist.bin";

    TracEon::Cache cache1;
    cache1.set("manual_key", "ACGT");
    cache1.set("another_key", "GGGG");
    REQUIRE(cache1.size() == 2);

    cache1.save(bin_path);

    TracEon::Cache cache2;
    cache2.restore(bin_path);

    REQUIRE(cache2.size() == 2);
    REQUIRE(cache2.get("manual_key") == "ACGT");
    REQUIRE(cache2.get("another_key") == "GGGG");

    std::filesystem::remove(bin_path);
}

TEST_CASE("Cache IndexMode selection", "[cache]") {
    SECTION("Default mode is GENOME") {
        TracEon::Cache cache;
        REQUIRE(cache.getIndexMode() == TracEon::IndexMode::GENOME);
    }

    SECTION("NGS mode is retained") {
        TracEon::Cache cache(TracEon::IndexMode::NGS);
        REQUIRE(cache.getIndexMode() == TracEon::IndexMode::NGS);

        std::string fasta_path = "tmp_ngs_cache.fasta";
        {
            std::ofstream out(fasta_path);
            out << ">read1\nACGT\n>read2\nTGCA\n";
        }
        cache.loadFile(fasta_path);
        REQUIRE(cache.size() == 2);
        REQUIRE(cache.get("read1") == "ACGT");
        REQUIRE(cache.get("read2") == "TGCA");
        REQUIRE(cache.getIndexMode() == TracEon::IndexMode::NGS);

        std::filesystem::remove(fasta_path);
    }
}
// ─────────────────────────────────────────────────────────────────────────────
// Cache facade completeness — the README documents cache.loadGzipFile(...)
// and the lifecycle contract names clearCache() as a reload path, but the
// public Cache facade previously did not expose them (only SmartStrategy did),
// so the documented examples did not compile. Regression coverage for the
// facade-completion PR: loadGzipFile, clearCache, getAllKeys, getQuality,
// getDetectedFormat.
// ─────────────────────────────────────────────────────────────────────────────
#include <zlib.h>

TEST_CASE("Cache facade: README-documented loadGzipFile", "[cache][gzip]") {
    const std::string gz_path = "facade_loadgzip.fastq.gz";

    SECTION("loadGzipFile loads a .gz FASTQ with quality intact") {
        const std::string record =
            "@read1\nACGTACGT\n+\nIIIIIIII\n"
            "@read2\nTTTT\n+\nJJJJ\n";
        {
            gzFile f = gzopen(gz_path.c_str(), "wb");
            REQUIRE(f != nullptr);
            int n = gzwrite(f, record.data(), static_cast<unsigned>(record.size()));
            REQUIRE(n == static_cast<int>(record.size()));
            gzclose(f);
        }

        TracEon::Cache cache;
        cache.loadGzipFile(gz_path);

        REQUIRE(cache.size() == 2);
        REQUIRE(cache.getView("read1") == "ACGTACGT");
        REQUIRE(cache.getView("read2") == "TTTT");
        REQUIRE(cache.getQuality("read1") == "IIIIIIII");
        REQUIRE(cache.getQuality("read2") == "JJJJ");

        auto rec = cache.getFastqRecord("read1");
        REQUIRE(rec.has_value());
        REQUIRE(rec->sequence == "ACGTACGT");
        REQUIRE(rec->quality == "IIIIIIII");
    }

    SECTION("loadGzipFile rejects a non-gzip file") {
        const std::string plain_path = "facade_notgzip.txt";
        {
            std::ofstream out(plain_path);
            out << ">seq1\nACGT\n";
        }
        TracEon::Cache cache;
        REQUIRE_THROWS_AS(cache.loadGzipFile(plain_path), std::runtime_error);
        // Failure atomicity: a failed gzip load leaves the cache empty.
        REQUIRE(cache.size() == 0);
        std::filesystem::remove(plain_path);
    }

    std::filesystem::remove(gz_path);
}

TEST_CASE("Cache facade: clearCache / getAllKeys / getQuality / getDetectedFormat", "[cache]") {
    SECTION("clearCache resets state and allows reuse") {
        TracEon::Cache cache;
        cache.set("k1", "ACGT");
        cache.set("k2", "TGCA");
        REQUIRE(cache.size() == 2);

        cache.clearCache();
        REQUIRE(cache.size() == 0);
        REQUIRE(cache.hasSequence("k1") == false);

        // After clearCache() the cache is reusable (build phase again).
        cache.set("k3", "GGGG");
        REQUIRE(cache.size() == 1);
        REQUIRE(cache.get("k3") == "GGGG");

        // And a full load→clear→load cycle works (lifecycle contract).
        const std::string fasta_path = "facade_clear.fasta";
        {
            std::ofstream out(fasta_path);
            out << ">seq1\nGATTACA\n";
        }
        cache.loadFile(fasta_path);
        REQUIRE(cache.size() == 1);
        REQUIRE(cache.getView("seq1") == "GATTACA");
        cache.clearCache();
        REQUIRE(cache.size() == 0);
        cache.loadFile(fasta_path);
        REQUIRE(cache.size() == 1);
        REQUIRE(cache.get("seq1") == "GATTACA");

        std::filesystem::remove(fasta_path);
    }

    SECTION("getAllKeys enumerates all keys") {
        TracEon::Cache cache;
        cache.set("alpha", "AAAA");
        cache.set("beta", "CCCC");
        cache.set("gamma", "GGGG");
        auto keys = cache.getAllKeys();
        REQUIRE(keys.size() == 3);
        REQUIRE(std::find(keys.begin(), keys.end(), "alpha") != keys.end());
        REQUIRE(std::find(keys.begin(), keys.end(), "beta") != keys.end());
        REQUIRE(std::find(keys.begin(), keys.end(), "gamma") != keys.end());
    }

    SECTION("getQuality returns empty for FASTA and stored quality for FASTQ") {
        const std::string fasta_path = "facade_qual.fasta";
        const std::string fastq_path = "facade_qual.fastq";
        {
            std::ofstream out(fasta_path);
            out << ">chr1\nACGTACGT\n";
            out.close();
            std::ofstream out2(fastq_path);
            out2 << "@r1\nACGT\n+\nIIII\n";
        }

        TracEon::Cache fasta_cache;
        fasta_cache.loadFile(fasta_path);
        REQUIRE(fasta_cache.getQuality("chr1").empty());
        REQUIRE(fasta_cache.getDetectedFormat() == TracEon::FileFormat::DNA_FASTA);

        TracEon::Cache fastq_cache;
        fastq_cache.loadFile(fastq_path);
        REQUIRE(fastq_cache.getQuality("r1") == "IIII");
        REQUIRE(fastq_cache.getDetectedFormat() == TracEon::FileFormat::DNA_FASTQ);

        std::filesystem::remove(fasta_path);
        std::filesystem::remove(fastq_path);
    }

    SECTION("getDetectedFormat is UNKNOWN on an empty cache") {
        TracEon::Cache cache;
        REQUIRE(cache.getDetectedFormat() == TracEon::FileFormat::UNKNOWN);
    }
}
