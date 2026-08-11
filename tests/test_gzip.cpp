#include "TestHelpers.h"

// test_gzip.cpp — GZIP domain (truncation, trailing garbage,
// Z_BUF_ERROR, parallel member decode, zip-bomb guard).
// Split from the monolithic SmartStrategyTests.cpp.


TEST_CASE("GZIP Decompression Correctness", "[strategy][gzip]") {
    std::string original = ">seq1\nACGT\n>seq2\nTGCA\n";
    std::string gz_filename = "temp_test.fasta.gz";

    gzFile file = gzopen(gz_filename.c_str(), "wb");
    REQUIRE(file != nullptr);
    int written = gzwrite(file, original.data(), (unsigned int)original.size());
    REQUIRE(written == (int)original.size());
    gzclose(file);

    TracEon::SmartStrategy strategy;
    
    SECTION("Can explicit load .gz file") {
        strategy.loadGzipFile(gz_filename);
        REQUIRE(strategy.getFileCacheSize() == 2);
        REQUIRE(strategy.getSequence("seq1") == "ACGT");
    }

    SECTION("Can auto-detect .gz extension") {
        strategy.loadFile(gz_filename);
        REQUIRE(strategy.getFileCacheSize() == 2);
        REQUIRE(strategy.getSequence("seq1") == "ACGT");
    }
    
    if (fs::exists(gz_filename)) fs::remove(gz_filename);
}


TEST_CASE("Parallel GZIP decompression — concatenated streams", "[strategy][gzip][parallel]") {
    std::string part1 = ">seq1\nACGTACGT\n>seq2\nTGCATGCA\n";
    std::string part2 = ">seq3\nAAAAAAAA\n>seq4\nCCCCCCCC\n";

    std::string tmp1 = "tmp_stream1.gz";
    std::string tmp2 = "tmp_stream2.gz";
    std::string concat_path = "tmp_concat.fasta.gz";

    auto bytes1 = compress_to_gzip_bytes(part1, tmp1);
    auto bytes2 = compress_to_gzip_bytes(part2, tmp2);

    {
        std::ofstream out(concat_path, std::ios::binary);
        out.write(bytes1.data(), bytes1.size());
        out.write(bytes2.data(), bytes2.size());
    }

    SECTION("All sequences present after concatenated load") {
        TracEon::SmartStrategy strategy;
        strategy.loadFile(concat_path);

        REQUIRE(strategy.getFileCacheSize() == 4);
        REQUIRE(strategy.getSequence("seq1") == "ACGTACGT");
        REQUIRE(strategy.getSequence("seq2") == "TGCATGCA");
        REQUIRE(strategy.getSequence("seq3") == "AAAAAAAA");
        REQUIRE(strategy.getSequence("seq4") == "CCCCCCCC");
    }

    fs::remove(tmp1);
    fs::remove(tmp2);
    fs::remove(concat_path);
}


TEST_CASE("Parallel GZIP — single-stream fallback unchanged", "[strategy][gzip][parallel]") {
    std::string data = ">chr1\nACGTACGTACGT\n>chr2\nTTTTGGGGCCCC\n";
    std::string gz_path = "tmp_single_stream.fasta.gz";

    gzFile f = gzopen(gz_path.c_str(), "wb");
    gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
    gzclose(f);

    TracEon::SmartStrategy strategy;
    strategy.loadFile(gz_path);

    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("chr1") == "ACGTACGTACGT");
    REQUIRE(strategy.getSequence("chr2") == "TTTTGGGGCCCC");

    fs::remove(gz_path);
}


TEST_CASE("Parallel GZIP — two large streams actually exercise loadGzipParallel", "[strategy][gzip][parallel]") {
    // Each compressed stream must exceed PARALLEL_GZIP_THRESHOLD/2 (512 KB) so
    // the concatenated file is > 1 MB, triggering loadGzipParallel.
    // FASTA nucleotide sequences compress very well with zlib (50-200:1 for
    // periodic patterns). Solution: use gzip level 0 (store mode), which wraps
    // data in a valid GZIP container without compression. This keeps the file
    // size close to the uncompressed size, guaranteeing the threshold is crossed.

    const int SEQ_LEN = 800;
    const char BASES[] = "ACGT";

    auto make_fasta = [&](const std::string& prefix, int count) {
        std::string s;
        for (int i = 0; i < count; ++i) {
            s += '>';
            s += prefix;
            s += std::to_string(i);
            s += '\n';
            for (int j = 0; j < SEQ_LEN; ++j)
                s += BASES[(i + j) % 4];
            s += '\n';
        }
        return s;
    };

    // Each stream: 1500 records × ~800 bytes ≈ 1.2 MB.
    // At gzip level 0 (no compression), file size ≈ uncompressed ≈ 1.2 MB per stream.
    // Total: ~2.4 MB → well above PARALLEL_GZIP_THRESHOLD (1 MB).
    std::string part1 = make_fasta("stream1_seq", 1500);
    std::string part2 = make_fasta("stream2_seq", 1500);

    std::string tmp1 = "tmp_large_stream1.gz";
    std::string tmp2 = "tmp_large_stream2.gz";
    std::string concat_path = "tmp_large_concat.fasta.gz";

    // "wb0" = gzip store mode (level 0): valid GZIP format, no compression.
    auto compress_stream = [](const std::string& data, const std::string& path) {
        gzFile f = gzopen(path.c_str(), "wb0");
        REQUIRE(f != nullptr);
        gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
        gzclose(f);
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        size_t sz = static_cast<size_t>(in.tellg());
        in.seekg(0);
        std::vector<char> bytes(sz);
        in.read(bytes.data(), sz);
        return bytes;
    };

    auto bytes1 = compress_stream(part1, tmp1);
    auto bytes2 = compress_stream(part2, tmp2);

    {
        std::ofstream out(concat_path, std::ios::binary);
        out.write(bytes1.data(), bytes1.size());
        out.write(bytes2.data(), bytes2.size());
    }

    // Guard: concatenated file must exceed threshold so this test is meaningful.
    size_t concat_size = static_cast<size_t>(
        std::filesystem::file_size(concat_path));
    REQUIRE(concat_size > 1024 * 1024); // must trigger loadGzipParallel

    TracEon::SmartStrategy strategy;
    strategy.loadFile(concat_path);

    REQUIRE(strategy.getFileCacheSize() == 3000);
    REQUIRE(strategy.hasSequence("stream1_seq0"));
    REQUIRE(static_cast<int>(strategy.getSequence("stream1_seq0").size()) == SEQ_LEN);
    REQUIRE(strategy.hasSequence("stream1_seq1499"));
    REQUIRE(strategy.hasSequence("stream2_seq0"));
    REQUIRE(static_cast<int>(strategy.getSequence("stream2_seq0").size()) == SEQ_LEN);
    REQUIRE(strategy.hasSequence("stream2_seq1499"));

    fs::remove(tmp1);
    fs::remove(tmp2);
    fs::remove(concat_path);
}


TEST_CASE("Parallel GZIP — coincidental magic bytes in a single-member stream don't cause a false split",
          "[strategy][gzip][parallel]") {
    // Regression test: scanGzipStreams() used to look for the raw byte
    // sequence 0x1f 0x8b 0x08 anywhere in the file to find concatenated
    // stream boundaries. That sequence can appear by chance *inside* a
    // single member's compressed payload (observed on real NCBI FASTA.gz
    // files). Splitting there feeds inflate() a misaligned/truncated
    // deflate stream, which can explode into gigabytes of garbage via LZ77
    // back-references — tripping the OOM guard on a perfectly valid file.
    //
    // Level 0 (store mode) gzip embeds input bytes close to verbatim, so
    // planting 0x1f 0x8b 0x08 in the plaintext reliably reproduces a
    // coincidental match in the compressed output.
    const int SEQ_LEN = 800;
    const char BASES[] = "ACGT";

    auto make_fasta = [&](const std::string& prefix, int count) {
        std::string s;
        for (int i = 0; i < count; ++i) {
            s += '>';
            s += prefix;
            s += std::to_string(i);
            s += '\n';
            for (int j = 0; j < SEQ_LEN; ++j)
                s += BASES[(i + j) % 4];
            s += '\n';
        }
        return s;
    };

    std::string content = make_fasta("bogus_magic_seq", 1500); // ~1.2 MB, single member

    // Plant a coincidental gzip-header-looking byte triple in the middle of
    // a sequence line (not on a newline boundary, so it can't be mistaken
    // for record structure).
    size_t plant_at = content.size() / 2;
    content[plant_at]     = static_cast<char>(0x1f);
    content[plant_at + 1] = static_cast<char>(0x8b);
    content[plant_at + 2] = static_cast<char>(0x08);

    std::string gz_path = "tmp_false_positive_magic.fasta.gz";
    gzFile f = gzopen(gz_path.c_str(), "wb0"); // store mode: preserves planted bytes verbatim
    REQUIRE(f != nullptr);
    gzwrite(f, content.data(), static_cast<unsigned>(content.size()));
    gzclose(f);

    size_t gz_size = static_cast<size_t>(std::filesystem::file_size(gz_path));
    REQUIRE(gz_size > 1024 * 1024); // must trigger loadGzipParallel's dispatch path

    TracEon::SmartStrategy strategy;
    REQUIRE_NOTHROW(strategy.loadFile(gz_path));

    REQUIRE(strategy.getFileCacheSize() == 1500);
    REQUIRE(strategy.hasSequence("bogus_magic_seq0"));
    REQUIRE(static_cast<int>(strategy.getSequence("bogus_magic_seq0").size()) == SEQ_LEN);
    REQUIRE(strategy.hasSequence("bogus_magic_seq1499"));
    REQUIRE(static_cast<int>(strategy.getSequence("bogus_magic_seq1499").size()) == SEQ_LEN);

    fs::remove(gz_path);
}


// ═══════════════════════════════════════════════════════════════════════════
// Bug C regression: truncated or garbage-tailed gzip must be rejected
// ═══════════════════════════════════════════════════════════════════════════
// loadGzipSingleStream() only threw on gzread() < 0 and treated bytes_read
// == 0 as clean EOF. zlib-ng's zlib-compat gzread returns 0 (not -1) for a
// stream truncated mid-deflate, so partial data was silently served as
// complete (gzerror == Z_BUF_ERROR 'unexpected end of file'). And gz_look()
// silently discards trailing garbage after the last member with err == Z_OK,
// so a gzerror check alone cannot catch 'valid.gz || garbage || valid2.gz'.
// The fix: consult gzerror at bytes_read == 0, reconcile gzoffset() against
// the file size, and have scanGzipStreams() reject leftover bytes after the
// final member (raw inflate knows the exact member end).

TEST_CASE("GZIP — valid single-member file still loads (control)",
          "[strategy][gzip][bugC]") {
    // Control: a well-formed file must load untouched.
    std::string data = ">chr1\nACGTACGT\n>chr2\nTTTTGGGG\n";
    std::string gz = "bugC_control.fasta.gz";
    gzFile f = gzopen(gz.c_str(), "wb");
    REQUIRE(f != nullptr);
    gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
    gzclose(f);

    TracEon::SmartStrategy strategy;
    REQUIRE_NOTHROW(strategy.loadFile(gz));
    REQUIRE(strategy.getFileCacheSize() == 2);
    REQUIRE(strategy.getSequence("chr1") == "ACGTACGT");
    REQUIRE(strategy.getSequence("chr2") == "TTTTGGGG");
    fs::remove(gz);
}


TEST_CASE("GZIP — truncated single-member stream rejected (no partial data served)",
          "[strategy][gzip][bugC][errors]") {
    std::string data = ">chr1\nACGTACGT\n>chr2\nTTTTGGGG\n";
    std::string gz = "bugC_trunc_single.gz";
    {
        gzFile f = gzopen(gz.c_str(), "wb");
        REQUIRE(f != nullptr);
        gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
        gzclose(f);
    }
    // Cut the file in half — the deflate stream can no longer reach
    // Z_STREAM_END. gzread returns 0 with gzerror Z_BUF_ERROR; the old code
    // accepted the partial payload as complete.
    std::string full;
    {
        std::ifstream in(gz, std::ios::binary);
        full.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    std::ofstream out(gz, std::ios::binary);
    out.write(full.data(), static_cast<std::streamsize>(full.size() / 2));
    out.close();

    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadFile(gz), std::runtime_error);
    // No partial data may be served as the loaded cache.
    REQUIRE(strategy.getFileCacheSize() == 0);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::UNKNOWN);
    fs::remove(gz);
}


TEST_CASE("GZIP — trailing garbage rejected via gzoffset (<1MB file, large tail)",
          "[strategy][gzip][bugC][errors]") {
    // A small valid member followed by a tail larger than zlib-ng's internal
    // 128 KiB read buffer: gzoffset() then reports a position short of the
    // file size and the file is rejected. (Tails ≤ 128 KiB are absorbed by
    // that buffer during the final read and are indistinguishable from a
    // clean EOF through the gzread() API — a documented limitation; files
    // ≥ 1 MiB get exact coverage from scanGzipStreams().)
    std::string data = ">chr1\nACGTACGT\n>chr2\nTTTTGGGG\n";
    std::string gz = "bugC_tail_gzoffset.gz";
    {
        gzFile f = gzopen(gz.c_str(), "wb");
        REQUIRE(f != nullptr);
        gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
        gzclose(f);
        std::ofstream out(gz, std::ios::binary | std::ios::app);
        const std::string garbage(200 * 1024, 'X'); // 200 KiB of junk
        out.write(garbage.data(), static_cast<std::streamsize>(garbage.size()));
    }
    REQUIRE(fs::file_size(gz) < 1024u * 1024u); // must stay on the single-stream path

    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadFile(gz), std::runtime_error);
    REQUIRE(strategy.getFileCacheSize() == 0);
    fs::remove(gz);
}


TEST_CASE("GZIP — trailing garbage rejected via scan (≥1MB single member)",
          "[strategy][gzip][bugC][errors]") {
    // Level-0 (store mode) gzip keeps the compressed size ≈ raw size, so the
    // member comfortably exceeds PARALLEL_GZIP_THRESHOLD (1 MiB). scanGzip-
    // Streams() then decodes the member with raw inflate and must reject the
    // leftover garbage after the exact member end — including the case where
    // a SECOND valid member sits beyond the garbage ('valid || garbage ||
    // valid2'), which gzread() would silently truncate at the garbage.
    const size_t SEQ_LEN = 800;
    const char BASES[] = "ACGT";
    std::string content;
    for (size_t i = 0; i < 1500; ++i) {
        content += '>';
        content += "tail_seq" + std::to_string(i);
        content += '\n';
        for (size_t j = 0; j < SEQ_LEN; ++j) content += BASES[(i + j) % 4];
        content += '\n';
    }
    std::string gz = "bugC_tail_scan.gz";
    std::string tmp2member = "bugC_tail_scan_member2.gz";
    {
        gzFile f = gzopen(gz.c_str(), "wb0");
        REQUIRE(f != nullptr);
        gzwrite(f, content.data(), static_cast<unsigned>(content.size()));
        gzclose(f);
        // Append: garbage (NOT a gzip member start) followed by a SECOND
        // valid member. gzread() would silently truncate at the garbage and
        // drop member 2; the raw-inflate scan must reject the whole file.
        std::string member2 = ">chr9\nACGT\n";
        gzFile f2 = gzopen(tmp2member.c_str(), "wb");
        REQUIRE(f2 != nullptr);
        gzwrite(f2, member2.data(), static_cast<unsigned>(member2.size()));
        gzclose(f2);
        std::ifstream in2(tmp2member, std::ios::binary | std::ios::ate);
        size_t sz2 = static_cast<size_t>(in2.tellg());
        in2.seekg(0);
        std::vector<char> bytes2(sz2);
        in2.read(bytes2.data(), static_cast<std::streamsize>(sz2));

        std::ofstream out(gz, std::ios::binary | std::ios::app);
        const std::string garbage = "THIS IS NOT GZIP DATA AT ALL";
        out.write(garbage.data(), static_cast<std::streamsize>(garbage.size()));
        out.write(bytes2.data(), static_cast<std::streamsize>(bytes2.size()));
    }
    REQUIRE(fs::file_size(gz) > 1024u * 1024u);

    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadFile(gz), std::runtime_error);
    REQUIRE(strategy.getFileCacheSize() == 0);
    fs::remove(tmp2member);
    fs::remove(gz);
}


TEST_CASE("GZIP — truncated multi-member file rejected (parallel path)",
          "[strategy][gzip][bugC][errors][parallel]") {
    // Two level-0 members concatenated (> 1 MiB total → parallel path), with
    // the second member cut mid-stream. scanGzipStreams() decodes the first
    // member cleanly, then fails to reach Z_STREAM_END on the second and must
    // reject the file with an accurate truncation error (previously the
    // parallel decoder ballooned memory on Z_BUF_ERROR before failing with a
    // misleading OOM message).
    const size_t SEQ_LEN = 800;
    const char BASES[] = "ACGT";
    auto make_fasta = [&](const std::string& prefix, size_t count) {
        std::string s;
        for (size_t i = 0; i < count; ++i) {
            s += '>';
            s += prefix + std::to_string(i);
            s += '\n';
            for (size_t j = 0; j < SEQ_LEN; ++j) s += BASES[(i + j) % 4];
            s += '\n';
        }
        return s;
    };
    auto compress_member = [](const std::string& data, const std::string& path) {
        gzFile f = gzopen(path.c_str(), "wb0");
        REQUIRE(f != nullptr);
        gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
        gzclose(f);
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        size_t sz = static_cast<size_t>(in.tellg());
        in.seekg(0);
        std::vector<char> bytes(sz);
        in.read(bytes.data(), static_cast<std::streamsize>(sz));
        return bytes;
    };

    std::string part1 = make_fasta("m1_", 1500); // ~1.2 MB raw
    std::string part2 = make_fasta("m2_", 1500);
    std::string tmp1 = "bugC_trunc_tmp1.gz";
    std::string tmp2 = "bugC_trunc_tmp2.gz";
    auto bytes1 = compress_member(part1, tmp1);
    auto bytes2 = compress_member(part2, tmp2);

    std::string gz = "bugC_trunc_multi.gz";
    {
        std::ofstream out(gz, std::ios::binary);
        out.write(bytes1.data(), static_cast<std::streamsize>(bytes1.size()));
        // Write only the first 90% of member 2 → truncated mid-deflate.
        size_t cut = bytes2.size() * 9 / 10;
        out.write(bytes2.data(), static_cast<std::streamsize>(cut));
    }
    REQUIRE(fs::file_size(gz) > 1024u * 1024u);

    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadFile(gz), std::runtime_error);
    REQUIRE(strategy.getFileCacheSize() == 0);

    fs::remove(tmp1);
    fs::remove(tmp2);
    fs::remove(gz);
}


TEST_CASE("GZIP — concatenated valid + truncated tail member rejected promptly (tiny fixture)",
          "[strategy][gzip][vuln][errors]") {
    // vuln-0002: a tiny concatenated gzip (valid member + truncated tail
    // member) must throw promptly with a truncation message — the tail must
    // not be accepted as a clean EOF, and no unbounded buffer growth may
    // happen. This exercises the single-stream path (< 1 MiB).
    auto compress_bytes = [](const std::string& data, const std::string& path) {
        gzFile f = gzopen(path.c_str(), "wb");
        REQUIRE(f != nullptr);
        gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
        gzclose(f);
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        size_t sz = static_cast<size_t>(in.tellg());
        in.seekg(0);
        std::vector<char> bytes(sz);
        in.read(bytes.data(), static_cast<std::streamsize>(sz));
        return bytes;
    };

    std::string tmp1 = "vuln0002_m1.gz";
    std::string tmp2 = "vuln0002_m2.gz";
    std::string content1 = ">m1\nACGTACGTACGT\n>m1b\nTTTT\n";
    // Large enough that a 50% cut lands mid-deflate, never on a member boundary.
    std::string content2 = ">m2\n" + std::string(4096, 'G') + "\n";
    auto bytes1 = compress_bytes(content1, tmp1);
    auto bytes2 = compress_bytes(content2, tmp2);

    std::string gz = "vuln0002_trunc_tail_tiny.gz";
    {
        std::ofstream out(gz, std::ios::binary);
        out.write(bytes1.data(), static_cast<std::streamsize>(bytes1.size()));
        size_t cut = bytes2.size() / 2;
        out.write(bytes2.data(), static_cast<std::streamsize>(cut));
    }
    REQUIRE(fs::file_size(gz) < 1024u * 1024u); // tiny fixture, no ballooning

    TracEon::SmartStrategy strategy;
    require_throws_with_msg([&] { strategy.loadFile(gz); }, "truncated");
    REQUIRE(strategy.getFileCacheSize() == 0);

    fs::remove(tmp1);
    fs::remove(tmp2);
    fs::remove(gz);
}


TEST_CASE("GZIP — concatenated valid + truncated tail member rejected promptly (parallel path)",
          "[strategy][gzip][vuln][errors][parallel]") {
    // vuln-0002 on the ≥ 1 MiB path: scanGzipStreams() decodes every member
    // with raw inflate and must reject the truncated tail BEFORE any parallel
    // per-stream buffers are allocated (the old code reached the parallel
    // decoder, which grew the stream buffer geometrically on Z_BUF_ERROR —
    // treating input-exhaustion as output-full — until the aggregate budget
    // guard tripped, ~3.4 GB from a 23-byte tail, with a misleading OOM
    // message). loadGzipParallel() additionally distinguishes the two cases
    // (Z_BUF_ERROR + avail_in == 0 → "GZIP stream truncated: unexpected end
    // of input") so a truncated stream can never balloon memory even if it is
    // reached with crafted/stale offsets.
    const size_t SEQ_LEN = 800;
    const char BASES[] = "ACGT";
    auto make_fasta = [&](const std::string& prefix, size_t count) {
        std::string s;
        for (size_t i = 0; i < count; ++i) {
            s += '>';
            s += prefix + std::to_string(i);
            s += '\n';
            for (size_t j = 0; j < SEQ_LEN; ++j) s += BASES[(i + j) % 4];
            s += '\n';
        }
        return s;
    };
    auto compress_member = [](const std::string& data, const std::string& path) {
        gzFile f = gzopen(path.c_str(), "wb0"); // level-0: compressed size ≈ raw size
        REQUIRE(f != nullptr);
        gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
        gzclose(f);
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        size_t sz = static_cast<size_t>(in.tellg());
        in.seekg(0);
        std::vector<char> bytes(sz);
        in.read(bytes.data(), static_cast<std::streamsize>(sz));
        return bytes;
    };

    std::string part1 = make_fasta("m1_", 1300); // ~1.06 MB raw → member > 1 MiB
    std::string part2 = make_fasta("m2_", 20);
    std::string tmp1 = "vuln0002_tmp1.gz";
    std::string tmp2 = "vuln0002_tmp2.gz";
    auto bytes1 = compress_member(part1, tmp1);
    auto bytes2 = compress_member(part2, tmp2);

    std::string gz = "vuln0002_trunc_tail_parallel.gz";
    {
        std::ofstream out(gz, std::ios::binary);
        out.write(bytes1.data(), static_cast<std::streamsize>(bytes1.size()));
        // Truncate member 2 at 50% — its deflate stream can never reach
        // Z_STREAM_END.
        size_t cut = bytes2.size() / 2;
        out.write(bytes2.data(), static_cast<std::streamsize>(cut));
    }
    REQUIRE(fs::file_size(gz) > 1024u * 1024u); // parallel path is exercised

    TracEon::SmartStrategy strategy;
    require_throws_with_msg([&] { strategy.loadFile(gz); }, "truncated");
    REQUIRE(strategy.getFileCacheSize() == 0);

    fs::remove(tmp1);
    fs::remove(tmp2);
    fs::remove(gz);
}


TEST_CASE("GZIP — highly-compressible single-stream file still grows and loads (zip-bomb guard control)",
          "[strategy][gzip][vuln]") {
    // vuln-0005: the growth loop in loadGzipSingleStream() now pre-checks every
    // geometric resize against available memory. A legitimate file that expands
    // well past 1 MiB (CHUNK_SIZE) from a tiny compressed input must still load
    // — the guard only trips on allocations above half of available memory,
    // which this 2 MiB expansion stays far below.
    const size_t SEQ_LEN = 2 * 1024 * 1024; // 2 MiB of 'A'
    std::string data = ">chr1\n" + std::string(SEQ_LEN, 'A') + "\n";
    std::string gz = "vuln0005_growth_control.gz";
    {
        gzFile f = gzopen(gz.c_str(), "wb");
        REQUIRE(f != nullptr);
        gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
        gzclose(f);
    }
    REQUIRE(fs::file_size(gz) < 1024u * 1024u); // must stay on the single-stream path

    TracEon::SmartStrategy strategy;
    REQUIRE_NOTHROW(strategy.loadFile(gz));
    REQUIRE(strategy.getFileCacheSize() == 1);
    const std::string seq = strategy.getSequence("chr1");
    REQUIRE(seq.size() == SEQ_LEN);
    REQUIRE(std::all_of(seq.begin(), seq.end(), [](char c) { return c == 'A'; }));
    fs::remove(gz);
}
