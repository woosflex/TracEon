#include "TestHelpers.h"

// test_binary_cache.cpp — `.traceon` v4 migration + v5 binary cache domain
// (format, CRC32C, mutation/truncation hardening, count
// bounds, legacy rejection, smart compression round-trips).
// Split from the monolithic SmartStrategyTests.cpp.


TEST_CASE("addEntry() round-trips through saveBinary/loadBinary",
          "[strategy][addentry][smart_compression]") {
    // Under the immutable-after-load contract (ADR-001, Bug 3 fix), addEntry()
    // is only legal before any load or after clearCache(), so a manual
    // build can no longer be mixed with a loadFile() into one object. This
    // still exercises addEntry()'s incremental serialized_size_estimate_
    // updates on top of the sizeof(uint64_t) record-count base, and the
    // loadBinary→saveBinary path (refreshPayloadEstimate seeding) is covered
    // by the Smart Compression round-trip tests below.
    const std::string bin = "addentry_roundtrip.bin";

    TracEon::SmartStrategy strategy;
    strategy.addEntry("manual1", "GATTACA", "IIIIIII");
    strategy.addEntry("manual2", "CATCATCAT", "");
    REQUIRE(strategy.getFileCacheSize() == 2);

    strategy.saveBinary(bin);

    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getFileCacheSize() == 2);
    REQUIRE(restored.getSequence("manual1") == "GATTACA");
    REQUIRE(restored.getQuality("manual1") == "IIIIIII");
    REQUIRE(restored.getSequence("manual2") == "CATCATCAT");

    fs::remove(bin);
}


// ── Smart Compression tests ───────────────────────────────────────────────────

TEST_CASE("Smart Compression — large DNA payload uses LZ4_HC (v5, better ratio)",
          "[strategy][smart_compression][lz4hc]") {
    const std::string src = "sc_large_dna.fasta";
    const std::string bin = "sc_large_dna.bin";
    {
        std::ofstream out(src);
        // 3000 × 4000-char sequences ≈ 12 MiB payload — above the 10 MiB HC threshold
        const std::string seq(4000, 'A');
        for (int i = 1; i <= 3000; ++i)
            out << ">seq" << i << '\n' << seq << '\n';
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTA);
    strategy.saveBinary(bin);

    // All save paths now write the portable streaming v5 format.
    {
        std::ifstream f(bin, std::ios::binary);
        char magic[4];
        f.read(magic, 4);
        REQUIRE(magic[3] == '\x05');
    }

    // LZ4_HC on homopolymer runs achieves extreme ratios — expect very small file
    REQUIRE(fs::file_size(bin) < 500'000); // < 500 KB for ~12 MiB of AAAA...

    // Full round-trip correctness
    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getFileCacheSize() == 3000);
    REQUIRE(restored.getSequence("seq1")    == std::string(4000, 'A'));
    REQUIRE(restored.getSequence("seq3000") == std::string(4000, 'A'));

    fs::remove(src);
    fs::remove(bin);
}


TEST_CASE("Smart Compression — large protein payload uses LZ4_default (not HC)",
          "[strategy][smart_compression][lz4default_protein]") {
    const std::string src = "sc_large_protein.fasta";
    const std::string bin = "sc_large_protein.bin";
    {
        std::ofstream out(src);
        // 3000 × 4000-char protein sequences ≈ 12 MiB payload
        // 'E' (glutamic acid) is not in the nucleotide set — forces PROTEIN_FASTA
        const std::string seq(4000, 'E');
        for (int i = 1; i <= 3000; ++i)
            out << ">prot" << i << '\n' << seq << '\n';
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::PROTEIN_FASTA);
    strategy.saveBinary(bin);

    // Portable v5 streaming LZ4 Frame format.
    {
        std::ifstream f(bin, std::ios::binary);
        char magic[4];
        f.read(magic, 4);
        REQUIRE(magic[3] == '\x05');
    }

    // Round-trip correctness
    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getFileCacheSize() == 3000);
    REQUIRE(restored.getSequence("prot1")    == std::string(4000, 'E'));
    REQUIRE(restored.getSequence("prot3000") == std::string(4000, 'E'));

    fs::remove(src);
    fs::remove(bin);
}


TEST_CASE("Smart Compression — small DNA payload uses LZ4_default (below HC threshold)",
          "[strategy][smart_compression][lz4default_small]") {
    const std::string src = "sc_small_dna.fasta";
    const std::string bin = "sc_small_dna.bin";
    {
        std::ofstream out(src);
        out << ">chr1\nACGT\n>chr2\nTGCA\n";
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src);
    REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::DNA_FASTA);
    strategy.saveBinary(bin);

    // LZ4_default also writes the portable v5 format.
    {
        std::ifstream f(bin, std::ios::binary);
        char magic[4];
        f.read(magic, 4);
        REQUIRE(magic[3] == '\x05');
    }

    // Round-trip correctness
    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getFileCacheSize() == 2);
    REQUIRE(restored.getSequence("chr1") == "ACGT");
    REQUIRE(restored.getSequence("chr2") == "TGCA");

    fs::remove(src);
    fs::remove(bin);
}


TEST_CASE("loadBinary throws on missing file", "[strategy][errors]") {
    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadBinary("nonexistent_cache_xyz.bin"), std::runtime_error);
}


TEST_CASE("loadBinary rejects corrupt/truncated file", "[strategy][errors]") {
    std::string bad_file = "corrupt.bin";
    {
        std::ofstream out(bad_file, std::ios::binary);
        // Valid v4 magic + codec flags + mode, then truncate before the
        // length fields / CRC / frame — deliberately truncated header.
        out.write("TRO\x04", 4);
        out.write("\x01", 1); // codec flags: LZ4 Frame
        out.write("\x00", 1); // mode byte (GENOME)
    }
    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadBinary(bad_file), std::runtime_error);
    fs::remove(bad_file);
}


// ── v4 round-trips ───────────────────────────────────────────────────────────

TEST_CASE("v5 round-trip GENOME mode through saveBinary/loadBinary", "[strategy][binary][v5]") {
    const std::string src = "v5_rt_genome.fasta";
    const std::string bin = "v5_rt_genome.bin";
    {
        std::ofstream out(src);
        out << ">chr1\nACGTACGT\n>chr2\nTTTTGGGG\n>chr3\nNNNN\n";
    }
    TracEon::SmartStrategy strategy;
    strategy.loadFile(src);
    strategy.saveBinary(bin);

    // Header sanity: v5 magic, LZ4 Frame codec flag, GENOME mode byte.
    {
        std::ifstream f(bin, std::ios::binary);
        char magic[4]; f.read(magic, 4);
        REQUIRE(std::memcmp(magic, "TRO\x05", 4) == 0);
        char flags; f.read(&flags, 1);
        REQUIRE(flags == '\x01');
        char mode; f.read(&mode, 1);
        REQUIRE(mode == '\x00');
    }

    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getIndexMode() == TracEon::IndexMode::GENOME);
    REQUIRE(restored.getFileCacheSize() == 3);
    REQUIRE(restored.getSequence("chr1") == "ACGTACGT");
    REQUIRE(restored.getSequence("chr2") == "TTTTGGGG");
    REQUIRE(restored.getSequence("chr3") == "NNNN");
    fs::remove(src);
    fs::remove(bin);
}


TEST_CASE("v5 round-trip NGS mode through saveBinary/loadBinary", "[strategy][binary][v5][ngs]") {
    const std::string src = "v5_rt_ngs.fastq";
    const std::string bin = "v5_rt_ngs.bin";
    {
        std::ofstream out(src);
        out << "@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\n####\n";
    }
    TracEon::SmartStrategy strategy(TracEon::IndexMode::NGS);
    strategy.loadFile(src);
    strategy.saveBinary(bin);

    {
        std::ifstream f(bin, std::ios::binary);
        char magic[4]; f.read(magic, 4);
        REQUIRE(std::memcmp(magic, "TRO\x05", 4) == 0);
        char flags; f.read(&flags, 1);
        REQUIRE(flags == '\x01');
        char mode; f.read(&mode, 1);
        REQUIRE(mode == '\x01'); // NGS mode byte
    }

    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getIndexMode() == TracEon::IndexMode::NGS);
    REQUIRE(restored.getFileCacheSize() == 2);
    REQUIRE(restored.getSequence("read1") == "ACGT");
    REQUIRE(restored.getQuality("read1") == "IIII");
    REQUIRE(restored.getSequence("read2") == "TGCA");
    REQUIRE(restored.getQuality("read2") == "####");
    fs::remove(src);
    fs::remove(bin);
}


// ── Legacy rejection (v1/v2/v3 unsupported; v4 migration reader retained) ──

TEST_CASE("loadBinary rejects legacy v1/v2/v3 magic (regenerate with v2.0.0)", "[strategy][binary][v4][reject]") {
    std::vector<char> payload;
    write_genome_payload(payload, {{"seq1", "ACGT", ""}});
    const auto frame = lz4f_compress(payload);
    for (uint8_t ver : {0x01, 0x02, 0x03}) {
        const std::string bin = "legacy_reject_" + std::to_string(ver) + ".bin";
        const auto blob = build_v4_blob(0, payload, frame, /*magic_ver=*/ver);
        write_blob(bin, blob);
        TracEon::SmartStrategy strategy;
        require_throws_with_msg([&] { strategy.loadBinary(bin); },
                                "unsupported cache version; expected v4 or v5");
        // Failure atomicity: a rejected load leaves the cache empty.
        REQUIRE(strategy.getFileCacheSize() == 0);
        fs::remove(bin);
    }
}


TEST_CASE("loadBinary rejects unknown version and non-TRO magic", "[strategy][binary][v4][reject]") {
    std::vector<char> payload;
    write_genome_payload(payload, {{"seq1", "ACGT", ""}});
    const auto frame = lz4f_compress(payload);

    // "TRO\x09" — a TRO magic with an unknown future version.
    {
        const std::string bin = "magic_ver9.bin";
        const auto blob = build_v4_blob(0, payload, frame, /*magic_ver=*/0x09);
        write_blob(bin, blob);
        TracEon::SmartStrategy strategy;
        require_throws_with_msg([&] { strategy.loadBinary(bin); },
                                "unsupported cache version; expected v4 or v5");
        fs::remove(bin);
    }
    // "XXXX" — not a TRO magic at all.
    {
        const std::string bin = "magic_xxxx.bin";
        auto blob = build_v4_blob(0, payload, frame);
        blob[0] = 'X'; blob[1] = 'X'; blob[2] = 'X'; blob[3] = 'X';
        write_blob(bin, blob);
        TracEon::SmartStrategy strategy;
        require_throws_with_msg([&] { strategy.loadBinary(bin); },
                                "Invalid binary cache magic bytes");
        fs::remove(bin);
    }
}


TEST_CASE("loadBinary rejects unsupported codec flags and invalid index mode", "[strategy][binary][v4][reject]") {
    std::vector<char> payload;
    write_genome_payload(payload, {{"seq1", "ACGT", ""}});
    const auto frame = lz4f_compress(payload);

    // flags == 0 (no codec) and an unknown reserved bit are both malformed.
    for (uint8_t flags : {0x00, 0x02, 0xFF}) {
        const std::string bin = "bad_flags.bin";
        const auto blob = build_v4_blob(0, payload, frame, std::nullopt, flags);
        write_blob(bin, blob);
        TracEon::SmartStrategy strategy;
        require_throws_with_msg([&] { strategy.loadBinary(bin); },
                                "unsupported codec flags");
        fs::remove(bin);
    }
    // mode byte 2+ is malformed (only GENOME=0 / NGS=1 exist).
    {
        const std::string bin = "bad_mode.bin";
        const auto blob = build_v4_blob(2, payload, frame);
        write_blob(bin, blob);
        TracEon::SmartStrategy strategy;
        require_throws_with_msg([&] { strategy.loadBinary(bin); },
                                "invalid index mode");
        fs::remove(bin);
    }
}


// ── CRC32C unit tests ────────────────────────────────────────────────────────

TEST_CASE("CRC32C known-answer (Castagnoli check value)", "[crc32c][v4]") {
    // Standard CRC-32C check value: crc32c("123456789") == 0xE3069283.
    const char* data = "123456789";
    REQUIRE(TracEon::crc32c(data, 9) == 0xE3069283u);
    REQUIRE(TracEon::crc32c_table_only(data, 9) == 0xE3069283u);
    REQUIRE(TracEon::crc32c_hw_only(data, 9) == 0xE3069283u);

    // Streaming chunking must equal the one-shot result.
    TracEon::Crc32c streamed;
    streamed.update(data, 4);
    streamed.update(data + 4, 5);
    REQUIRE(streamed.finalize() == 0xE3069283u);

    // Empty input: init ^ xorout == 0.
    REQUIRE(TracEon::crc32c(data, 0) == 0u);
}


TEST_CASE("CRC32C fallback vs hardware equivalence over many lengths", "[crc32c][v4]") {
    // Deterministic pseudo-random payload (xorshift64), lengths chosen to
    // exercise every tail branch (0,1,2,3,4,5,7,8,9,15,16,17,31,32,33,
    // 63,64,65, 1000, 1 MiB).
    const std::vector<size_t> lengths = {0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17,
                                         31, 32, 33, 63, 64, 65, 1000, 1u << 20};
    uint64_t state = 0x9E3779B97F4A7C15ull;
    for (size_t len : lengths) {
        std::vector<char> buf(len);
        for (auto& b : buf) {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b = static_cast<char>(state);
        }
        const uint32_t sw = TracEon::crc32c_table_only(buf.data(), len);
        const uint32_t hw = TracEon::crc32c_hw_only(buf.data(), len);
        const uint32_t disp = TracEon::crc32c(buf.data(), len);
        REQUIRE(sw == hw);
        REQUIRE(disp == sw);
    }
}


// ── v5 integrity: mutation / truncation / wrong-field tests ─────────────────

TEST_CASE("v5 single-byte mutation in header → checksum mismatch", "[strategy][binary][v5][crc]") {
    // Flip the index-mode byte (offset 5): it is covered by the CRC, and the
    // checksum is verified BEFORE parsing, so the reader must reject with a
    // checksum error (deterministic).
    const std::string bin = "v4_mut_header.bin";
    TracEon::SmartStrategy strategy;
    strategy.addEntry("a", "ACGT", "");
    strategy.saveBinary(bin);
    auto bytes = read_blob(bin);
    REQUIRE(bytes.size() > 26);
    bytes[5] ^= 0x01; // GENOME -> NGS in the header only
    write_blob(bin, bytes);
    TracEon::SmartStrategy restored;
    require_throws_with_msg([&] { restored.loadBinary(bin); }, "checksum mismatch");
    REQUIRE(restored.getFileCacheSize() == 0);
    fs::remove(bin);
}


TEST_CASE("v5 single-byte mutation of the stored CRC field → checksum mismatch", "[strategy][binary][v5][crc]") {
    const std::string bin = "v4_mut_crc.bin";
    TracEon::SmartStrategy strategy;
    strategy.addEntry("a", "ACGT", "");
    strategy.saveBinary(bin);
    auto bytes = read_blob(bin);
    bytes[23] ^= 0xFF; // corrupt the stored checksum itself
    write_blob(bin, bytes);
    TracEon::SmartStrategy restored;
    require_throws_with_msg([&] { restored.loadBinary(bin); }, "checksum mismatch");
    fs::remove(bin);
}


TEST_CASE("v5 single-byte mutation in the middle of the frame → checksum mismatch", "[strategy][binary][v5][crc]") {
    // LZ4 stores literals raw, so flipping a literal VALUE byte changes one
    // decompressed byte without changing the length: the frame survives and
    // the whole-payload CRC catches the change. Flip candidate bytes inside
    // the frame until one lands on such a literal (a token/length byte flip
    // yields a different rejection, which the scan skips). With unique
    // non-repeating records the frame is almost entirely literals, so at
    // least one literal position is guaranteed to exist.
    const std::string bin = "v4_mut_middle.bin";
    TracEon::SmartStrategy strategy;
    for (int i = 0; i < 40; ++i)
        strategy.addEntry("id" + std::to_string(i), "ACGT" + std::to_string(i) + "TGCAT", "");
    strategy.saveBinary(bin);
    const auto original = read_blob(bin);
    REQUIRE(original.size() > 26);

    bool saw_checksum_error = false;
    for (size_t pos = 27; pos + 4 < original.size() && !saw_checksum_error; ++pos) {
        auto mutated = original;
        mutated[pos] ^= 0x01;
        write_blob(bin, mutated);
        TracEon::SmartStrategy r;
        try {
            r.loadBinary(bin);
        } catch (const std::runtime_error& e) {
            if (std::string(e.what()).find("checksum mismatch") != std::string::npos)
                saw_checksum_error = true;
        }
    }
    REQUIRE(saw_checksum_error);
    fs::remove(bin);
}


TEST_CASE("v5 single-byte mutation at the tail (end of frame) → rejected", "[strategy][binary][v5][crc]") {
    // The final byte of a v5 file is the tail of the LZ4 Frame (the frame's
    // own content checksum), so a tail mutation is caught by the frame
    // validator (or, if the frame survives, by the payload CRC). Either way
    // the file must be rejected — never silently accepted.
    const std::string bin = "v4_mut_tail.bin";
    TracEon::SmartStrategy strategy;
    strategy.addEntry("a", "ACGT", "");
    strategy.saveBinary(bin);
    auto bytes = read_blob(bin);
    REQUIRE(bytes.size() > 26);
    bytes.back() ^= 0x01;
    write_blob(bin, bytes);
    TracEon::SmartStrategy restored;
    REQUIRE_THROWS_AS(restored.loadBinary(bin), std::runtime_error);
    REQUIRE(restored.getFileCacheSize() == 0);
    fs::remove(bin);
}


TEST_CASE("v5 truncation at frame boundaries → rejected", "[strategy][binary][v5][truncate]") {
    // Cut a valid file at several points: header end (26), a few bytes into
    // the frame, mid-frame, and one byte short of the end. Every truncated
    // variant must be rejected — truncation is detected by the exact-length /
    // exact-frame-termination requirements, never silently accepted.
    const std::string bin = "v4_trunc.bin";
    TracEon::SmartStrategy strategy;
    strategy.addEntry("a", "ACGT", "");
    strategy.addEntry("b", "TTTT", "");
    strategy.saveBinary(bin);
    const auto original = read_blob(bin);
    REQUIRE(original.size() > 40);

    const size_t cuts[] = {27, 31, original.size() / 2, original.size() - 1};
    for (size_t cut : cuts) {
        write_blob(bin, std::vector<char>(original.begin(), original.begin() + static_cast<long>(cut)));
        TracEon::SmartStrategy r;
        REQUIRE_THROWS_AS(r.loadBinary(bin), std::runtime_error);
        REQUIRE(r.getFileCacheSize() == 0);
    }
    fs::remove(bin);
}


TEST_CASE("v5 wrong logical length → exact-length requirement rejects", "[strategy][binary][v5][length]") {
    const std::string bin = "v4_wrong_len.bin";
    TracEon::SmartStrategy strategy;
    strategy.addEntry("a", "ACGT", "");
    strategy.saveBinary(bin);
    auto bytes = read_blob(bin);

    // Patch logical_len (offset 7..14) to +1: the frame decompresses to the
    // real size, so dst_remaining ends at 1 → size mismatch. (Checked before
    // the CRC, so the rejection is deterministic.)
    {
        auto bigger = bytes;
        uint64_t len = 0;
        for (int i = 0; i < 8; ++i) len |= static_cast<uint64_t>(static_cast<uint8_t>(bigger[7 + i])) << (8 * i);
        for (int i = 0; i < 8; ++i) bigger[7 + i] = static_cast<char>(((len + 1) >> (8 * i)) & 0xFFu);
        write_blob(bin, bigger);
        TracEon::SmartStrategy r;
        require_throws_with_msg([&] { r.loadBinary(bin); }, "decompressed size mismatch");
    }
    // Patch logical_len to a smaller value: the frame cannot fully decompress
    // into the undersized arena → frame not terminated at declared boundary.
    {
        auto smaller = bytes;
        uint64_t len = 0;
        for (int i = 0; i < 8; ++i) len |= static_cast<uint64_t>(static_cast<uint8_t>(smaller[7 + i])) << (8 * i);
        for (int i = 0; i < 8; ++i) smaller[7 + i] = static_cast<char>(((len - 1) >> (8 * i)) & 0xFFu);
        write_blob(bin, smaller);
        TracEon::SmartStrategy r;
        REQUIRE_THROWS_AS(r.loadBinary(bin), std::runtime_error);
    }
    fs::remove(bin);
}


TEST_CASE("v4 wrong record count → rejected (implausible / corrupt)", "[strategy][binary][v4][count]") {
    // (a) count = 3 with one record serialized: remaining payload bytes
    // cannot cover 3 records → count-bounded reserve rejects up front.
    {
        std::vector<char> payload;
        write_genome_payload(payload, {{"seq1", "ACGT", ""}}); // count == 1
        // Patch the count field (first 8 payload bytes) to 3.
        uint64_t count = 3;
        for (int i = 0; i < 8; ++i) payload[static_cast<size_t>(i)] = static_cast<char>((count >> (8 * i)) & 0xFFu);
        const auto frame = lz4f_compress(payload);
        const std::string bin = "v4_bad_count1.bin";
        write_blob(bin, build_v4_blob(0, payload, frame));
        TracEon::SmartStrategy r;
        require_throws_with_msg([&] { r.loadBinary(bin); }, "record count implausible");
        fs::remove(bin);
    }
    // (b) count = 2 with one record + 8 junk trailing bytes: the byte-bound
    // reserve passes (24 remaining ≥ 2×12), so the second record's length
    // field is read from junk → payload-corruption rejection mid-parse.
    {
        std::vector<char> payload;
        write_genome_payload(payload, {{"seq1", "ACGT", ""}});
        uint64_t count = 2;
        for (int i = 0; i < 8; ++i) payload[static_cast<size_t>(i)] = static_cast<char>((count >> (8 * i)) & 0xFFu);
        payload.insert(payload.end(), 8, static_cast<char>(0xFF)); // junk: len = 0xFFFFFFFF
        const auto frame = lz4f_compress(payload);
        const std::string bin = "v4_bad_count2.bin";
        write_blob(bin, build_v4_blob(0, payload, frame));
        TracEon::SmartStrategy r;
        require_throws_with_msg([&] { r.loadBinary(bin); }, "payload is corrupt");
        fs::remove(bin);
    }
}


// ── Hardening regression tests, ported from v1/v2/v3 to v4 semantics ────────
// The v2 class of INT_MAX int-truncation bugs cannot occur in v4 (the reader
// is 64-bit-clean LZ4F streaming), so the ported probes target the equivalent
// guards: implausible logical length rejected before allocation, oversized
// frame length rejected by the subtraction-form bounds check, and the
// count-bomb rejected by count-bounded reserve.

TEST_CASE("loadBinary v4 rejects 2^64-1 logical length before allocating (truncation-bypass port)",
          "[strategy][binary][vuln][v4]") {
    // vuln-0007 analogue: the old v2 reader could be tricked by original_size
    // = 2^32+N truncating in the int LZ4 API. v4 declares logical_len as u64;
    // an absurd value (2^64-1) is rejected by the OOM guard (it always
    // exceeds half of available memory on any real host) BEFORE any
    // allocation attempt.
    std::vector<char> payload;
    write_genome_payload(payload, {{"seq1", "ACGT", ""}});
    const auto frame = lz4f_compress(payload);
    const std::string bin = "vuln0007_v4_logical.bin";
    write_blob(bin, build_v4_blob(0, payload, frame, std::nullopt, std::nullopt,
                                  /*logical_len=*/UINT64_MAX));
    TracEon::SmartStrategy strategy;
    require_throws_with_msg([&] { strategy.loadBinary(bin); },
                            "OOM guard: refusing to decompress");
    REQUIRE(strategy.getFileCacheSize() == 0);
    fs::remove(bin);
}


TEST_CASE("loadBinary v4 rejects undersized logical length (< 8: no record-count field)",
          "[strategy][binary][vuln][v4]") {
    std::vector<char> payload;
    write_genome_payload(payload, {{"seq1", "ACGT", ""}});
    const auto frame = lz4f_compress(payload);
    const std::string bin = "vuln_v4_small_logical.bin";
    // logical_len = 4: every payload begins with the 8-byte count, so this is
    // always malformed.
    write_blob(bin, build_v4_blob(0, payload, frame, std::nullopt, std::nullopt,
                                  /*logical_len=*/4));
    TracEon::SmartStrategy strategy;
    require_throws_with_msg([&] { strategy.loadBinary(bin); },
                            "implausible decompressed size");
    fs::remove(bin);
}


TEST_CASE("loadBinary v4 rejects oversized frame length (2^64-1) via subtraction-form bounds check",
          "[strategy][binary][vuln][v4]") {
    // Bug C analogue: frame_len = 2^64-1 used to pass `ptr + n > end` (the
    // addition wraps) in the old readers. The subtraction-form check
    // `frame_len > end - ptr` rejects it with no pointer arithmetic.
    std::vector<char> payload;
    write_genome_payload(payload, {{"seq1", "ACGT", ""}});
    const auto frame = lz4f_compress(payload);
    const std::string bin = "bugc_v4_frame.bin";
    write_blob(bin, build_v4_blob(0, payload, frame, std::nullopt, std::nullopt,
                                  std::nullopt, /*frame_len=*/UINT64_MAX));
    TracEon::SmartStrategy strategy;
    require_throws_with_msg([&] { strategy.loadBinary(bin); },
                            "Binary cache is truncated");
    REQUIRE(strategy.getFileCacheSize() == 0);
    fs::remove(bin);
}


TEST_CASE("loadBinary v4 bounds map.reserve() by remaining payload bytes (count=5e8)",
          "[strategy][binary][vuln][v4]") {
    // vuln-0009 port: count = 500M with a tiny payload must throw 'record
    // count implausible' instead of attempting a ~40 GB map.reserve().
    std::vector<char> payload;
    uint64_t count = 500'000'000ULL;
    for (int i = 0; i < 8; ++i) payload.push_back(static_cast<char>((count >> (8 * i)) & 0xFFu));
    uint32_t len = 1;
    payload.insert(payload.end(), reinterpret_cast<const char*>(&len), reinterpret_cast<const char*>(&len) + 4);
    payload.push_back('a');
    len = 4;
    payload.insert(payload.end(), reinterpret_cast<const char*>(&len), reinterpret_cast<const char*>(&len) + 4);
    const char seq[] = "ACGT";
    payload.insert(payload.end(), seq, seq + 4);
    len = 0;
    payload.insert(payload.end(), reinterpret_cast<const char*>(&len), reinterpret_cast<const char*>(&len) + 4);

    const auto frame = lz4f_compress(payload);
    const std::string bin = "vuln0009_v4_count.bin";
    write_blob(bin, build_v4_blob(0, payload, frame));
    TracEon::SmartStrategy strategy;
    require_throws_with_msg([&] { strategy.loadBinary(bin); },
                            "Binary cache record count implausible");
    REQUIRE(strategy.getFileCacheSize() == 0);
    fs::remove(bin);
}


TEST_CASE("loadBinary rejects old MMAP format", "[strategy][errors]") {
    std::string old_file = "old_format.bin";
    {
        std::ofstream out(old_file, std::ios::binary);
        out.write("MMAP", 4); // old v1.0 magic
        uint8_t mode = 0; out.write(reinterpret_cast<const char*>(&mode), 1);
        uint64_t count = 0; out.write(reinterpret_cast<const char*>(&count), 8);
    }
    TracEon::SmartStrategy strategy;
    REQUIRE_THROWS_AS(strategy.loadBinary(old_file), std::runtime_error);
    fs::remove(old_file);
}


// ─── v5 Binary Format — save/restore round-trip + lifecycle ─────────────────

TEST_CASE("v5 binary format save → restore round-trip via public API", "[strategy][binary][v5]") {
    const std::string src = "tmp_v5_api_src.fasta";
    const std::string bin = "tmp_v5_api.bin";
    {
        std::ofstream out(src);
        out << ">chr1\nACGTACGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAA\n";
    }
    {
        TracEon::Cache cache;
        cache.loadFile(src);
        REQUIRE(cache.size() == 2);
        cache.save(bin);
    }
    {
        TracEon::Cache cache;
        cache.restore(bin);
        REQUIRE(cache.size() == 2);
        REQUIRE(cache.get("chr1") == "ACGTACGTACGTACGT");
        REQUIRE(cache.get("chr2") == "TTTTGGGGCCCCAAAA");
        REQUIRE(cache.getView("chr1") == "ACGTACGTACGTACGT");
    }
    fs::remove(src);
    fs::remove(bin);
}

TEST_CASE("binary restore rejects bytes after the declared frame", "[strategy][binary][eof]") {
    const std::string bin = "v5_trailing_bytes.bin";
    TracEon::SmartStrategy writer;
    writer.addEntry("x", "ACGT", "");
    writer.saveBinary(bin);
    {
        std::ofstream out(bin, std::ios::binary | std::ios::app);
        out << "TRAIL";
        REQUIRE(out.good());
    }
    TracEon::SmartStrategy restored;
    require_throws_with_msg([&] { restored.loadBinary(bin); }, "trailing bytes");
    REQUIRE(restored.getFileCacheSize() == 0);
    REQUIRE(restored.getDetectedFormat() == TracEon::FileFormat::UNKNOWN);
    fs::remove(bin);
}

TEST_CASE("legacy v4 restore rejects verified appended bytes", "[strategy][binary][v4][eof]") {
    std::vector<char> payload;
    write_genome_payload(payload, {{"x", "ACGT", ""}});
    const auto frame = lz4f_compress(payload);
    auto blob = build_v4_blob(0, payload, frame);
    blob.insert(blob.end(), {'T','R','A','I','L'});
    const std::string bin = "v4_trailing_bytes.bin";
    write_blob(bin, blob);
    REQUIRE(fs::file_size(bin) == blob.size());
    TracEon::SmartStrategy restored;
    require_throws_with_msg([&] { restored.loadBinary(bin); }, "trailing bytes");
    REQUIRE(restored.getFileCacheSize() == 0);
    REQUIRE(restored.getArenaBytes() == 0);
    fs::remove(bin);
}

TEST_CASE("failed restore does not change later text-load index mode", "[strategy][binary][mode]") {
    const std::string bad = "bad_ngs_mode_cache.bin";
    std::vector<char> payload;
    write_genome_payload(payload, {{"seq1", "ACGT", ""}});
    const auto frame = lz4f_compress(payload);
    auto blob = build_v4_blob(1, payload, frame);
    blob.back() ^= 0x5a;
    write_blob(bad, blob);

    TracEon::SmartStrategy strategy(TracEon::IndexMode::GENOME);
    REQUIRE_THROWS_AS(strategy.loadBinary(bad), std::runtime_error);
    REQUIRE(strategy.getIndexMode() == TracEon::IndexMode::GENOME);

    const std::string fasta = "mode_after_failed_restore.fasta";
    { std::ofstream out(fasta); out << ">id\nACGT\n"; }
    strategy.loadFile(fasta);
    REQUIRE(strategy.getIndexMode() == TracEon::IndexMode::GENOME);
    REQUIRE(strategy.getSequence("id") == "ACGT");
    fs::remove(bad);
    fs::remove(fasta);
}

TEST_CASE("legacy v4 NGS cache is read by identifier rather than persisted hash",
          "[strategy][binary][v4][ngs][compat]") {
    std::vector<char> payload;
    auto append_native = [&](const auto& value) {
        const char* p = reinterpret_cast<const char*>(&value);
        payload.insert(payload.end(), p, p + sizeof(value));
    };
    const uint64_t count = 1;
    const uint64_t deliberately_wrong_hash = 0;
    const uint32_t id_len = 4, seq_len = 4, qual_len = 4;
    append_native(count);
    append_native(deliberately_wrong_hash);
    append_native(id_len); payload.insert(payload.end(), {'r','e','a','d'});
    append_native(seq_len); payload.insert(payload.end(), {'A','C','G','T'});
    append_native(qual_len); payload.insert(payload.end(), {'I','I','I','I'});
    const auto frame = lz4f_compress(payload);
    const std::string bin = "legacy_v4_ngs.traceon";
    write_blob(bin, build_v4_blob(1, payload, frame));

    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getIndexMode() == TracEon::IndexMode::NGS);
    REQUIRE(restored.getSequence("read") == "ACGT");
    REQUIRE(restored.getQuality("read") == "IIII");
    fs::remove(bin);
}

TEST_CASE("v5 persists exact detected FASTQ format", "[strategy][binary][v5][format]") {
    const std::string src = "protein_format.fastq";
    const std::string bin = "protein_format.traceon";
    { std::ofstream out(src); out << "@p\nEFILP\n+\nIIIII\n"; }
    TracEon::SmartStrategy writer(TracEon::IndexMode::NGS);
    writer.loadFile(src);
    REQUIRE(writer.getDetectedFormat() == TracEon::FileFormat::PROTEIN_FASTQ);
    writer.saveBinary(bin);
    TracEon::SmartStrategy restored;
    restored.loadBinary(bin);
    REQUIRE(restored.getDetectedFormat() == TracEon::FileFormat::PROTEIN_FASTQ);
    REQUIRE(restored.getIndexMode() == TracEon::IndexMode::NGS);
    REQUIRE(restored.getSequence("p") == "EFILP");
    fs::remove(src);
    fs::remove(bin);
}

TEST_CASE("CRC-valid duplicate identifiers are rejected in v4 and v5",
          "[strategy][binary][duplicates]") {
    const std::string root = std::string(TRACEON_SOURCE_DIR) +
        "/fuzz/corpus/v4_loader/";
    for (const std::string& name : {
             "v4_duplicate_ids.traceon", "v4_duplicate_ids_ngs.traceon",
             "v5_duplicate_ids.traceon", "v5_duplicate_ids_ngs.traceon"}) {
        CAPTURE(name);
        TracEon::SmartStrategy strategy;
        strategy.addEntry("old", "ACGT", "");
        REQUIRE_THROWS_AS(strategy.loadBinary(root + name), std::runtime_error);
        REQUIRE(strategy.getFileCacheSize() == 0);
        REQUIRE(strategy.getDetectedFormat() == TracEon::FileFormat::UNKNOWN);
    }
}

TEST_CASE("mapped v5 preserves records, lifecycle and conversions", "[binary][mapped]") {
    for (auto mode : {TracEon::IndexMode::GENOME, TracEon::IndexMode::NGS}) {
        const std::string bin = "mapped_roundtrip.traceon", compressed = "mapped_to_lz4.traceon";
        TracEon::SmartStrategy writer(mode);
        writer.addEntry("read", "ACGTN", "!@III");
        writer.addEntry("empty", "", "");
        writer.saveBinary(bin, TracEon::CacheStorage::Mapped);
        REQUIRE(read_blob(bin)[4] == 0);
        TracEon::SmartStrategy mapped;
        mapped.loadBinary(bin);
        REQUIRE(mapped.getIndexMode() == mode);
        REQUIRE(mapped.getSequence("read") == "ACGTN");
        REQUIRE(mapped.getQuality("read") == "!@III");
        REQUIRE(mapped.hasSequence("empty"));
        REQUIRE_FALSE(mapped.hasSequence("missing"));
        REQUIRE(mapped.getArenaBytes() == fs::file_size(bin));
        REQUIRE_THROWS(mapped.addEntry("late", "A", ""));
        mapped.saveBinary(compressed);
        mapped.loadBinary(compressed);
        REQUIRE(mapped.getQuality("read") == "!@III");
        mapped.saveBinary(bin, TracEon::CacheStorage::Mapped);
        mapped.loadBinary(bin);
        mapped.clearCache();
        REQUIRE(mapped.getArenaBytes() == 0);
        REQUIRE(mapped.getFileCacheSize() == 0);
        fs::remove(bin); fs::remove(compressed);
    }
    TracEon::Cache empty, restored;
    empty.save("mapped_empty.traceon", TracEon::CacheStorage::Mapped);
    restored.restore("mapped_empty.traceon");
    REQUIRE(restored.size() == 0);
    restored.clearCache(); fs::remove("mapped_empty.traceon");
}

TEST_CASE("mapped v5 rejects corruption before publication", "[binary][mapped]") {
    const std::string bin = "mapped_bad.traceon";
    TracEon::Cache writer;
    writer.set("id", "ACGT"); writer.save(bin, TracEon::CacheStorage::Mapped);
    const auto good = read_blob(bin);
    auto reject = [&](const std::vector<char>& bytes) {
        write_blob(bin, bytes);
        TracEon::Cache reader; reader.set("old", "TT");
        REQUIRE_THROWS(reader.restore(bin));
        REQUIRE(reader.size() == 0);
        REQUIRE(reader.getArenaBytes() == 0);
        REQUIRE_FALSE(reader.hasSequence("old"));
    };
    // Every truncated prefix, every changed byte, and a physical suffix.
    for (size_t i = 0; i < good.size(); ++i) {
        CAPTURE(i);
        reject(std::vector<char>(good.begin(), good.begin() + i));
        auto changed = good; changed[i] ^= 0x40; reject(changed);
    }
    auto suffix = good; suffix.push_back('x'); reject(suffix);
    // Independently assemble CRC-valid payloads to exercise structural validation.
    auto blob = [](const std::vector<char>& payload) {
        std::vector<char> b={'T','R','O',5,0,0,char(0xff)};
        append_le64(b,payload.size()); append_le64(b,payload.size());
        TracEon::Crc32c crc; crc.update(payload.data(),payload.size()); crc.update(b.data(),b.size());
        append_le32(b,crc.finalize()); b.insert(b.end(),payload.begin(),payload.end()); return b;
    };
    std::vector<char> payload; append_le64(payload,2);
    for (int i=0;i<2;++i) {
        append_le32(payload,1); payload.push_back('x');
        append_le32(payload,1); payload.push_back('A'); append_le32(payload,0);
    }
    reject(blob(payload)); // duplicate IDs with a valid checksum
    payload.clear(); append_le64(payload,1); append_le32(payload,0xffffffffu);
    append_le32(payload,0); append_le32(payload,0); reject(blob(payload));
    payload.clear(); append_le64(payload,0); payload.push_back('x'); reject(blob(payload));
    payload.clear(); append_le64(payload,0xffffffffffffffffULL); reject(blob(payload));
    fs::remove(bin);
}

#ifndef _WIN32
TEST_CASE("mapped views survive atomic cache replacement", "[binary][mapped][atomic]") {
    const std::string bin="mapped_replace.traceon";
    TracEon::Cache first; first.set("id","OLD");
    first.save(bin,TracEon::CacheStorage::Mapped);
    TracEon::Cache reader; reader.restore(bin); auto view=reader.getView("id");
    // Re-saving onto the same pathname must not truncate its mapped inode.
    reader.save(bin,TracEon::CacheStorage::Mapped);
    REQUIRE(view=="OLD");
    TracEon::Cache second; second.set("id","NEW");
    second.save(bin,TracEon::CacheStorage::Mapped);
    TracEon::Cache latest; latest.restore(bin);
    REQUIRE(latest.getView("id")=="NEW");
    REQUIRE(view=="OLD"); REQUIRE(reader.getView("id")=="OLD");
    reader.clearCache(); latest.clearCache(); fs::remove(bin);
}
#endif

TEST_CASE("mapped FASTQ preserves format, quality and empty reads", "[binary][mapped][fastq]") {
    const std::string input="mapped.fastq", bin="mapped_fastq.traceon";
    { std::ofstream out(input); out << "@p\nEFILP\n+\n!@+IJ\n@empty\n\n+\n\n"; }
    TracEon::Cache cache(TracEon::IndexMode::NGS);
    cache.loadFile(input); cache.save(bin,TracEon::CacheStorage::Mapped);
    cache.restore(bin);
    REQUIRE(cache.getDetectedFormat()==TracEon::FileFormat::PROTEIN_FASTQ);
    REQUIRE(cache.getIndexMode()==TracEon::IndexMode::NGS);
    REQUIRE(cache.getFastqRecord("p")->sequence=="EFILP");
    REQUIRE(cache.getFastqRecord("p")->quality=="!@+IJ");
    REQUIRE(cache.getFastqRecord("empty").has_value());
    REQUIRE(cache.getFastqRecord("empty")->sequence.empty());
    REQUIRE_FALSE(cache.getFastqRecord("missing").has_value());
    cache.clearCache(); fs::remove(input);fs::remove(bin);
}
