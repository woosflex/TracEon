// TestHelpers.h — shared includes + fixtures/helpers for the
// per-domain TracEon unit tests (split from SmartStrategyTests.cpp).
// All helpers are `inline` so every test TU gets its own copy and
// no separate .cpp is required.
#pragma once

#include <catch2/catch_test_macros.hpp>
#include "SmartStrategy.h"
#include "Cache.h"
#include "SimdUtils.h"
#include "Crc32c.h"
#include <fstream>
#include <filesystem>
#include <thread>
#include <atomic>
#include <zlib.h>
#include <lz4.h>
#include <lz4frame.h>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <functional>
#include <string>
#include <optional>
#include <tuple>
#include <vector>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace fs = std::filesystem;

// ── Shared v4 binary-cache crafting helpers ──────────────────────
// ── `.traceon` v4 binary format helpers ─────────────────────────────────────
// v2.0.0 ships v4 ("TRO\x04") as the ONLY readable binary format: the
// v1/v2/v3 readers were removed (design review Q1/Q5 — no legacy readers).
// These helpers hand-build v4 blobs for the corruption / truncation /
// hardening tests below, mirroring the header layout in src/SmartStrategy.cpp:
//
//   magic "TRO\x04" | codec flags(1, 0x01=LZ4 Frame) | index mode(1) |
//   logical length(u64 LE) | frame length(u64 LE) | CRC32C(u32 LE) |
//   LZ4 Frame (frame length bytes)
//
// The payload layout is unchanged from v1–v3 (count u64 + records), so
// write_genome_payload() below is shared with the v4 crafting helpers.

// Helper: assert `fn` throws std::runtime_error whose message contains `substr`.
static void require_throws_with_msg(const std::function<void()>& fn,
                                    const std::string& substr) {
    bool threw = false;
    try {
        fn();
    } catch (const std::runtime_error& e) {
        threw = true;
        REQUIRE(std::string(e.what()).find(substr) != std::string::npos);
    }
    REQUIRE(threw);
}

static void write_genome_payload(std::vector<char>& buf, const std::vector<std::tuple<std::string, std::string, std::string>>& records) {
    uint64_t count = records.size();
    buf.insert(buf.end(), reinterpret_cast<const char*>(&count), reinterpret_cast<const char*>(&count) + sizeof(count));
    for (const auto& [id, seq, qual] : records) {
        auto append_field = [&](const std::string& field) {
            uint32_t len = static_cast<uint32_t>(field.size());
            buf.insert(buf.end(), reinterpret_cast<const char*>(&len), reinterpret_cast<const char*>(&len) + sizeof(len));
            buf.insert(buf.end(), field.data(), field.data() + len);
        };
        append_field(id);
        append_field(seq);
        append_field(qual);
    }
}

static void append_le64(std::vector<char>& buf, uint64_t v) {
    for (int i = 0; i < 8; ++i) buf.push_back(static_cast<char>((v >> (8 * i)) & 0xFFu));
}
static void append_le32(std::vector<char>& buf, uint32_t v) {
    for (int i = 0; i < 4; ++i) buf.push_back(static_cast<char>((v >> (8 * i)) & 0xFFu));
}

// Compress a payload into an LZ4 Frame (same defaults as saveBinary()).
static std::vector<char> lz4f_compress(const std::vector<char>& payload) {
    LZ4F_preferences_t prefs;
    std::memset(&prefs, 0, sizeof(prefs));
    const size_t cap = LZ4F_compressFrameBound(payload.size(), &prefs);
    std::vector<char> frame(cap);
    const size_t fsz = LZ4F_compressFrame(frame.data(), cap, payload.data(), payload.size(), &prefs);
    REQUIRE(!LZ4F_isError(fsz));
    frame.resize(fsz);
    return frame;
}

// Build a v4 blob. Overrides let tests corrupt individual header fields; the
// CRC is (re)computed over the final header + payload unless explicitly given,
// so `crc` overrides are how tests corrupt the checksum field itself.
static std::vector<char> build_v4_blob(uint8_t mode, const std::vector<char>& payload,
                                       const std::vector<char>& frame,
                                       std::optional<uint8_t> magic_ver = std::nullopt,
                                       std::optional<uint8_t> flags = std::nullopt,
                                       std::optional<uint64_t> logical_len = std::nullopt,
                                       std::optional<uint64_t> frame_len = std::nullopt,
                                       std::optional<uint32_t> crc = std::nullopt) {
    std::vector<char> blob;
    const uint8_t ver = magic_ver.value_or(0x04);
    blob.insert(blob.end(), {'T', 'R', 'O', static_cast<char>(ver)});
    blob.push_back(static_cast<char>(flags.value_or(0x01)));
    blob.push_back(static_cast<char>(mode));
    append_le64(blob, logical_len.value_or(payload.size()));
    append_le64(blob, frame_len.value_or(frame.size()));
    // CRC32C over the ENTIRE uncompressed payload followed by the canonical
    // header fields [0..22) (the checksum field itself excluded) — matching
    // saveBinary()/loadBinary()'s streaming order (payload first, header
    // last, since frame_len is only known after the frame is complete).
    const uint32_t c = crc.value_or([&]() {
        TracEon::Crc32c c32;
        c32.update(payload.data(), payload.size());
        c32.update(blob.data(), 22);
        return c32.finalize();
    }());
    append_le32(blob, c);
    blob.insert(blob.end(), frame.begin(), frame.end());
    return blob;
}

static void write_blob(const std::string& path, const std::vector<char>& blob) {
    std::ofstream out(path, std::ios::binary);
    out.write(blob.data(), static_cast<std::streamsize>(blob.size()));
}

static std::vector<char> read_blob(const std::string& path) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    const size_t sz = static_cast<size_t>(in.tellg());
    in.seekg(0);
    std::vector<char> bytes(sz);
    in.read(bytes.data(), static_cast<std::streamsize>(sz));
    return bytes;
}

// ── Multithreaded parser test helpers (>10 MB fixtures) ──────────
// ═══════════════════════════════════════════════════════════════════════════
// Multithreaded SIMD Parser Tests
// ═══════════════════════════════════════════════════════════════════════════
//
// These tests exercise parseFastaMultithreadedTemplate() and
// parseFastqMultithreadedTemplate() via the public loadFile() API by creating
// content larger than MULTITHREAD_THRESHOLD (10 MB).
//
// The threshold check lives in parseArena():
//   if (content.size() > 10*1024*1024) {
//       if (isFastq) parseFastqMultithreadedTemplate(...);
//       else         parseFastaMultithreadedTemplate(...);
//   }
//
// ── Helpers ───────────────────────────────────────────────────────────────────

namespace {
    // Build a FASTA string larger than 10 MB with record_count records.
    // Each record: ">record_N\n{pattern}\n"
    // pattern is repeated_chunk repeated to fill approx seq_len bytes per record.
    std::string make_large_fasta(size_t target_bytes, size_t seq_len = 4000) {
        const size_t REPEATS = seq_len / 4;
        const std::string chunk = "ACGT";
        std::string seq;
        seq.reserve(seq_len);
        for (size_t i = 0; i < REPEATS; ++i) seq += chunk;

        std::string result;
        result.reserve(target_bytes + 4096);
        int n = 1;
        while (result.size() < target_bytes) {
            result += ">record_" + std::to_string(n) + "\n";
            result += seq;
            result += '\n';
            ++n;
        }
        return result;
    }

    // Build a FASTQ string larger than 10 MB with record_count records.
    std::string make_large_fastq(size_t target_bytes, size_t seq_len = 4000) {
        const size_t REPEATS = seq_len / 4;
        const std::string seq = std::string(REPEATS, 'A') + std::string(REPEATS, 'C') +
                                std::string(REPEATS, 'G') + std::string(REPEATS, 'T');
        const std::string qual(seq_len, 'I');

        std::string result;
        result.reserve(target_bytes + 4096);
        int n = 1;
        while (result.size() < target_bytes) {
            result += "@record_" + std::to_string(n) + "\n";
            result += seq;
            result += '\n';
            result += "+\n";
            result += qual;
            result += '\n';
            ++n;
        }
        return result;
    }
} // anonymous namespace

// ── GZIP byte helpers ────────────────────────────────────────────
// Helper: compress bytes with zlib-ng gzwrite into a temp file, read back raw bytes
static std::vector<char> compress_to_gzip_bytes(const std::string& data, const std::string& tmp_path) {
    gzFile f = gzopen(tmp_path.c_str(), "wb");
    gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
    gzclose(f);

    std::ifstream in(tmp_path, std::ios::binary | std::ios::ate);
    size_t sz = in.tellg();
    in.seekg(0);
    std::vector<char> bytes(sz);
    in.read(bytes.data(), sz);
    return bytes;
}
