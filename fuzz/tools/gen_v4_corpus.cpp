// gen_v4_corpus.cpp — regenerates the `.traceon` v4-migration/v5 seed corpus.
// Real v5 files go through SmartStrategy::saveBinary(); crafted v4 GENOME and
// NGS files retain migration coverage alongside malformed variants.
//
// The v5 files go through the real SmartStrategy::saveBinary() path; the
// crafted variants are built with the same helpers the tests use (LZ4 Frame
// via LZ4F_compressFrame + whole-payload CRC32C via the project's Crc32c.h).
//
// Build (from the repo root; links the project's lz4 static lib):
//   g++ -std=c++20 -O2 -I include -I third_party/lz4/lib \
//       fuzz/tools/gen_v4_corpus.cpp \
//       build/third_party/lz4/build/cmake/liblz4.a -o /tmp/gen_v4_corpus
//   /tmp/gen_v4_corpus fuzz/corpus/v4_loader
//
// Run from the repo root so the include paths resolve.

#include "Crc32c.h"
#include "SmartStrategy.h"

#include <lz4frame.h>

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace {

void write_blob(const std::string& path, const std::vector<char>& blob) {
    std::ofstream out(path, std::ios::binary);
    out.write(blob.data(), static_cast<std::streamsize>(blob.size()));
}

void append_le64(std::vector<char>& buf, uint64_t v) {
    for (int i = 0; i < 8; ++i) buf.push_back(static_cast<char>((v >> (8 * i)) & 0xFFu));
}
void append_le32(std::vector<char>& buf, uint32_t v) {
    for (int i = 0; i < 4; ++i) buf.push_back(static_cast<char>((v >> (8 * i)) & 0xFFu));
}

void write_genome_payload(std::vector<char>& buf,
                          const std::vector<std::tuple<std::string, std::string, std::string>>& records) {
    uint64_t count = records.size();
    buf.insert(buf.end(), reinterpret_cast<const char*>(&count),
               reinterpret_cast<const char*>(&count) + sizeof(count));
    for (const auto& [id, seq, qual] : records) {
        auto append_field = [&](const std::string& field) {
            uint32_t len = static_cast<uint32_t>(field.size());
            buf.insert(buf.end(), reinterpret_cast<const char*>(&len),
                       reinterpret_cast<const char*>(&len) + sizeof(len));
            buf.insert(buf.end(), field.data(), field.data() + len);
        };
        append_field(id);
        append_field(seq);
        append_field(qual);
    }
}

std::vector<char> lz4f_compress(const std::vector<char>& payload) {
    LZ4F_preferences_t prefs;
    std::memset(&prefs, 0, sizeof(prefs));
    const size_t cap = LZ4F_compressFrameBound(payload.size(), &prefs);
    std::vector<char> frame(cap);
    const size_t fsz = LZ4F_compressFrame(frame.data(), cap, payload.data(),
                                          payload.size(), &prefs);
    if (LZ4F_isError(fsz)) { std::abort(); }
    frame.resize(fsz);
    return frame;
}

// Mirror of build_v4_blob in tests/TestHelpers.h.
std::vector<char> build_v4_blob(uint8_t mode, const std::vector<char>& payload,
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

void write_v5_payload(std::vector<char>& buf,
                      const std::vector<std::tuple<std::string, std::string, std::string>>& records) {
    append_le64(buf, records.size());
    for (const auto& [id, seq, qual] : records) {
        for (const auto* field : {&id, &seq, &qual}) {
            append_le32(buf, static_cast<uint32_t>(field->size()));
            buf.insert(buf.end(), field->begin(), field->end());
        }
    }
}

std::vector<char> build_v5_blob(uint8_t mode, uint8_t format,
                                const std::vector<char>& payload, bool mapped = false) {
    const auto frame = mapped ? payload : lz4f_compress(payload);
    std::vector<char> blob{'T','R','O',static_cast<char>(0x05),
                           static_cast<char>(mapped ? 0x00 : 0x01), static_cast<char>(mode),
                           static_cast<char>(format)};
    append_le64(blob, payload.size());
    append_le64(blob, frame.size());
    TracEon::Crc32c crc;
    crc.update(payload.data(), payload.size());
    crc.update(blob.data(), 23);
    append_le32(blob, crc.finalize());
    blob.insert(blob.end(), frame.begin(), frame.end());
    return blob;
}

} // namespace

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : "fuzz/corpus/v4_loader";

    // ── Valid v5 round-trip caches via the real saveBinary() path ──────────
    {
        TracEon::SmartStrategy s;
        s.addEntry("chr1", "ACGTACGTACGTACGT", "");
        s.addEntry("chr2", "TTTTGGGGCCCCAAAA", "");
        s.saveBinary(dir + "/v5_valid_genome.traceon");
        s.saveBinary(dir + "/v5_mapped_genome.traceon", TracEon::CacheStorage::Mapped);
    }
    {
        TracEon::SmartStrategy s(TracEon::IndexMode::NGS);
        s.addEntry("read1", "ACGT", "IIII");
        s.addEntry("read2", "TGCA", "####");
        s.saveBinary(dir + "/v5_valid_ngs.traceon");
        s.saveBinary(dir + "/v5_mapped_ngs.traceon", TracEon::CacheStorage::Mapped);
    }

    // Directed v5 seeds with valid CRCs reach payload validation.
    std::vector<char> v5_one;
    write_v5_payload(v5_one, {{"seq1", "ACGT", ""}});
    write_blob(dir + "/v5_invalid_format.traceon", build_v5_blob(0, 0x7f, v5_one));
    std::vector<char> v5_dup;
    write_v5_payload(v5_dup, {{"dup", "ACGT", ""}, {"dup", "TGCA", ""}});
    write_blob(dir + "/v5_duplicate_ids.traceon", build_v5_blob(0, 0x00, v5_dup));
    write_blob(dir + "/v5_mapped_duplicates.traceon", build_v5_blob(0, 0x00, v5_dup, true));
    write_blob(dir + "/v5_duplicate_ids_ngs.traceon", build_v5_blob(1, 0x03, v5_dup));
    auto trailing_payload = v5_one;
    trailing_payload.push_back('X');
    write_blob(dir + "/v5_payload_trailing.traceon", build_v5_blob(0, 0x00, trailing_payload));
    std::vector<char> count_boundary(8, static_cast<char>(0xff));
    write_blob(dir + "/v5_count_boundary.traceon", build_v5_blob(0, 0x00, count_boundary));
    std::vector<char> length_boundary;
    append_le64(length_boundary, 1);
    append_le32(length_boundary, 0xffffffffu);
    write_blob(dir + "/v5_length_boundary.traceon", build_v5_blob(0, 0x00, length_boundary));
    auto v5_truncated = build_v5_blob(0, 0x00, v5_one);
    v5_truncated.resize(v5_truncated.size() - 1);
    write_blob(dir + "/v5_truncated_frame.traceon", v5_truncated);

    // ── Valid migration seeds in the real legacy v4 layouts ───────────────
    std::vector<char> payload;
    write_genome_payload(payload, {{"seq1", "ACGT", ""}});
    std::vector<char> v4_dup;
    write_genome_payload(v4_dup, {{"dup", "ACGT", ""}, {"dup", "TGCA", ""}});
    write_blob(dir + "/v4_duplicate_ids.traceon",
               build_v4_blob(0, v4_dup, lz4f_compress(v4_dup)));
    std::vector<char> v4_ngs_dup;
    append_le64(v4_ngs_dup, 2);
    for (uint64_t legacy_hash : {uint64_t{1}, uint64_t{2}}) {
        append_le64(v4_ngs_dup, legacy_hash);
        for (const std::string& field : {std::string("dup"), std::string("ACGT"), std::string("IIII")}) {
            append_le32(v4_ngs_dup, static_cast<uint32_t>(field.size()));
            v4_ngs_dup.insert(v4_ngs_dup.end(), field.begin(), field.end());
        }
    }
    write_blob(dir + "/v4_duplicate_ids_ngs.traceon",
               build_v4_blob(1, v4_ngs_dup, lz4f_compress(v4_ngs_dup)));
    auto frame = lz4f_compress(payload);
    write_blob(dir + "/v4_valid_genome.traceon", build_v4_blob(0, payload, frame));
    std::vector<char> ngs_payload;
    const uint64_t one = 1;
    ngs_payload.insert(ngs_payload.end(), reinterpret_cast<const char*>(&one),
                       reinterpret_cast<const char*>(&one) + sizeof(one));
    const uint64_t legacy_hash = 0x123456789abcdef0ULL;
    ngs_payload.insert(ngs_payload.end(), reinterpret_cast<const char*>(&legacy_hash),
                       reinterpret_cast<const char*>(&legacy_hash) + sizeof(legacy_hash));
    for (const std::string field : {std::string("read1"), std::string("ACGT"), std::string("IIII")}) {
        const uint32_t len = static_cast<uint32_t>(field.size());
        ngs_payload.insert(ngs_payload.end(), reinterpret_cast<const char*>(&len),
                           reinterpret_cast<const char*>(&len) + sizeof(len));
        ngs_payload.insert(ngs_payload.end(), field.begin(), field.end());
    }
    write_blob(dir + "/v4_valid_ngs.traceon",
               build_v4_blob(1, ngs_payload, lz4f_compress(ngs_payload)));

    // ── Crafted malformed variants (mirroring tests/test_binary_cache.cpp) ──

    // Legacy v1/v2/v3 magic.
    for (uint8_t ver : {0x01, 0x02, 0x03}) {
        const std::string name = dir + "/v4_legacy_v" + std::to_string(ver) + ".traceon";
        write_blob(name, build_v4_blob(0, payload, frame, /*magic_ver=*/ver));
    }
    // Unknown future version and non-TRO magic.
    write_blob(dir + "/v4_magic_ver9.traceon",
               build_v4_blob(0, payload, frame, /*magic_ver=*/0x09));
    {
        auto blob = build_v4_blob(0, payload, frame);
        blob[0] = 'X'; blob[1] = 'X'; blob[2] = 'X'; blob[3] = 'X';
        write_blob(dir + "/v4_magic_xxxx.traceon", blob);
    }
    // Unsupported codec flags and invalid index mode.
    for (uint8_t flags : {0x00, 0x02, 0xFF}) {
        write_blob(dir + "/v4_bad_flags_" + std::to_string(static_cast<int>(flags)) + ".traceon",
                   build_v4_blob(0, payload, frame, std::nullopt, flags));
    }
    write_blob(dir + "/v4_bad_mode2.traceon", build_v4_blob(2, payload, frame));

    // Old v1.0 MMAP format.
    {
        std::vector<char> old{'M', 'M', 'A', 'P', 0, 0, 0, 0, 0, 0, 0, 0, 0};
        write_blob(dir + "/v4_old_mmap.traceon", old);
    }

    // Truncations of a valid v5 file: header end (27), a few bytes into the
    // frame, mid-frame, and one byte short of the end.
    {
        TracEon::SmartStrategy s;
        s.addEntry("a", "ACGT", "");
        s.addEntry("b", "TTTT", "");
        s.saveBinary(dir + "/_tmp_valid.traceon");
        std::ifstream in(dir + "/_tmp_valid.traceon", std::ios::binary | std::ios::ate);
        const size_t sz = static_cast<size_t>(in.tellg());
        in.seekg(0);
        std::vector<char> bytes(sz);
        in.read(bytes.data(), static_cast<std::streamsize>(sz));
        in.close();
        std::remove((dir + "/_tmp_valid.traceon").c_str());
        const size_t cuts[] = {27, 31, sz / 2, sz - 1};
        for (size_t i = 0; i < 4; ++i)
            write_blob(dir + "/v5_trunc_" + std::to_string(i) + ".traceon",
                       std::vector<char>(bytes.begin(),
                                         bytes.begin() + static_cast<long>(cuts[i])));
    }

    // Absurd / wrong logical lengths and frame length.
    write_blob(dir + "/v4_logical_max.traceon",
               build_v4_blob(0, payload, frame, std::nullopt, std::nullopt,
                             /*logical_len=*/UINT64_MAX));
    write_blob(dir + "/v4_logical_small.traceon",
               build_v4_blob(0, payload, frame, std::nullopt, std::nullopt,
                             /*logical_len=*/4));
    write_blob(dir + "/v4_frame_max.traceon",
               build_v4_blob(0, payload, frame, std::nullopt, std::nullopt,
                             std::nullopt, /*frame_len=*/UINT64_MAX));
    write_blob(dir + "/v4_empty_frame.traceon",
               build_v4_blob(0, payload, frame, std::nullopt, std::nullopt,
                             std::nullopt, /*frame_len=*/0));

    // Corrupted checksum field (offset 22) and corrupted header byte (offset 5).
    {
        auto blob = build_v4_blob(0, payload, frame);
        blob[22] ^= 0xFF;
        write_blob(dir + "/v4_crc_mutated.traceon", blob);
    }
    {
        auto blob = build_v4_blob(0, payload, frame);
        blob[5] ^= 0x01;
        write_blob(dir + "/v4_header_mutated.traceon", blob);
    }

    // Count bombs (payload count field patched; CRC kept valid so the parse
    // step is reached): count=5e8 with a tiny payload (vuln-0009 port),
    // count=3 with one record, count=2 with one record + 8 junk bytes.
    {
        std::vector<char> p;
        uint64_t count = 500'000'000ULL;
        for (int i = 0; i < 8; ++i) p.push_back(static_cast<char>((count >> (8 * i)) & 0xFFu));
        uint32_t len = 1;
        p.insert(p.end(), reinterpret_cast<const char*>(&len), reinterpret_cast<const char*>(&len) + 4);
        p.push_back('a');
        len = 4;
        p.insert(p.end(), reinterpret_cast<const char*>(&len), reinterpret_cast<const char*>(&len) + 4);
        const char seq[] = "ACGT";
        p.insert(p.end(), seq, seq + 4);
        len = 0;
        p.insert(p.end(), reinterpret_cast<const char*>(&len), reinterpret_cast<const char*>(&len) + 4);
        write_blob(dir + "/v4_count_bomb_5e8.traceon",
                   build_v4_blob(0, p, lz4f_compress(p)));
    }
    {
        std::vector<char> p;
        write_genome_payload(p, {{"seq1", "ACGT", ""}});
        uint64_t count = 3;
        for (int i = 0; i < 8; ++i) p[static_cast<size_t>(i)] = static_cast<char>((count >> (8 * i)) & 0xFFu);
        write_blob(dir + "/v4_bad_count_3.traceon", build_v4_blob(0, p, lz4f_compress(p)));
    }
    {
        std::vector<char> p;
        write_genome_payload(p, {{"seq1", "ACGT", ""}});
        uint64_t count = 2;
        for (int i = 0; i < 8; ++i) p[static_cast<size_t>(i)] = static_cast<char>((count >> (8 * i)) & 0xFFu);
        p.insert(p.end(), 8, static_cast<char>(0xFF));
        write_blob(dir + "/v4_bad_count_2_junk.traceon", build_v4_blob(0, p, lz4f_compress(p)));
    }

    std::cout << "v4 migration + v5 corpus written under " << dir << std::endl;
    return 0;
}
