// gen_v4_corpus.cpp — regenerates the `.traceon` v4 seed corpus for
// fuzz_v4_loader: a valid GENOME + a valid NGS round-trip cache plus every
// malformed variant from tests/test_binary_cache.cpp (legacy magic, bad
// flags/mode, truncations, wrong lengths, mutated CRC, count bombs).
//
// The valid files go through the real SmartStrategy::saveBinary() path; the
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

} // namespace

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : "fuzz/corpus/v4_loader";

    // ── Valid round-trip caches via the real saveBinary() path ──────────────
    {
        TracEon::SmartStrategy s;
        s.addEntry("chr1", "ACGTACGTACGTACGT", "");
        s.addEntry("chr2", "TTTTGGGGCCCCAAAA", "");
        s.saveBinary(dir + "/v4_valid_genome.traceon");
    }
    {
        TracEon::SmartStrategy s(TracEon::IndexMode::NGS);
        s.addEntry("read1", "ACGT", "IIII");
        s.addEntry("read2", "TGCA", "####");
        s.saveBinary(dir + "/v4_valid_ngs.traceon");
    }

    // ── Crafted malformed variants (mirroring tests/test_binary_cache.cpp) ──
    std::vector<char> payload;
    write_genome_payload(payload, {{"seq1", "ACGT", ""}});
    const auto frame = lz4f_compress(payload);

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

    // Truncations of a valid file: header end (26), a few bytes into the
    // frame, mid-frame, one byte short of the end.
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
        const size_t cuts[] = {26, 30, sz / 2, sz - 1};
        for (size_t i = 0; i < 4; ++i)
            write_blob(dir + "/v4_trunc_" + std::to_string(i) + ".traceon",
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

    std::cout << "v4 corpus written under " << dir << std::endl;
    return 0;
}
