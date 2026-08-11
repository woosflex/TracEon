// fuzz_trki.cpp — libFuzzer target for the TRKI mmap-persistent k-mer
// reference index loader (KmerReferenceIndex::mmap_open).
//
// Contract under test (vuln-0006): every file-controlled field (magic,
// version, k, n_slots, n_keys, n_positions, slot offset/count pairs) is
// validated with overflow-safe arithmetic before any pointer/product is
// computed; malformed/crafted files must be rejected with
// std::runtime_error — never an OOB read, never a wild pointer, never a
// SEGV. For the (astronomically unlikely) file that passes validation, the
// harness additionally exercises query() with derived k-mers — query()
// re-validates each slot's (offset, count) before building the return span.

#include "KmerReferenceIndex.h"
#include "fuzz_common.h"

#include <cstdint>
#include <string>
#include <vector>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    const std::string path = traceon_fuzz::write_input_to_tmp(data, size);
    try {
        auto idx = TracEon::KmerReferenceIndex::mmap_open(path);
        uint64_t base_kmer = 0;
        if (size >= 8) {
            for (int i = 0; i < 8; ++i)
                base_kmer = (base_kmer << 8) | data[static_cast<size_t>(i)];
        }
        volatile size_t sink = 0;
        for (int i = 0; i < 16; ++i) {
            const auto span = idx.query(base_kmer + static_cast<uint64_t>(i) * 2654435761ull);
            sink += span.size();
            for (uint64_t p : span) sink += static_cast<size_t>(p & 0xFF);
        }
        sink += idx.n_keys() + static_cast<size_t>(idx.k());
        (void)sink;
    } catch (const std::exception&) {
        // Expected: every malformed TRKI must be rejected cleanly.
    }
    return 0;
}
