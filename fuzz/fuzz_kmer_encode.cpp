// fuzz_kmer_encode.cpp — libFuzzer target for encode_kmer() (the vuln-0001
// contract: never read past seq.size(), never shift by k<=0 or k>32).
//
// Input layout: the first 4 bytes are a little-endian int32 candidate k;
// the remaining bytes are the sequence view. The harness additionally probes
// the guard boundaries (0, 1, 31, 32, 33, -1, INT32_MAX, and a k derived
// from the first sequence byte) so a single input exercises both the raw
// attacker-controlled k and the edges. Any OOB read / UB here is a finding.

#include "KmerEncoding.h"
#include "fuzz_common.h"

#include <cstdint>
#include <string_view>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    if (size == 0) return 0;

    int k = 0;
    size_t seq_off = 0;
    if (size >= 4) {
        k = static_cast<int>(static_cast<uint32_t>(data[0]) |
                             (static_cast<uint32_t>(data[1]) << 8) |
                             (static_cast<uint32_t>(data[2]) << 16) |
                             (static_cast<uint32_t>(data[3]) << 24));
        seq_off = 4;
    } else {
        k = static_cast<int>(data[0]);
        seq_off = 1;
    }
    const std::string_view seq(reinterpret_cast<const char*>(data + seq_off),
                               size - seq_off);

    const int candidates[] = {
        k,
        0,          // k == 0: must return false, never shift
        1,          // minimal valid k
        31, 32,     // packing limits
        33,         // one past the uint64_t packing limit
        -1,         // negative k
        2147483647, // INT32_MAX: absurd k
        seq.empty() ? 0 : static_cast<int>(static_cast<unsigned char>(seq[0])),
    };

    uint64_t out = 0;
    volatile uint64_t sink = 0;
    for (int c : candidates) {
        const bool ok = TracEon::encode_kmer(seq, c, out);
        sink += ok ? out : static_cast<uint64_t>(static_cast<unsigned int>(c));
    }
    (void)sink;
    return 0;
}
