#ifndef TRACEON_KMER_ENCODING_H
#define TRACEON_KMER_ENCODING_H

#include <cstdint>
#include <string_view>

namespace TracEon {

// 2-bit-per-base packing for k-mers up to k=32 (fits in a uint64_t).
// Matches the representation minimap2/khash-based tools use internally.
inline int kmer_base_code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

// Encodes the k-mer at seq[0..k) into a packed uint64_t. Returns false
// (without touching `out`) if:
//   - k is outside [1, 32] (32 bases x 2 bits is the uint64_t packing limit),
//   - the view is shorter than k (never reads past seq.size(),
//     see vuln-0001), or
//   - the window contains a non-ACGT base.
// On success `out` receives the packed 2-bit-per-base code.
inline bool encode_kmer(std::string_view seq, int k, uint64_t& out) {
    if (k <= 0 || k > 32 || (int)seq.size() < k) return false;
    uint64_t code = 0;
    for (int i = 0; i < k; ++i) {
        int b = kmer_base_code(seq[i]);
        if (b < 0) return false;
        code = (code << 2) | (uint64_t)b;
    }
    out = code;
    return true;
}

} // namespace TracEon

#endif // TRACEON_KMER_ENCODING_H
