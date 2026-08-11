#ifndef TRACEON_KMER_INDEX_H
#define TRACEON_KMER_INDEX_H

#include "MapDefs.h"
#include "KmerEncoding.h"

#include <cstdint>
#include <string_view>
#include <vector>

namespace TracEon {

// Monomorphic packed-k-mer -> position index. Deliberately independent of
// SmartStrategy/file_cache_: no std::variant discriminant, no string keys,
// no SequenceView, no id re-compare -- the packed uint64_t key IS the
// k-mer, so a hash-bucket hit is a confirmed match (subject to the usual
// 2^-64 collision odds of any 64-bit hash key, same as khash/minimap2).
class KmerIndex {
public:
    using iterator = typename HashMap<uint64_t, uint64_t>::iterator;
    using const_iterator = typename HashMap<uint64_t, uint64_t>::const_iterator;

    void reserve(size_t n) { map_.reserve(n); }

    void insert(uint64_t kmer, uint64_t pos) { map_[kmer] = pos; }

    // Returns nullptr on miss; otherwise a pointer to the stored position.
    const uint64_t* find(uint64_t kmer) const {
        auto it = map_.find(kmer);
        return it == map_.end() ? nullptr : &it->second;
    }

    // khash kh_put()-equivalent: get-or-create the slot for `kmer`, report
    // via `is_absent` whether it was freshly inserted (value-initialized to
    // 0, caller fills it in) or already present (existing value returned
    // untouched) -- lets C callers (e.g. a patched minimap2) mirror khash's
    // put-then-fill call pattern almost mechanically.
    uint64_t& put(uint64_t kmer, bool& is_absent) {
        auto [it, inserted] = map_.try_emplace(kmer, 0);
        is_absent = inserted;
        return it->second;
    }

    size_t size() const { return map_.size(); }

    iterator begin() { return map_.begin(); }
    iterator end() { return map_.end(); }
    const_iterator begin() const { return map_.begin(); }
    const_iterator end() const { return map_.end(); }

    // Builds the index over every valid k-mer window in `seq`, using a
    // rolling 2-bit encode (shift+mask) instead of re-encoding all k bases
    // per position -- O(1) amortized per position instead of O(k).
    void build_from_sequence(std::string_view seq, int k) {
        if (k <= 0 || (int)seq.size() < k) return;  // k<=0: negative shift UB; short seq: no window
        reserve(seq.size() - k + 1);

        const uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
        uint64_t code = 0;
        int valid_run = 0; // consecutive valid bases accumulated into `code`

        for (size_t i = 0; i < seq.size(); ++i) {
            int b = kmer_base_code(seq[i]);
            if (b < 0) {
                valid_run = 0;
                code = 0;
                continue;
            }
            code = ((code << 2) | (uint64_t)b) & mask;
            valid_run++;
            if (valid_run >= k) {
                insert(code, i - k + 1);
            }
        }
    }

private:
    HashMap<uint64_t, uint64_t> map_;
};

} // namespace TracEon

#endif // TRACEON_KMER_INDEX_H
