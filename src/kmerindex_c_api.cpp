#include "kmerindex_c_api.h"

#include <ankerl/unordered_dense.h>

#include <cstdint>
#include <cstddef>
#include <new>
#include <string>
#include <type_traits>

// Custom hash/equality replicating minimap2's exact khash instantiation
// for its minimizer table (index.c: idx_hash/idx_eq) -- bit 0 of the key
// is a flag both functions ignore, so a query with bit0=0 correctly finds
// an entry stored with bit0=1. This is deliberately NOT the same as
// TracEon::KmerIndex (Phase 1's validated generic packed-kmer index,
// include/KmerIndex.h), which uses plain full-key equality -- reusing
// that class here would silently produce wrong lookups on any minimizer
// that was ever flagged as a singleton.
struct MinimizerKeyHash {
    using is_avalanching = void;
    size_t operator()(uint64_t k) const noexcept {
        return ankerl::unordered_dense::hash<uint64_t>{}(k >> 1);
    }
};
struct MinimizerKeyEq {
    bool operator()(uint64_t a, uint64_t b) const noexcept { return (a >> 1) == (b >> 1); }
};

using MinimizerMap = ankerl::unordered_dense::map<uint64_t, uint64_t, MinimizerKeyHash, MinimizerKeyEq>;

struct kmerindex_t {
    MinimizerMap map;
    // Set once the read phase begins (get/iter_begin/freeze). Blocks
    // insert/reserve so interior pointers returned by get() can never be
    // invalidated by a later rehash. `mutable` because get()/iter_begin()
    // take a const handle (reads) but still transition the index phase.
    mutable bool frozen = false;
};

// Internal state stored inside the caller-owned opaque kmerindex_iter_t.
// `magic` guards against iter_next() being called on an uninitialized
// iterator; the ankerl const_iterator over a map is a pointer into the
// map's value vector, so it stays valid as long as the (frozen) index is
// alive and unmodified -- the caller must not destroy the index while an
// iterator over it is live.
struct KmerIndexIterState {
    static constexpr uint64_t kMagic = 0x4B4D455249544552ull; // "KMERITER"
    uint64_t magic = 0;
    const MinimizerMap* map = nullptr;
    MinimizerMap::const_iterator it{};
};
static_assert(sizeof(KmerIndexIterState) <= sizeof(kmerindex_iter_t),
              "kmerindex_iter_t opaque storage too small for the internal iterator");

namespace {

// Thread-local diagnostic slot: each C thread sees its own last-error
// string, so concurrent callers never race on it. Never left dangling:
// set on every failure, cleared on every success.
thread_local std::string g_last_error;

void set_error(const char* msg) noexcept { g_last_error = msg ? msg : ""; }
void clear_error() noexcept { g_last_error.clear(); }

} // namespace

extern "C" {

kmerindex_t* kmerindex_create(void) noexcept {
    try {
        return new kmerindex_t();
    } catch (const std::bad_alloc&) {
        set_error("kmerindex_create: allocation failed");
        return nullptr;
    } catch (...) {
        set_error("kmerindex_create: unknown error");
        return nullptr;
    }
}

void kmerindex_destroy(kmerindex_t* h) noexcept {
    // Destructors here are noexcept (map/vector teardown cannot throw);
    // delete of nullptr is a no-op.
    delete h;
}

int kmerindex_reserve(kmerindex_t* h, size_t n) noexcept {
    if (!h) {
        set_error("kmerindex_reserve: null handle");
        return 0;
    }
    if (h->frozen) {
        set_error("kmerindex_reserve: index is frozen");
        return 0;
    }
    try {
        h->map.reserve(n);
        clear_error();
        return 1;
    } catch (const std::length_error&) {
        set_error("kmerindex_reserve: requested size too large");
        return 0;
    } catch (const std::bad_alloc&) {
        set_error("kmerindex_reserve: allocation failed");
        return 0;
    } catch (...) {
        set_error("kmerindex_reserve: unknown error");
        return 0;
    }
}

size_t kmerindex_size(const kmerindex_t* h) noexcept {
    if (!h) {
        set_error("kmerindex_size: null handle");
        return 0;
    }
    clear_error();
    return h->map.size();
}

int kmerindex_insert(kmerindex_t* h, uint64_t final_key, uint64_t value) noexcept {
    if (!h) {
        set_error("kmerindex_insert: null handle");
        return 0;
    }
    if (h->frozen) {
        set_error("kmerindex_insert: index is frozen");
        return 0;
    }
    try {
        auto [it, inserted] = h->map.try_emplace(final_key, value);
        if (!inserted) {
            // Base-key collision: a defined (should-never-happen) outcome,
            // not an error -- clear last_error so callers can tell it
            // apart from frozen/error returns (empty = plain collision).
            clear_error();
            return 0;
        }
        clear_error();
        return 1;
    } catch (const std::bad_alloc&) {
        set_error("kmerindex_insert: allocation failed");
        return -1;
    } catch (...) {
        set_error("kmerindex_insert: unknown error");
        return -1;
    }
}

const uint64_t* kmerindex_get(const kmerindex_t* h, uint64_t query_key,
                               uint64_t* matched_key_out) noexcept {
    if (!h) {
        set_error("kmerindex_get: null handle");
        return nullptr;
    }
    if (!matched_key_out) {
        set_error("kmerindex_get: null matched_key_out");
        return nullptr;
    }
    try {
        // Freeze BEFORE lookup: the returned pointer is an interior pointer
        // that must remain stable; no inserts may happen from this point on.
        h->frozen = true;
        auto it = h->map.find(query_key);
        if (it == h->map.end()) {
            clear_error();
            return nullptr;
        }
        *matched_key_out = it->first;
        clear_error();
        return &it->second;
    } catch (...) {
        set_error("kmerindex_get: unknown error");
        return nullptr;
    }
}

int kmerindex_freeze(kmerindex_t* h) noexcept {
    if (!h) {
        set_error("kmerindex_freeze: null handle");
        return 0;
    }
    h->frozen = true;
    clear_error();
    return 1;
}

int kmerindex_iter_begin(const kmerindex_t* h, kmerindex_iter_t* out) noexcept {
    if (!h) {
        set_error("kmerindex_iter_begin: null handle");
        return 0;
    }
    if (!out) {
        set_error("kmerindex_iter_begin: null iterator");
        return 0;
    }
    try {
        auto* st = new (out->opaque) KmerIndexIterState();
        st->magic = KmerIndexIterState::kMagic;
        st->map = &h->map;
        st->it = h->map.begin();
        h->frozen = true; // iteration is part of the read phase
        clear_error();
        return 1;
    } catch (...) {
        set_error("kmerindex_iter_begin: unknown error");
        return 0;
    }
}

int kmerindex_iter_next(kmerindex_iter_t* it, uint64_t* key_out, uint64_t* val_out) noexcept {
    if (!it) {
        set_error("kmerindex_iter_next: null iterator");
        return 0;
    }
    if (!key_out || !val_out) {
        set_error("kmerindex_iter_next: null output pointer");
        return 0;
    }
    auto* st = reinterpret_cast<KmerIndexIterState*>(it->opaque);
    if (st->magic != KmerIndexIterState::kMagic) {
        set_error("kmerindex_iter_next: iterator not initialized (call kmerindex_iter_begin first)");
        return 0;
    }
    if (!st->map || st->it == st->map->end()) {
        clear_error();
        return 0;
    }
    *key_out = st->it->first;
    *val_out = st->it->second;
    ++st->it;
    clear_error();
    return 1;
}

const char* kmerindex_last_error(const kmerindex_t* h) noexcept {
    (void)h; // diagnostic is thread-local; handle presence is irrelevant
    return g_last_error.c_str();
}

} // extern "C"
