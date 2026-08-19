#include "kmerindex_c_api.h"

#include <ankerl/unordered_dense.h>

#include <atomic>
#include <cstdint>
#include <cstddef>
#include <new>
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
    // Phase-transition flag: set once the read phase begins
    // (get/iter_begin/freeze), blocking insert/reserve so interior
    // pointers returned by get() can never be invalidated by a later
    // rehash. ATOMIC (not a plain bool): minimap2 performs concurrent
    // lookups from mapping worker threads, and every get() WRITES the
    // flag while a straggler insert()/reserve() READS it -- concurrent
    // read/write on a plain bool is UB. The freeze paths use a RELEASE
    // store (publishing all map writes done before the freeze) and the
    // check paths an ACQUIRE load, so an insert that observes frozen==true
    // happens-after every write published by the freeze. `mutable` because
    // get()/iter_begin() take a const handle (reads) but still transition
    // the index phase.
    mutable std::atomic<bool> frozen{false};
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

// Fixed-size, NON-ALLOCATING thread-local diagnostic slot. Each C thread
// sees its own last-error message, so concurrent callers never race on it.
// The whole point of the fixed buffer is the noexcept boundary: assigning
// a std::string here (as before) can allocate, and an allocation that
// throws std::bad_alloc from within a noexcept function calls
// std::terminate -- so the advertised no-throw boundary was not total.
// Every error write goes through set_error(), which copies into this
// buffer with bounds checks and never allocates, and is therefore safe
// to call from catch blocks even under OOM. kmerindex_last_error()
// returns a pointer into this buffer; the content is valid until the
// next error (or success, which clears it) on the SAME thread.
static thread_local char g_last_error[256];

// Bounded non-allocating copy: never allocates, never leaves the buffer
// unterminated, truncates (with NUL) if the message does not fit.
void set_error(const char* msg) noexcept {
    size_t i = 0;
    const size_t cap = sizeof(g_last_error) - 1;
    if (msg) {
        while (i < cap && msg[i] != '\0') {
            g_last_error[i] = msg[i];
            ++i;
        }
    }
    g_last_error[i] = '\0';
}

void clear_error() noexcept { g_last_error[0] = '\0'; }

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
    // ACQUIRE: observe a freeze published by get()/iter_begin()/freeze()
    // (their RELEASE stores). If the index is frozen, fail BEFORE touching
    // the map -- no mutation may race with the concurrent readers.
    if (h->frozen.load(std::memory_order_acquire)) {
        set_error("kmerindex_reserve: index is frozen");
        return 0;
    }
    // Deterministic pre-check (GitHub issue #9): ankerl::unordered_dense
    // 4.4.0 reserve() CLAMPS to max_size() instead of throwing, so a
    // request above max_size() (~2^32 entries for this uint32 value_idx)
    // silently succeeds on machines where the huge allocation lands -- the
    // C boundary would return 1 instead of the documented 0. Reject such
    // sizes here, BEFORE any allocation, so the contract is
    // platform-independent. The try/catch below still handles GENUINE
    // allocation failures (bad_alloc / length_error) for sizes at or
    // below max_size() under real memory pressure.
    if (n > h->map.max_size()) {
        set_error("kmerindex_reserve: requested size too large");
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
    // ACQUIRE: observe a freeze published by get()/iter_begin()/freeze()
    // (their RELEASE stores). If the index is frozen, fail BEFORE touching
    // the map -- no mutation may race with the concurrent readers.
    if (h->frozen.load(std::memory_order_acquire)) {
        set_error("kmerindex_insert: index is frozen");
        return 0;
    }
    try {
        auto [it, inserted] = h->map.try_emplace(final_key, value);
        if (!inserted) {
            // `inserted` here is the map's own try_emplace bool (0/1 only,
            // never -1), NOT the C API status code -- this branch is a
            // base-key collision, not an error. Only the RETURN below is
            // the public contract (1 inserted / 0 not / -1 exception).
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
        // RELEASE store: any insert that later ACQUIRE-loads the flag and
        // fails has happened-after this point, so the interior pointer it
        // might otherwise have invalidated can never dangle.
        h->frozen.store(true, std::memory_order_release);
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
    // RELEASE store (see kmerindex_get): publishes the freeze so any
    // concurrent/future insert/reserve ACQUIRE-load sees it.
    h->frozen.store(true, std::memory_order_release);
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
        // Iteration is part of the read phase; RELEASE store (see
        // kmerindex_get) so inserts/reserves observe the freeze.
        h->frozen.store(true, std::memory_order_release);
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
    return g_last_error;
}

} // extern "C"
