#include "kmerindex_c_api.h"

#include <ankerl/unordered_dense.h>

#include <atomic>
#include <cstdint>
#include <cstddef>
#include <new>
#include <type_traits>
#include <mutex>

// Custom hash/equality replicating minimap2's exact khash instantiation
// for its minimizer table (index.c: idx_hash/idx_eq) -- bit 0 of the key
// is a flag both functions ignore, so a query with bit0=0 correctly finds
// an entry stored with bit0=1. This is deliberately NOT plain full-key
// equality (as a generic map would use): that would silently produce wrong
// lookups on any minimizer that was ever flagged as a singleton.
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
    // Phase-transition flag: set once the read phase begins (get),
    // blocking insert/reserve so interior pointers returned by get() can
    // never be invalidated by a later rehash. ATOMIC (not a plain bool):
    // minimap2 performs concurrent lookups from mapping worker threads,
    // and every get() WRITES the flag while a straggler insert()/reserve()
    // READS it -- concurrent read/write on a plain bool is UB. The freeze
    // path uses a RELEASE store (publishing all map writes done before the
    // freeze) and the check paths an ACQUIRE load, so an insert that
    // observes frozen==true happens-after every write published by the
    // freeze. `mutable` because get() takes a const handle (a read) but
    // still transitions the index phase.
    mutable std::atomic<bool> frozen{false};
    // Serializes only the build-to-read transition and build-phase mutations.
    // Once frozen is observed true, lookups skip this mutex and the map is
    // immutable, preserving lock-free concurrent reads.
    mutable std::mutex phase_mu;
};

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
        kmerindex_t* h = new kmerindex_t();
        clear_error();
        return h;
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
    clear_error();
}

int kmerindex_reserve(kmerindex_t* h, size_t n) noexcept {
    if (!h) {
        set_error("kmerindex_reserve: null handle");
        return 0;
    }
    try {
        // Admission and mutation share the transition mutex with the first
        // get(). The lock itself is inside the noexcept translation boundary.
        std::lock_guard<std::mutex> phase_lock(h->phase_mu);
        if (h->frozen.load(std::memory_order_acquire)) {
            set_error("kmerindex_reserve: index is frozen");
            return 0;
        }
        if (n > h->map.max_size()) {
            set_error("kmerindex_reserve: requested size too large");
            return 0;
        }
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
        set_error("kmerindex_reserve: synchronization or unknown failure");
        return 0;
    }
}

size_t kmerindex_size(const kmerindex_t* h) noexcept {
    if (!h) {
        set_error("kmerindex_size: null handle");
        return 0;
    }
    try {
        if (!h->frozen.load(std::memory_order_acquire)) {
            std::lock_guard<std::mutex> phase_lock(h->phase_mu);
            const size_t result = h->map.size();
            clear_error();
            return result;
        }
        clear_error();
        return h->map.size();
    } catch (...) {
        set_error("kmerindex_size: synchronization or unknown failure");
        return 0;
    }
}

int kmerindex_insert(kmerindex_t* h, uint64_t final_key, uint64_t value) noexcept {
    if (!h) {
        set_error("kmerindex_insert: null handle");
        return 0;
    }
    try {
        std::lock_guard<std::mutex> phase_lock(h->phase_mu);
        if (h->frozen.load(std::memory_order_acquire)) {
            set_error("kmerindex_insert: index is frozen");
            return 0;
        }
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
        set_error("kmerindex_insert: synchronization or unknown failure");
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
        // The first lookup closes the build phase under the same mutex used by
        // every mutator. Later lookups observe frozen with acquire semantics
        // and remain lock-free because the map can no longer change.
        if (!h->frozen.load(std::memory_order_acquire)) {
            std::lock_guard<std::mutex> phase_lock(h->phase_mu);
            h->frozen.store(true, std::memory_order_release);
        }
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

const char* kmerindex_last_error(const kmerindex_t* h) noexcept {
    (void)h; // diagnostic is thread-local; handle presence is irrelevant
    return g_last_error;
}

} // extern "C"
