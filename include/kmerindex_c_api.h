#ifndef TRACEON_KMERINDEX_C_API_H
#define TRACEON_KMERINDEX_C_API_H

/* Pure-C shim exposing a TracEon-backed hash table (C++) to C code, so a
 * patched minimap2 (compiled as C) can use it as a drop-in replacement for
 * its per-bucket khash_t(idx) minimizer table (index.c).
 *
 * This mirrors minimap2's EXACT khash instantiation semantics, not a
 * generic uint64->uint64 map:
 *   #define idx_hash(a) ((a)>>1)
 *   #define idx_eq(a, b) ((a)>>1 == (b)>>1)
 *   KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
 * i.e. the key's bit 0 is a "single-occurrence" flag that both hash and
 * equality ignore, so a lookup with bit0=0 finds an entry that was stored
 * with bit0=1. index.c never inserts a duplicate base-key (every kh_put
 * there is asserted absent==true, since the array being indexed is
 * pre-sorted and grouped by distinct minimizer) -- so this shim only
 * needs a one-shot "insert with final key" call, never khash's in-place
 * kh_key(...) |= 1 post-insert mutation (which isn't safely expressible
 * over a map that doesn't guarantee stable in-place key mutation).
 *
 * ---------------------------------------------------------------------------
 * ERROR MODEL (C ABI boundary — no C++ exceptions ever cross this API)
 * ---------------------------------------------------------------------------
 * Every exported function is noexcept at the ABI boundary: all C++
 * exceptions (std::bad_alloc from reserve/try_emplace, etc.) are caught
 * and translated to status codes, so a C caller can never observe a
 * terminate/abort from this library.
 *
 *   * kmerindex_create()  returns NULL on allocation failure.
 *   * kmerindex_reserve() returns 1 on success, 0 on failure.
 *   * kmerindex_insert()  returns  1 = inserted, 0 = not inserted
 *                         (collision, frozen index, or null handle),
 *                         -1 = a C++ exception was caught (e.g. bad_alloc).
 *   * Pointer-returning functions return NULL on error/miss.
 *   * NULL handle / NULL output-pointer arguments never crash: they return
 *     the function's failure value and record a diagnostic.
 *
 * A thread-local human-readable diagnostic for the LAST failed call in the
 * calling thread is available via kmerindex_last_error() (never NULL;
 * empty string when the last call succeeded). Each C thread has its own
 * error slot, so concurrent callers never race on it.
 *
 * ---------------------------------------------------------------------------
 * FREEZE / POINTER-STABILITY CONTRACT
 * ---------------------------------------------------------------------------
 * kmerindex_get() returns an interior pointer into the index. To prevent
 * dangling-pointer hazards from a later rehash, the index is FROZEN once
 * any of the following happens (whichever comes first):
 *   * kmerindex_get() is called,
 *   * kmerindex_iter_begin() is called, or
 *   * kmerindex_freeze() is called explicitly.
 * After freezing, kmerindex_insert() and kmerindex_reserve() fail with
 * return 0 and last_error "index is frozen". Read operations (get, size,
 * iteration) remain valid and safe after freezing.
 *
 * ---------------------------------------------------------------------------
 * ITERATION (caller-owned, reentrant)
 * ---------------------------------------------------------------------------
 * kmerindex_iter_t is an opaque, FIXED-SIZE structure owned by the caller
 * (stack-allocatable in C). Each iter_begin/iter_next pair uses its own
 * private cursor -- no shared mutable state lives in the index, so two
 * iterators over the same index can be interleaved freely, and iteration
 * is safe for concurrent readers of a frozen index. An iterator must be
 * initialized with kmerindex_iter_begin() before kmerindex_iter_next().
 * iter_begin() freezes the index (see above).
 */

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
/* noexcept at the ABI boundary: C++ exceptions never cross it (vuln-0003).
 * C callers see the plain (non-throwing) declaration. */
#define TRACEON_CAPI_NOEXCEPT noexcept
#else
#define TRACEON_CAPI_NOEXCEPT
#endif

typedef struct kmerindex_t kmerindex_t;

/* Caller-owned opaque iterator. The fixed size (4 pointers) is a public
 * ABI commitment: the internal C++ cursor must fit inside it (checked by
 * static_assert in the implementation). */
typedef struct kmerindex_iter_t {
    void* opaque[4];
} kmerindex_iter_t;

/* Creates a new empty index. Returns NULL (and records an error) on
 * allocation failure. */
kmerindex_t* kmerindex_create(void) TRACEON_CAPI_NOEXCEPT;

/* Destroys the index and frees all resources. NULL-safe. */
void kmerindex_destroy(kmerindex_t* h) TRACEON_CAPI_NOEXCEPT;

/* Pre-reserves capacity for `n` entries (avoids mid-build rehash). Returns
 * 1 on success, 0 on failure (allocation failure, frozen index, or NULL
 * handle). */
int kmerindex_reserve(kmerindex_t* h, size_t n) TRACEON_CAPI_NOEXCEPT;

/* Number of stored entries. Returns 0 on a NULL handle. */
size_t kmerindex_size(const kmerindex_t* h) TRACEON_CAPI_NOEXCEPT;

/* Inserts `final_key` (already carrying the single-occurrence flag bit,
 * i.e. base_key | (is_singleton ? 1 : 0)) -> value. Returns:
 *    1  inserted
 *    0  NOT inserted -- base_key collision (mirrors `assert(absent)` in
 *       the original code; should never happen given index.c's
 *       pre-grouped-by-distinct-minimizer insert order), OR the index is
 *       frozen, OR the handle is NULL. Distinguish via
 *       kmerindex_last_error() (empty string = plain collision).
 *   -1  a C++ exception was caught (e.g. std::bad_alloc on rehash).
 * Fails (0, "index is frozen") once kmerindex_get()/iter_begin()/freeze()
 * has begun the read phase. */
int kmerindex_insert(kmerindex_t* h, uint64_t final_key, uint64_t value) TRACEON_CAPI_NOEXCEPT;

/* Looks up `query_key` (bit0 always 0, as index.c constructs it fresh per
 * query) using minimap2's ignore-bit0 equality. On hit, returns a pointer
 * to the STORED VALUE (stable for the lifetime of the index -- the index
 * is frozen by this call, matching mm_idx_get's contract of returning
 * `&kh_val(h,k)` directly), and writes the ACTUAL stored key (whose bit0
 * tells single- vs multi-occurrence, matching `kh_key(h,k)&1` in
 * mm_idx_get) to *matched_key_out. Returns NULL on miss or on error
 * (NULL handle or NULL matched_key_out). FREEZES the index. */
const uint64_t* kmerindex_get(const kmerindex_t* h, uint64_t query_key,
                               uint64_t* matched_key_out) TRACEON_CAPI_NOEXCEPT;

/* Explicitly freezes the index: subsequent insert/reserve fail. Returns 1
 * on success, 0 on a NULL handle. Read operations remain available. */
int kmerindex_freeze(kmerindex_t* h) TRACEON_CAPI_NOEXCEPT;

/* Initializes the caller-owned iterator `out` over the index `h`. Returns
 * 1 on success, 0 on error (NULL handle/iterator). FREEZES the index. */
int kmerindex_iter_begin(const kmerindex_t* h, kmerindex_iter_t* out) TRACEON_CAPI_NOEXCEPT;

/* Advances the caller-owned iterator; yields one (stored_key, value) pair
 * per call, in no particular order. Returns 1 while entries remain, 0 at
 * end or on error (NULL iterator/outputs). Reentrant: independent cursors
 * per kmerindex_iter_t. */
int kmerindex_iter_next(kmerindex_iter_t* it, uint64_t* key_out, uint64_t* val_out) TRACEON_CAPI_NOEXCEPT;

/* Thread-local diagnostic for the last failed call in the calling thread.
 * Never NULL; empty string if the last call succeeded. Passing a NULL
 * handle is safe and still returns the diagnostic. */
const char* kmerindex_last_error(const kmerindex_t* h) TRACEON_CAPI_NOEXCEPT;

#ifdef __cplusplus
}
#endif

#endif /* TRACEON_KMERINDEX_C_API_H */
