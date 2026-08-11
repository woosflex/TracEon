# ADR-001: Lock-Free Read Optimization

**Date:** December 2025  
**Applies to:** TracEon v1.0.0 "Avalon"  
**Deciders:** Adnan Raza (Woosflex)

## Context

TracEon's initial implementation used `std::shared_mutex` for all cache access operations. While safe, profiling revealed that lock overhead consumed 25-35% of lookup time on large datasets (100MB+), capping throughput at ~9M OPS/s.

## Decision

We implement a **lock-free read path** using C++20 atomics with acquire/release semantics:

1. After `loadFile()` or `loadBinary()` completes, all data structures become **immutable**
2. An `std::atomic<bool> data_ready_` flag signals when reads are safe
3. Read operations (`getView()`, `hasSequence()`) bypass the mutex once `data_ready_` is true
4. Write operations (`loadFile()`, `clearCache()`) continue using exclusive locks

### Enforcing Immutability (v1.4.1 / Bug 3 fix)

Immutability is not just a convention — it is **enforced by the API**:

- A separate `std::atomic<bool> data_loaded_` flag is set **only** by the load
  paths (`parseArena()` via `loadFile()`/`loadGzipFile()`, and `loadBinary()`)
  and cleared by `clearInternal()` (`clearCache()`).
- `addEntry()` (the `Cache::set()` path) checks `data_loaded_` first and throws
  `std::logic_error` if a load has happened — so the index can never be mutated
  after it has been published to lock-free readers.
- `data_loaded_` is distinct from `data_ready_` because `addEntry()` itself sets
  `data_ready_` so readers can see manually-built entries during a build phase.
  Using `data_ready_` alone would wrongly reject the second `set()` of a manual
  build.

## Reader Quiescence (v2.0.0 — explicit external contract)

**Date:** v2.0.0 (2026)
**Decision:** Retain the documented **external-quiescence contract** for the
naked `std::string_view` API. Reads are safe concurrently with other reads
**only while the same loaded snapshot is installed**; a returned view is
**invalid once `clearCache()`, reload, or destruction begins**.

### Why not an in-library reclamation mechanism

- `getView()` acquire-loads `data_ready_` and reads `file_cache_` without
  taking `cache_mutex_`; `clearInternal()` clears the map and then releases
  `text_arena_`, `manual_store_`, and the mmap. The release/acquire flag
  publishes initialization; it does **not** wait for a reader that already
  observed `true` and is not a reclamation mechanism.
- `std::string_view` is non-owning: even a reader-side generation check
  during lookup cannot extend the lifetime of the bytes **after** the call
  returns. A real epoch/RCU design requires a grace period and would need
  the whole immutable snapshot (map + arena/mmap ownership) as the
  protected object — a separate API design (`Snapshot`/`SnapshotView`).
- An atomic `shared_ptr<const Snapshot>` would require changing the return
  type to an owner-bearing handle; C++ also does not promise
  `atomic<shared_ptr>` is lock-free.

### What IS guaranteed

1. Concurrent `getView()`/`hasSequence()`/`getQuality()` against a *loaded*
   snapshot are safe and lock-free (no mutex on the read path).
2. Writers (`loadFile`, `loadGzipFile`, `loadBinary`, `clearCache`,
   `addEntry`) are serialized by `cache_mutex_`. The mutex does **not**
   protect readers — do not market it as reader reclamation.
3. The application must achieve quiescence (stop using views from a
   snapshot) before `clearCache()`, reload, or destruction. Sequential
   clear/reload is supported and tested.
4. `TRACEON_DEBUG_LIFECYCLE` (opt-in) adds a reader counter/scope that
   warns in `clearInternal()` when a `getView()` lookup overlaps a teardown.
   It is a diagnostic for coordinated misuse, **not** synchronization, and
   cannot detect post-return use of a retained view.

**Contract:** `set()`/`addEntry()` is allowed only **before** any load (build
phase) or **after** `clearCache()`. A loaded cache is truly immutable, which is
what makes the lock-free read guarantee sound. Build-phase `set()` calls are
single-threaded by contract (no concurrent readers during a build).

### Memory Ordering Contract
```cpp
// WRITE SIDE (happens-before relationship)
load_data_into_structures();
data_ready_.store(true, std::memory_order_release); // [A]

// READ SIDE (synchronizes-with)
if (data_ready_.load(std::memory_order_acquire)) { // [B]
    // All writes before [A] are visible after [B]
    access_immutable_structures();
}
```

## Consequences

### Positive
- **2x throughput improvement** on concurrent workloads (13M → 26M OPS/s potential)
- **Zero contention** on read-heavy workloads (typical for genomics)
- Enables horizontal scaling across cores

### Negative
- Subtle correctness requirement: structures must remain truly immutable
- Cannot support "hot reload" or incremental updates without major refactoring
- Requires careful auditing of any future write operations

### Mitigation
- Comprehensive documentation of memory ordering guarantees
- Unit tests validating concurrent access patterns
- **Enforced immutable-after-load contract**: there is **no** static assertion
  that `file_cache_` is never modified after `data_ready_` becomes true — such
  an assertion cannot exist, because `addEntry()` deliberately sets
  `data_ready_` during the build phase so lock-free readers can see
  manually-built entries (that is precisely why the distinct `data_loaded_`
  flag exists). Enforcement is runtime: the separate `data_loaded_` flag, set
  only by the load paths and cleared by `clearCache()`, gates
  `addEntry()`/`set()` — once a load has happened, mutation throws
  `std::logic_error` (see "Enforcing Immutability" above).
- **Write-side exclusivity**: `clearCache()` and reload must not overlap
  lock-free readers — they tear down `text_arena_` / the map while readers hold
  views into them. Write operations are exclusive by contract (write-side
  only); only reads of a *loaded* cache are lock-free and safe concurrently.

## Alternatives Considered

1. **RW Lock Upgrade**: Still requires atomic operations; minimal improvement
2. **RCU (Read-Copy-Update)**: Overly complex for immutable-after-load semantics
3. **Thread-Local Caching**: Memory overhead unacceptable for genomic workloads

## References
- [C++ Memory Model](https://en.cppreference.com/w/cpp/atomic/memory_order)
- Preshing, *Acquire and Release Semantics* (2012)
- Boehm & Adve, *Foundations of the C++ Concurrency Memory Model* (2008)