# ADR-001: Lock-Free Read Optimization

**Date:** December 2025  
**Applies to:** TracEon v1.0.0 "Avalon"  
**Deciders:** Adnan Raza

## Context

TracEon's initial implementation used `std::shared_mutex` for all cache access operations. While safe, profiling revealed that lock overhead consumed 25-35% of lookup time on large datasets (100MB+), capping throughput at ~9M OPS/s.

## Decision

We implement a **lock-free read path** using C++20 atomics with acquire/release semantics:

1. After `loadFile()` or `loadBinary()` completes, all data structures become **immutable**
2. An `std::atomic<bool> data_ready_` flag signals when reads are safe
3. Read operations (`getView()`, `hasSequence()`) bypass the mutex once `data_ready_` is true
4. Write operations (`loadFile()`, `clearCache()`) continue using exclusive locks

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
- Static assertion that `file_cache_` is never modified after `data_ready_` becomes true

## Alternatives Considered

1. **RW Lock Upgrade**: Still requires atomic operations; minimal improvement
2. **RCU (Read-Copy-Update)**: Overly complex for immutable-after-load semantics
3. **Thread-Local Caching**: Memory overhead unacceptable for genomic workloads

## References
- [C++ Memory Model](https://en.cppreference.com/w/cpp/atomic/memory_order)
- Preshing, *Acquire and Release Semantics* (2012)
- Boehm & Adve, *Foundations of the C++ Concurrency Memory Model* (2008)