# ADR-003: SIMD Parsing & Hash Map Optimization

**Status:** Accepted (v1.2.0 "Caladbolg")  
**Date:** 2026-06-09  
**Deciders:** Adnan Raza (Woosflex)  
**Related:** [ADR-001: Lock-Free Reads](ADR-001-lock-free-reads.md), [ADR-002: GZIP Integration](ADR-002-gzip-integration.md)

---

## Context and Problem Statement

After v1.1.0 optimized GZIP decompression, profiling revealed the load pipeline bottleneck shifted to CPU-bound work:
- **FASTA/FASTQ parsing**: 46% of load time (scalar byte-scanning)
- **Hash map building**: 48% of load time (string key copies + rehashing)
- **GZIP decompression**: 5-9% of load time (already optimized)

**Problem:** For a 100MB GZIP file (599MB uncompressed), total load time was 1.843s — dominated by parsing (0.85s) and indexing (0.88s), not I/O.

**Decision drivers:**
- Parsing and indexing are CPU-bound, not I/O-bound
- Existing `simd_find_char()` infrastructure was unused in parser hot loops
- `robin_hood::unordered_flat_map` development has stalled
- Thread-local maps lacked pre-reservation, causing mid-parse rehashing

---

## Decision Outcome

**Chosen (v1.2.0):** Three-pronged optimization: SIMD-accelerated boundary scanning, `ankerl::unordered_dense` hash map, and pre-reserved thread-local maps.

### 1. SIMD-Accelerated Record Boundary Scanning

**Problem:** Parsers used scalar byte-scanning loops (`while (ptr < end && *ptr != '\n') ++ptr`) — 1 byte/iteration on hot paths.

**Solution:** Replace scalar loops with `simd_find_char()` (32 bytes/iter AVX2, 16 bytes/iter NEON):

```cpp
// Before (scalar): O(total_bytes)
while (ptr < end && *ptr != '\n' && *ptr != '\r') ++ptr;

// After (SIMD): O(records) — jumps to next delimiter
ptr = simd_find_char(ptr, end, '\n');
```

**Scope:**
- `parseFastaInternal()`: newline scanning, `>` record boundary detection
- `parseFastqInternal()`: `@`, `\n`, `+` scanning
- `parseFastaMultithreadedTemplate()`: sequence boundary and newline loops
- `parseFastqMultithreadedTemplate()`: newline scanning with `\r`-trimming

**Trade-offs:**
- ID scans remain scalar (10-30 bytes, SIMD overhead exceeds benefit)
- CRLF handling requires explicit `\r`-trimming after `simd_find_char` calls
- `simd_find_char` is read-only, safe for concurrent multithreaded use

### 2. ankerl::unordered_dense Hash Map

**Problem:** `robin_hood::unordered_flat_map` — good but development stalled, `std::string` key copies per record (500K heap allocations for 100MB WGS).

**Solution:** Replace with `ankerl::unordered_dense::map` (header-only, MIT license):

```cpp
// MapDefs.h
template<typename Key, typename Value>
using HashMap = ankerl::unordered_dense::map<Key, Value>;

template<typename Value>
using StringHashMap = ankerl::unordered_dense::map<
    std::string, Value, TransparentStringHash, std::equal_to<>>;
```

**Integration:** CMake FetchContent with `v4.4.0` pin.

**Trade-offs:**
- Transparent lookup preserved (`is_transparent` tag + `std::equal_to<>`)
- Binary serialization format unchanged (stores raw keys, not internal hashes)
- Development actively maintained (Martinus Ankerl, Swiss-table inspired)

### 3. Pre-Reserved Thread-Local Maps

**Problem:** Multithreaded parsers created `std::vector<MapType> thread_caches(num_threads)` without reserving capacity, causing mid-parse rehashing.

**Solution:** Pre-reserve with conservative heuristic before worker spawn:

```cpp
// FASTA: assume ~500 bytes/record
const size_t est_per_thread = chunk_size / 500;
for (auto& cache : thread_caches) {
    cache.reserve(static_cast<size_t>(est_per_thread * 1.25));
}
```

**Heuristic rationale:**
- FASTA: `chunk_size / 500` (header ~50B + sequence ~450B)
- FASTQ: `chunk_size / 600` (header ~50B + seq ~250B + qual ~250B + newlines)
- 1.25x safety margin absorbs boundary alignment variance

**Trade-offs:**
- Initial heuristic `chunk_size / 100` caused 10x over-reservation with large records (seq_len=1000)
- Fixed to `/500` and `/600` — conservative but avoids cache thrashing
- Future: content sampling for dynamic estimation

---

## Alternatives Considered

### Hash Maps
| Library | Chosen | Rationale |
|---------|--------|-----------|
| `robin_hood::unordered_flat_map` | ❌ Replaced | Development stalled, no active maintenance |
| `absl::flat_hash_map` | ❌ Not chosen | Google-specific, complex build integration |
| `tsl::robin_map` | ❌ Not chosen | No transparent hash support |
| `parallel-hashmap` | ❌ Not chosen | Header-only but bloated, slower inserts |
| **ankerl::unordered_dense** | ✅ Chosen | Header-only, MIT, actively maintained, fastest inserts, transparent lookup |

### SIMD Strategies
| Strategy | Chosen | Rationale |
|----------|--------|-----------|
| Compiler auto-vectorization | ❌ Not chosen | Unreliable across compilers, poor NEON support |
| Manual SIMD intrinsic loops | ✅ Chosen | `simd_find_char()` already existed, predictable performance |
| External library (ISPC) | ❌ Not chosen | Extra build dependency, no C++ interop |

### Pre-Reservation Heuristics
| Heuristic | Chosen | Rationale |
|-----------|--------|-----------|
| `chunk_size / 100` | ❌ Abandoned | 10x over-reservation with large records (22% MT regression) |
| `chunk_size / 500` (FASTA) | ✅ Chosen | Conservative, avoids cache thrashing |
| `chunk_size / 600` (FASTQ) | ✅ Chosen | Accounts for 4 lines per record |

---

## Consequences

### Positive
- **86% load time reduction** (1.843s → 0.245s for 100MB GZIP)
- **30% memory reduction** (263MB → 185MB peak RSS)
- Zero API changes — transparent to users
- Cross-platform (AVX2/NEON/scalar dispatch)

### Negative
- `simd_find_char` only scans for `\n` — CRLF requires explicit `\r`-trimming
- `ankerl::unordered_dense` adds external dependency (FetchContent, header-only)
- Pre-reservation heuristic is conservative — may under-reserve for very small records

### Neutral
- Binary cache format unchanged — no migration needed
- `normalizeFastaArena()` bug discovered and fixed (resize zeroing trailing `\n`)
- WGS lookup throughput (12M OPS/s) below 15M target but exceeds 8M regression threshold

---

## Performance Validation

| Metric | v1.1.0 | v1.2.0 | Δ |
|--------|--------|--------|---|
| 100MB GZIP load time | 1.843s | 0.245s | **-86%** |
| Peak RSS | 263MB | 185MB | **-30%** |
| Lookup throughput (WGS) | — | 12M OPS/s | >8M threshold |
| Lookup throughput (long-read) | — | 28-81M OPS/s | >15M threshold |

---

## References

- `simd_find_char()`: `include/SimdUtils.h` — AVX2/NEON/scalar runtime dispatch
- `ankerl::unordered_dense`: https://github.com/martinus/unordered_dense (v4.4.0)
- `TransparentStringHash`: `include/MapDefs.h` — heterogeneous lookup support
- Benchmark data: `benchmarks/README.md`

---

## Contributors

**ADR Authors:**
- Adnan Raza (Woosflex)

---

**Status:** ✅ Accepted (v1.2.0)  
**Last Updated:** 2026-06-09  
**Version:** 1.2.0 "Caladbolg"

*"Trace On" - Decompressing genomic data at the speed of thought.*
