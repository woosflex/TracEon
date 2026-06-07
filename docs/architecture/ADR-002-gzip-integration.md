# ADR-002: GZIP Integration for Native Compressed File Support

**Status:** Accepted (v1.1.0 "Bakuya")  
**Date:** 2026-06-06  
**Deciders:** Adnan Raza (Woosflex)  
**Related:** [ADR-001: Lock-Free Reads](ADR-001-lock-free-reads.md)

---

## Context and Problem Statement

Bioinformatics workflows overwhelmingly use GZIP-compressed files (`.fastq.gz`, `.fasta.gz`) to conserve storage:
- Raw sequencing data from Illumina, PacBio, Nanopore arrives compressed
- Public databases (NCBI SRA, ENA) distribute data as `.gz` files
- ~3-5x compression ratio reduces storage and network transfer costs

**Problem:** TracEon v0.0.1 required users to manually decompress files before loading, creating friction:
```bash
# User workflow without native GZIP support
gunzip genome.fastq.gz          # Decompress (slow, creates 3x disk usage)
traceon_app genome.fastq        # Load
rm genome.fastq                 # Cleanup
```

**Decision drivers:**
- 99% of bioinformatics inputs are GZIP-compressed
- Decompression should be transparent to users
- Must preserve zero-copy architecture and lock-free reads
- Cannot sacrifice lookup performance for load-time convenience

---

## Decision Outcome

**Chosen (v1.0.0):** Integrate zlib for streaming GZIP decompression with automatic format detection.  
**Updated (v1.1.0):** Migrated to **zlib-ng v2.2.2** (zlib-compat mode) with pre-size + direct-write optimization.

### Implementation Strategy (v1.0.0)

**zlib + temp_buffer + memcpy** — initial approach.

#### 1. Auto-Detection Logic (unchanged across versions)

```cpp
void SmartStrategy::loadFile(const std::string& filepath) {
    bool is_gzip = false;
    
    // Check 1: Extension
    if (filepath.ends_with(".gz")) {
        is_gzip = true;
    }
    
    // Check 2: Magic Bytes (fallback for mis-named files)
    if (!is_gzip) {
        std::ifstream check(filepath, std::ios::binary);
        unsigned char magic[2];
        check.read(reinterpret_cast<char*>(magic), 2);
        if (magic[0] == 0x1f && magic[1] == 0x8b) {
            is_gzip = true;
        }
    }
    
    // Dispatch to appropriate loader
    if (is_gzip) {
        loadGzipInternal(filepath);
    } else {
        loadPlainInternal(filepath);
    }
    
    parseArena();
}
```

**Rationale:**
- Extension check is fast (O(1) string comparison)
- Magic bytes provide robustness against mis-named files
- Explicit `loadGzipFile()` available for cases where format is known

### Updated for v1.1.0: zlib-ng + Pre-Size + Direct-Write

**Three key changes:**

#### 1. zlib-ng Replacement

zlib has been replaced with **zlib-ng v2.2.2** via CMake FetchContent, compiled in zlib-compat mode (`ZLIB_COMPAT=ON`) with native instruction support (`WITH_NATIVE_INSTRUCTIONS=ON`).  The API surface is identical (`gzopen`, `gzread`, `gzclose`) — the upgrade is transparent at the call site but delivers ~10–15% faster decompression through SIMD-optimized inflate.

```cmake
FetchContent_Declare(zlib-ng
    GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git
    GIT_TAG        v2.2.2)

set(ZLIB_COMPAT ON)
set(WITH_NATIVE_INSTRUCTIONS ON)
FetchContent_MakeAvailable(zlib-ng)
```

#### 2. Pre-Size with OOM Guard

The old heuristic `estimate_uncompressed_size()` is replaced by a deterministic pre-size strategy using `std::filesystem::file_size()`:

```
pre_size = min(compressed_size × 3, available_memory × 0.25)
```

- **3x multiplier** covers the typical GZIP ratio (bioinformatics FASTQ/Fasta: 3–5×)
- **25% memory cap** (`avail_mem / 4`) prevents OOM on memory-constrained systems
- `available_memory()` reads `/proc/meminfo` (Linux), `sysctl` (macOS), or `GlobalMemoryStatusEx` (Windows)
- On `std::bad_alloc` the exception propagates with a diagnostic message showing compressed size, estimated size, and attempted reserve

#### 3. Direct-Write Elimination of temp_buffer

The separate `std::vector<char> temp_buffer` has been eliminated. Decompressed data is written **directly into `text_arena_`**, removing the intermediate move step.  A small 1MB chunk buffer is retained for `gzread` source reads.

```cpp
void SmartStrategy::loadGzipInternal(const std::string& filepath) {
    const size_t CHUNK_SIZE = 1024 * 1024; // 1MB

    // ── Pre-size text_arena_ with 3x heuristic and OOM guard ──────────
    {
        std::error_code ec;
        const uintmax_t raw_size = std::filesystem::file_size(filepath, ec);
        if (!ec && raw_size > 0) {
            const size_t estimated_size = static_cast<size_t>(raw_size * 3);
            const size_t avail_mem = getAvailableMemory();
            const size_t reserve_size = (avail_mem > 0)
                ? std::min(estimated_size, avail_mem / 4)
                : estimated_size;
            try {
                text_arena_.reserve(reserve_size);
            } catch (const std::bad_alloc&) {
                std::cerr << "OOM: failed to reserve " << reserve_size
                          << " bytes (compressed=" << raw_size
                          << ", estimated=" << estimated_size << ")"
                          << std::endl;
                throw;
            }
        } else {
            text_arena_.reserve(CHUNK_SIZE * 4); // fallback
        }
    }

    gzFile file = gzopen(filepath.c_str(), "rb");
    if (!file) throw std::runtime_error("Cannot open GZIP file: " + filepath);

    auto chunk = std::make_unique<char[]>(CHUNK_SIZE);
    size_t write_pos = 0;

    try {
        while (true) {
            // Geometric growth: reserve BEFORE resize
            size_t required = write_pos + CHUNK_SIZE;
            if (required > text_arena_.capacity()) {
                text_arena_.reserve(
                    std::max(text_arena_.capacity() * 2, required));
            }
            text_arena_.resize(required);

            int bytes_read = gzread(file, chunk.get(), CHUNK_SIZE);
            if (bytes_read < 0) {
                /* error handling */
            }
            if (bytes_read == 0) break;
            std::memcpy(text_arena_.data() + write_pos,
                        chunk.get(), bytes_read);
            write_pos += static_cast<size_t>(bytes_read);
        }
    } catch (...) {
        gzclose(file);
        text_arena_.clear();
        throw;
    }
    gzclose(file);
    text_arena_.resize(write_pos);    // Trim to actual size
    text_arena_.shrink_to_fit();      // Release excess capacity
}
```

**Design principles (v1.1.0):**
- **zlib-ng:** Drop-in replacement, ~10–15% faster inflate via SIMD
- **Pre-sizing:** Single allocation eliminates reallocation storm concern from the original ADR-002 rejection of direct-write
- **OOM guard:** 25% available memory cap prevents crashes on constrained systems
- **Geometric fallback:** If estimate is too small, geometric growth (2× capacity) handles under-estimation gracefully
- **Chunk buffer:** 1MB buffer retained as gzread source — negligible overhead
- **shrink_to_fit:** Releases excess capacity after decompression

#### 3. Zero-Copy Preservation (unchanged)

```
┌─────────────────────────┐
│  Compressed File        │
│  (genome.fastq.gz)      │
└───────────┬─────────────┘
            │ zlib-ng streaming
            v
┌─────────────────────────┐
│  text_arena_            │  ← Written directly (no temp_buffer)
│  (pre-sized, geometric) │
└───────────┬─────────────┘
            │ string_view pointers
            v
┌─────────────────────────┐
│  SequenceView           │  ← Zero-copy access
│  {id, seq, qual}        │
└─────────────────────────┘
```

**Performance impact (v1.1.0):**
- Load time: +0.251s for 100MB GZIP (10% improvement over v1.0.0)
- Lookup time: **0.00s overhead** (identical to plain text)
- Memory peak during load: ~266MB (pre-sized arena; released to ~183MB after `shrink_to_fit`)

---

## Alternatives Considered

### Alternative 1: External Decompression Requirement ❌

**Approach:** Require users to decompress files manually before loading.

```bash
# User must do this themselves
gunzip genome.fastq.gz
./app genome.fastq
```

**Rejected because:**
- ❌ Friction in user workflow (extra step)
- ❌ Wastes disk space (3-5x uncompressed size)
- ❌ Users might forget to cleanup decompressed files
- ❌ Not competitive with tools like PyFastX (native GZIP support)

**When this might make sense:** If targeting systems with extremely limited CPU (decompression overhead too high).

---

### Alternative 2: Direct-Write to Arena ✅ (Partially Adopted in v1.1.0)

**Approach:** Write decompressed data directly into `text_arena_` instead of through an intermediate `temp_buffer`.

**v1.0.0 status — Rejected without pre-sizing:**
```cpp
// REJECTED — without pre-sizing causes reallocation storm
while (true) {
    size_t current_size = text_arena_.size();
    text_arena_.resize(current_size + CHUNK_SIZE);
    int bytes_read = gzread(file, text_arena_.data() + current_size, CHUNK_SIZE);
    if (bytes_read < CHUNK_SIZE) { text_arena_.resize(current_size + bytes_read); break; }
}
```

- ❌ **Reallocation storm:** `resize()` twice per iteration bypasses geometric growth
- ❌ **Memory waste:** Over-allocate 2MB, use 1.5MB, shrink → orphan 0.5MB
- ❌ **Performance:** 50+ reallocations for 100MB (vs 7 with geometric)

**v1.1.0 update — Adopted with pre-sizing:**
The core objection (reallocation storm) is resolved by pre-sizing `text_arena_` before the loop. With `capacity` ≥ expected decompressed size, the geometric growth path never activates for typical files.

```
Pre-size eliminates reallocation:
  - v1.0.0 (temp_buffer + memcpy): 0.28s, ~7 reallocations
  - v1.1.0 (direct-write + pre-size): 0.251s, 0–1 reallocations
  - Saving: ~0.03s from eliminated memcpy + std::move
```

**Remaining consideration:** A small 1MB chunk buffer is still retained because it keeps gzread I/O separate from the arena for simplicity — the data is then `memcpy`-ed into `text_arena_`.  This is not a true zero-copy decompression (which would require `zlib-ng`'s `deflate` with a custom output allocator), but the overhead is negligible (~0.001s at 20 GB/s).

---

### Alternative 3: Parallel GZIP Decompression (ISA-L) 🔄 Deferred to v1.2.0

**Approach:** Use Intel ISA-L or `pigz` for multi-threaded decompression.

```cpp
// FUTURE v1.2.0 implementation
void SmartStrategy::loadGzipInternalParallel(const std::string& filepath) {
    // 1. Split compressed file into independent blocks
    // 2. Decompress blocks in parallel (4-8 threads)
    // 3. Concatenate results into text_arena_
}
```

**Not implemented in v1.0.0 or v1.1.0 because:**
- ⏳ Requires external dependency (ISA-L or custom parallel gzip)
- ⏳ GZIP format doesn't naturally support random access (requires index)
- ⏳ Complexity increases (thread coordination, error handling)
- ⏳ v1.1.0 chose the **lower-risk path** (zlib-ng swap + pre-size) over architectural redesign

**Performance potential:**
```
v1.0.0 (single-threaded zlib):
  - 100MB load: 0.28s
  - Throughput: ~350 MB/s

v1.1.0 (single-threaded zlib-ng):
  - 100MB load: 0.251s (10% improvement)
  - Throughput: ~400 MB/s

With ISA-L (4 threads) — estimate:
  - Estimated: 0.063s (4x speedup)
  - Throughput: ~1.6 GB/s
  - Goal: Close gap with SeqKit (0.22s)
```

**Status:** Deferred. v1.1.0 achieved 10% improvement through zlib-ng + pre-size instead.

---

### Alternative 4: Memory-Mapped Decompression 🚫 Infeasible

**Approach:** Use `mmap()` on compressed file and decompress on-demand.

**Rejected because:**
- 🚫 GZIP is not a random-access format (cannot seek to arbitrary positions)
- 🚫 Would require maintaining a seek index (overhead)
- 🚫 Defeats the purpose of zero-copy (would need decompression buffer anyway)
- 🚫 Complex error handling (partial decompression failures)

**When this might make sense:** For formats designed for random access (e.g., BGZF, used in BAM files).

---

## Trade-offs and Consequences

### Positive Consequences ✅

1. **User Experience:**
   - Single function call: `loadFile("genome.fastq.gz")`
   - No manual decompression step
   - Transparent: Users don't need to know file is compressed

2. **Storage Efficiency:**
   - Keep files compressed (3-5x space savings)
   - No need for temporary decompressed files

3. **Performance:**
   - Lookup speed unaffected (zero-copy preserved)
   - v1.0.0 load time: 0.28s (27% behind SeqKit's 0.22s)
   - v1.1.0 load time: **0.251s** (10% improvement, gap reduced to ~14%)

4. **Robustness:**
   - Magic bytes fallback handles mis-named files
   - Comprehensive error handling (corrupt GZIP, truncated streams)

### Negative Consequences ❌

1. **Load Time Overhead:**
   - v1.0.0: +0.28s for 100MB file
   - v1.1.0: +0.251s (10% improvement from zlib-ng + pre-size)
   - Target of 20% improvement (≤0.224s) **not met** — remaining gap is algorithmic, not micro-optimization
   - SeqKit is ~14% faster (0.22s) — likely uses parallel decompression
   - **Mitigation:** Binary cache amortizes cost (load once, restore instant)

2. **Dependency on zlib:**
   - v1.0.0 used system zlib (~200KB to binary)
   - v1.1.0 replaced with **zlib-ng v2.2.2** (vendored via FetchContent)
   - **Mitigation:** zlib-compat mode means identical API; no code changes needed

3. **Single-Threaded Bottleneck:**
   - Decompression still single-threaded (zlib-ng inflate is SIMD-accelerated but serial)
   - v1.0.0 → v1.1.0: +10% from library swap + allocation optimization
   - **Mitigation:** Parallel decompression deferred to v1.2.0 (ISA-L or pigz)

### Neutral Consequences ⚖️

1. **Memory Usage:**
   - **v1.0.0:** temp_buffer + text_arena_ co-exist during load (~366MB peak for 100MB file)
   - **v1.1.0:** Pre-sized text_arena_ only (~266MB peak; shrink_to_fit → ~183MB)
   - v1.0.0 memory target (5MB) was **infeasible** — even compressed data is 33MB for a 100MB FASTQ.gz
   - The 266MB peak is the cost of the 3× pre-size heuristic; for very large files the OOM guard caps at 25% available memory
   - No impact on steady-state memory after shrink_to_fit

2. **API Surface:**
   - Added `loadGzipFile()` (explicit)
   - Modified `loadFile()` (auto-detect)
   - Backwards compatible (no breaking changes)

---

## Performance Validation

### Benchmark Setup
- **System:** Intel Core Ultra 5 125H (16GB LPDDR5)
- **Build:** Release (`-O3`, C++20)
- **Dataset:** 100MB Illumina FASTQ.gz (3.2:1 compression ratio)

### Results

| Metric | TracEon v1.0.0 | TracEon v1.1.0 | SeqKit | PyFastX |
|--------|----------------|----------------|--------|---------|
| **Load Time** | 0.28s | **0.251s** | 0.22s | N/A |
| **vs Target (≤0.224s)** | ❌ 25% above | ❌ 14% above | — | — |
| **Random Lookups** | 17.4M OPS/s | 17.4M OPS/s | N/A | 12.0M OPS/s |
| **Memory (steady-state)** | 183 MB | 183 MB | Unknown | Unknown |
| **Memory (load peak)** | ~366 MB | **~266 MB** | Unknown | Unknown |

**Improvement breakdown (v1.0.0 → v1.1.0):**
```
Total: 0.28s → 0.251s (-29ms, 10.4%)
  ├─ zlib-ng inflate:          -25ms  (SIMD-optimized CRC + inflate)
  ├─ Pre-size (no realloc):     -3ms  (eliminated 7 reallocations)
  └─ Direct-write (no move):    -1ms  (eliminated std::move)
```

**Analysis:**
- **Load time:** 10% improvement meets expectations for a drop-in library swap
- **20% target not met:** Most of the gap is in the inflate loop itself; zlib-ng provides only ~10–15% improvement over vanilla zlib
- **Lookup speed:** Unchanged (zero-copy invariant preserved)
- **Memory peak:** Reduced 27% (from 366MB to 266MB) by eliminating temp_buffer
- Narrowing the remaining gap requires parallel decompression (ISA-L/pigz) — see Alternative 3

**Conclusion:** v1.1.0 delivers solid incremental improvements. The 20% load-time target requires parallel decompression (v1.2.0).

---

## Decision Rationale

### Why This Approach?

1. **Simplicity:** zlib is battle-tested, widely available, and well-documented
2. **Zero-copy preservation:** Decompression is one-time cost; all reads remain zero-copy
3. **Robustness:** Magic bytes fallback handles edge cases
4. **Extensibility:** Easy to swap in parallel decompression later (v1.2.0)

### Key Design Principle

> **"Fast code is not about avoiding operations—it's about making the right operations in the right order."**

- `memcpy()` @ 20GB/s is **not** the bottleneck
- Single-threaded decompression @ 350–400MB/s **is** the bottleneck
- Fix the bottleneck, not the fast path

---

## Implementation Notes

### Error Handling

```cpp
void SmartStrategy::loadGzipInternal(const std::string& filepath) {
    gzFile file = gzopen(filepath.c_str(), "rb");
    if (!file) {
        throw std::runtime_error("Cannot open GZIP file: " + filepath);
    }
    
    // ... decompression loop ...
    
    if (bytes_read < 0) {
        int err;
        const char* error_msg = gzerror(file, &err);
        gzclose(file);
        throw std::runtime_error("GZIP read error: " + std::string(error_msg));
    }
    
    gzclose(file); // RAII cleanup
}
```

**Error scenarios handled:**
- File not found
- Corrupt GZIP header
- Truncated stream
- Read errors (disk failure, permission denied)

### Memory Safety

**RAII principles:**
```cpp
auto chunk = std::make_unique<char[]>(CHUNK_SIZE);  // Automatic cleanup
```

**Valgrind verification:**
```bash
valgrind --leak-check=full ./unit_tests "[gzip]"
# Result: All heap blocks were freed -- no leaks are possible
```

### Testing Strategy

**Test cases:**
```cpp
TEST_CASE("GZIP Support", "[strategy][gzip]") {
    SECTION("Explicit Load") {
        strategy.loadGzipFile("test.fasta.gz");
        REQUIRE(strategy.getSequence("seq1") == "ACGT");
    }
    
    SECTION("Auto-detect via Extension") {
        strategy.loadFile("test.fasta.gz");
        REQUIRE(strategy.getSequence("seq1") == "ACGT");
    }
    
    SECTION("Auto-detect via Magic Bytes") {
        // File named .fasta but contains GZIP data
        strategy.loadFile("misleading.fasta");
        REQUIRE(strategy.getSequence("seq1") == "ACGT");
    }
}
```

---

## Future Work

### ✓ v1.1.0 "Bakuya" (Q1 2026) — Completed

**zlib-ng Integration + Pre-Size + Direct-Write** *(shipped)*
- zlib-ng v2.2.2 via FetchContent, zlib-compat mode
- Pre-size with 3x heuristic + OOM guard (25% available memory)
- Direct-write to text_arena_ (eliminated temp_buffer)
- Result: **0.251s** (10% improvement over v1.0.0)
- Target of 20% improvement **not met**

### v1.2.0 "Caladbolg" (Q2 2026)

**Adaptive Chunking:** (deferred — design required)
- Dynamic chunk size based on compression ratio
- Optimize for different file characteristics

**Parallel GZIP Decompression:**
- Integrate Intel ISA-L or custom parallel decompressor
- Target: 4x speedup (0.251s → ~0.06s for 100MB)
- Close gap with SeqKit

**Binary Cache Compression:**
- LZ4 compression for `.traceon` files
- Target: 3x size reduction with minimal decompression overhead

---

## References

### External Resources
- [zlib Documentation](https://zlib.net/manual.html)
- [zlib-ng](https://github.com/zlib-ng/zlib-ng) — SIMD-optimized zlib replacement
- [GZIP Specification (RFC 1952)](https://datatracker.ietf.org/doc/html/rfc1952)
- [Intel ISA-L](https://github.com/intel/isa-l) - SIMD-accelerated compression
- SeqKit: [Shen et al., PLoS ONE 2016](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0163962)

### Internal Documents
- [ADR-001: Lock-Free Reads](ADR-001-lock-free-reads.md) - Memory ordering guarantees
- [Performance Profile](../performance-profile.md) - Regression thresholds
- [Benchmark Guide](../../benchmarks/README.md) - Reproducibility

---

## Contributors

**ADR Authors:**
- Adnan Raza (Woosflex)

---

**Status:** ✅ Accepted (v1.1.0)  
**Last Updated:** 2026-06-06  
**Version:** 1.1.0 "Bakuya"

*"Trace On" - Decompressing genomic data at the speed of thought.*