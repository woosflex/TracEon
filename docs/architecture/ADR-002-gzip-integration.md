# ADR-002: GZIP Integration for Native Compressed File Support

**Status:** Accepted  
**Date:** 2025-12-18  
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

**Chosen:** Integrate zlib for streaming GZIP decompression with automatic format detection.

### Implementation Strategy

#### 1. Auto-Detection Logic
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
    
    parseArena(); // Parse after decompression
}
```

**Rationale:**
- Extension check is fast (O(1) string comparison)
- Magic bytes provide robustness against mis-named files
- Explicit `loadGzipFile()` available for cases where format is known

#### 2. Streaming Decompression
```cpp
void SmartStrategy::loadGzipInternal(const std::string& filepath) {
    gzFile file = gzopen(filepath.c_str(), "rb");
    
    const size_t CHUNK_SIZE = 1024 * 1024; // 1MB
    std::vector<char> temp_buffer;
    
    // Heuristic: Reserve 3x compressed size
    temp_buffer.reserve(estimate_uncompressed_size(filepath));
    
    // RAII buffer
    auto chunk = std::make_unique<char[]>(CHUNK_SIZE);
    
    while (true) {
        int bytes_read = gzread(file, chunk.get(), CHUNK_SIZE);
        if (bytes_read <= 0) break;
        
        // Geometric growth (O(log n) reallocations)
        size_t old_size = temp_buffer.size();
        temp_buffer.resize(old_size + bytes_read);
        std::memcpy(temp_buffer.data() + old_size, chunk.get(), bytes_read);
    }
    
    gzclose(file);
    
    // Zero-copy move into arena
    text_arena_ = std::move(temp_buffer);
}
```

**Design principles:**
- **Streaming:** Process file in chunks (doesn't load entire compressed file into RAM)
- **Geometric growth:** `std::vector` doubles capacity automatically (minimal reallocations)
- **Single copy:** `memcpy()` from decompression buffer to temp_buffer, then `std::move()` to arena
- **RAII:** `std::unique_ptr` ensures cleanup even on exceptions

#### 3. Zero-Copy Preservation

**Critical invariant:** Once decompression completes, all subsequent access is zero-copy.

```
┌─────────────────────────┐
│  Compressed File        │
│  (genome.fastq.gz)      │
└───────────┬─────────────┘
            │ zlib streaming
            v
┌─────────────────────────┐
│  temp_buffer            │  ← Temporary (during load)
│  (std::vector<char>)    │
└───────────┬─────────────┘
            │ std::move (zero-copy)
            v
┌─────────────────────────┐
│  text_arena_            │  ← Permanent (after load)
│  (zero allocations)     │
└───────────┬─────────────┘
            │ string_view pointers
            v
┌─────────────────────────┐
│  SequenceView           │  ← Zero-copy access
│  {id, seq, qual}        │
└─────────────────────────┘
```

**Performance impact:**
- Load time: +0.28s for 100MB GZIP (one-time cost)
- Lookup time: **0.00s overhead** (identical to plain text)
- Memory: +0 bytes (no additional structures needed)

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

### Alternative 2: Direct-Write to Arena ❌

**Approach:** Read directly into `text_arena_` by repeatedly calling `resize()`.

```cpp
// REJECTED IMPLEMENTATION
while (true) {
    size_t current_size = text_arena_.size();
    text_arena_.resize(current_size + CHUNK_SIZE);  // Pre-allocate
    
    int bytes_read = gzread(file, text_arena_.data() + current_size, CHUNK_SIZE);
    
    if (bytes_read < CHUNK_SIZE) {
        text_arena_.resize(current_size + bytes_read);  // Shrink
        break;
    }
}
```

**Rejected because:**
- ❌ **Reallocation storm:** Calling `resize()` twice per iteration bypasses geometric growth
- ❌ **Memory waste:** Over-allocate 2MB, use 1.5MB, shrink → orphan 0.5MB (cache pollution)
- ❌ **Performance:** 50+ reallocations for 100MB file (vs 7 with geometric growth)

**Benchmark evidence:**
```
Previous (temp_buffer + memcpy):
  - 100MB load: 0.28s
  - Reallocations: ~7 (log2(100MB/1MB))
  - memcpy overhead: ~0.1s (@ 1GB/s)

Direct-write (rejected):
  - 100MB load: ~1.5s (5.4x slower!)
  - Reallocations: ~50 (linear growth)
  - Overhead from realloc: ~1.25s
```

**Rule of thumb:** 1 reallocation ≈ 50-100 memcpy operations (same size).

**When this might make sense:** If `std::vector` didn't have geometric growth (e.g., custom allocator).

---

### Alternative 3: Parallel GZIP Decompression (ISA-L) 🔄 Future

**Approach:** Use Intel ISA-L or `pigz` for multi-threaded decompression.

```cpp
// FUTURE v1.1.0 implementation
void SmartStrategy::loadGzipInternalParallel(const std::string& filepath) {
    // 1. Split compressed file into independent blocks
    // 2. Decompress blocks in parallel (4-8 threads)
    // 3. Concatenate results into text_arena_
}
```

**Not implemented in v1.0.0 because:**
- ⏳ Requires external dependency (ISA-L or custom parallel gzip)
- ⏳ GZIP format doesn't naturally support random access (requires index)
- ⏳ Complexity increases (thread coordination, error handling)

**Performance potential:**
```
Current (single-threaded zlib):
  - 100MB load: 0.28s
  - Throughput: ~350 MB/s

With ISA-L (4 threads):
  - Estimated: 0.07s (4x speedup)
  - Throughput: ~1.4 GB/s
  - Goal: Match/beat SeqKit (0.22s)
```

**Status:** Planned for v1.1.0 "Bakuya" (Q1 2026)

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
   - Load time competitive (0.28s vs SeqKit's 0.22s)

4. **Robustness:**
   - Magic bytes fallback handles mis-named files
   - Comprehensive error handling (corrupt GZIP, truncated streams)

### Negative Consequences ❌

1. **Load Time Overhead:**
   - GZIP decompression adds ~0.28s for 100MB file
   - SeqKit is ~27% faster (0.22s) - likely uses parallel decompression
   - **Mitigation:** Binary cache amortizes cost (load once, restore instant)

2. **Dependency on zlib:**
   - Adds ~200KB to binary size
   - Requires zlib headers at compile time
   - **Mitigation:** zlib is ubiquitous (included in most systems)

3. **Single-Threaded Bottleneck:**
   - Decompression doesn't scale beyond 1 core
   - **Mitigation:** Planned for v1.1.0 (parallel GZIP)

### Neutral Consequences ⚖️

1. **Memory Usage:**
   - Temporary buffer during load (~3x compressed size)
   - Freed immediately after decompression
   - No impact on steady-state memory

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

| Metric | TracEon v1.0.0 | SeqKit | PyFastX |
|--------|----------------|--------|---------|
| **Load Time** | 0.28s | 0.22s | N/A |
| **Random Lookups** | 17.4M OPS/s | N/A | 12.0M OPS/s |
| **Memory** | 183 MB | Unknown | Unknown |

**Analysis:**
- Load time competitive but not best-in-class (SeqKit 27% faster)
- Lookup speed excellent (45% faster than PyFastX)
- Memory efficiency maintained (5x better than BioPython)

**Conclusion:** Performance trade-offs acceptable for v1.0.0. Will optimize in v1.1.0.

---

## Decision Rationale

### Why This Approach?

1. **Simplicity:** zlib is battle-tested, widely available, and well-documented
2. **Zero-copy preservation:** Decompression is one-time cost; all reads remain zero-copy
3. **Robustness:** Magic bytes fallback handles edge cases
4. **Extensibility:** Easy to swap in parallel decompression later (v1.1.0)

### Key Design Principle

> **"Fast code is not about avoiding operations—it's about making the right operations in the right order."**

- `memcpy()` @ 20GB/s is **not** the bottleneck
- Single-threaded decompression @ 350MB/s **is** the bottleneck
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

### v1.1.0 "Bakuya" (Q1 2026)

**Parallel GZIP Decompression:**
- Integrate Intel ISA-L or custom parallel decompressor
- Target: 4x speedup (0.28s → 0.07s for 100MB)
- Close gap with SeqKit

**Adaptive Chunking:**
- Dynamic chunk size based on compression ratio
- Optimize for different file characteristics

### v1.2.0 "Caladbolg" (Q2 2026)

**Binary Cache Compression:**
- LZ4 compression for `.traceon` files
- Target: 3x size reduction with minimal decompression overhead

---

## References

### External Resources
- [zlib Documentation](https://zlib.net/manual.html)
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

**Status:** ✅ Accepted  
**Last Updated:** 2025-12-04  
**Version:** 1.0.0 "Avalon"

*"Trace On" - Decompressing genomic data at the speed of thought.*