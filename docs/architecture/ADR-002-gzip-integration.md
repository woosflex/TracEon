# ADR-002: GZIP Integration for Native Compressed File Support

**Status:** Accepted (v1.2.0 "Caladbolg")  
**Date:** 2026-06-09  
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
**Updated (v1.2.0):** SIMD-accelerated boundary scanning + ankerl::unordered_dense hash map + pre-reserved thread-local maps — **86% load time reduction**.

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

### Updated for v1.2.0: SIMD Boundary Scanning + ankerl::unordered_dense

**Three key changes:**

#### 1. SIMD-Accelerated Record Boundary Scanning

The parsers now use `simd_find_char()` which scans 32 bytes per iteration (AVX2) or 16 bytes per iteration (NEON) instead of the previous byte-by-byte loop. Runtime dispatch via `__builtin_cpu_supports("avx2")` on GCC/Clang x86-64, with a `std::memchr` scalar fallback for older CPUs.

```cpp
// normalizeFastaArena() — SIMD everywhere
while (read < rec_end) {
    const char* line_nl = simd_find_char(read, rec_end, '\n');  // 32 bytes/iter
    // bulk copy non-newline run...
}
```

**Integrated in:**
- `normalizeFastaArena()` — newline scanning in post-process compaction
- `parseFastaInternal()` / `parseFastqInternal()` — `'>'`, `'@'`, `'+'`, `'\n'` boundary scanning
- `parseFastaMultithreadedTemplate()` / `parseFastqMultithreadedTemplate()` — same SIMD acceleration in worker threads

#### 2. ankerl::unordered_dense Hash Map

Replaced `robin_hood::unordered_flat_map` with `ankerl::unordered_dense::map` (Swiss-table design) via CMake FetchContent. The new map offers:

- **Better cache locality** — fewer cache misses during insertion and lookup
- **~0.8s insert time savings** for 100MB datasets (dominant component of parse time)
- **Lower memory overhead** — more compact table representation saves ~30% peak RSS
- **Same API surface** — header-only, transparent to calling code

```cmake
FetchContent_Declare(ankerl_unordered_dense
    GIT_REPOSITORY https://github.com/martinus/unordered_dense.git
    GIT_TAG        v4.4.0
)
FetchContent_MakeAvailable(ankerl_unordered_dense)
```

#### 3. Pre-Reserved Thread-Local Maps

Multithreaded parsers now pre-reserve thread-local hash maps before the worker phase:

```cpp
// Conservative heuristic to avoid over-reservation
const size_t est_per_thread = chunk_size / 500;  // FASTA
// chunk_size / 600 for FASTQ
for (auto& cache : thread_caches) {
    cache.reserve(static_cast<size_t>(est_per_thread * 1.25));
}
```

This eliminates mid-parse rehashing entirely, which was a significant source of overhead in the merge-heavy multithreaded path.

**Performance impact (v1.2.0):**
- Load+parse time: **0.245s** for 100MB GZIP (**86% reduction** from v1.1.0 baseline)
- Lookup time: **0.00s overhead** (unchanged — zero-copy invariant preserved)
- Memory peak during load: **~185MB** (30% reduction from 263MB)
- Improvement breakdown:

```
Total: 1.843s → 0.245s (-1.598s, 87%)
  ├─ simd_find_char() boundary scanning:   -0.6s
  ├─ ankerl::unordered_dense map insert:   -0.8s
  ├─ Pre-reserved thread-local maps:       -0.1s
  └─ normalizeFastaArena() trailing \n fix: -0.05s
```

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
   - v1.0.0 load time: 1.843s
   - v1.1.0 load time: 1.843s (dominant cost was hash map insertion + byte-by-byte parsing)
   - v1.2.0 load time: **0.245s** (86% reduction — SIMD + Swiss-table hash map)
   - Gap with SeqKit reduced from ~14% (v1.1.0) to **~9%** (v1.2.0)

4. **Robustness:**
   - Magic bytes fallback handles mis-named files
   - Comprehensive error handling (corrupt GZIP, truncated streams)

### Negative Consequences ❌

1. **Load Time Overhead:**
   - v1.0.0: 1.843s for 100MB file (full load+parse)
   - v1.1.0: 1.843s (no improvement in parse path)
   - v1.2.0: **0.245s** (86% reduction — SIMD + hash map optimization)
   - The v1.2.0 improvements addressed the true bottleneck: hash map insertion (48%) and FASTA parsing (46%)
   - **Mitigation:** Binary cache amortizes cost (load once, restore instant)

2. **Dependency on zlib:**
   - v1.0.0 used system zlib (~200KB to binary)
   - v1.1.0 replaced with **zlib-ng v2.2.2** (vendored via FetchContent)
   - **Mitigation:** zlib-compat mode means identical API; no code changes needed

3. **Single-Threaded Bottleneck:**
   - Decompression is still single-threaded, but overall load time is dominated by parsing + map build — not inflate
   - v1.2.0 analysis: decompression (zlib-ng inflate) accounts for only 5-9% of total load time
   - **Mitigation:** SIMD parsing + hash map optimization addressed the actual bottleneck; parallel decompression is now lower priority

### Neutral Consequences ⚖️

1. **Memory Usage:**
   - **v1.0.0:** temp_buffer + text_arena_ co-exist during load (~366MB peak for 100MB file)
   - **v1.1.0:** Pre-sized text_arena_ only (~266MB peak; shrink_to_fit → ~183MB)
   - **v1.2.0:** ankerl::unordered_dense compact table layout reduces peak to **~185MB** (30% reduction from 263MB)
   - v1.0.0 memory target (5MB) was **infeasible** — even compressed data is 33MB for a 100MB FASTQ.gz
   - The 185MB peak in v1.2.0 is only 2MB above steady-state — load peak nearly equals final footprint
   - No impact on lookup performance or API surface

2. **API Surface:**
   - Added `loadGzipFile()` (explicit)
   - Modified `loadFile()` (auto-detect)
   - Backwards compatible (no breaking changes)
   - Map typedefs in `MapDefs.h` updated from `robin_hood` to `ankerl` — fully transparent to callers

---

## Performance Validation

### Benchmark Setup
- **System:** Intel Core Ultra 5 125H (16GB LPDDR5)
- **Build:** Release (`-O3`, C++20)
- **Dataset:** 100MB Illumina FASTQ.gz (3.2:1 compression ratio)

### Results

| Metric | TracEon v1.0.0 | TracEon v1.1.0 | TracEon v1.2.0 | SeqKit | PyFastX |
|--------|----------------|----------------|----------------|--------|---------|
| **Load+Parse Time** | 1.843s | 1.843s | **0.245s** | 0.22s | N/A |
| **vs Target (≤0.224s)** | ❌ 8x above | ❌ 8x above | ⚠️ **9% above** | — | — |
| **Random Lookups** | 17.4M OPS/s | 17.4M OPS/s | 17.4M OPS/s | N/A | 12.0M OPS/s |
| **Memory (steady-state)** | 183 MB | 183 MB | **185 MB** | Unknown | Unknown |
| **Memory (load peak)** | ~366 MB | ~266 MB | **~185 MB** | Unknown | Unknown |

**Improvement breakdown (v1.1.0 → v1.2.0):**
```
Total: 1.843s → 0.245s (-1.598s, 86.7%)
  ├─ simd_find_char() boundary scanning:   -0.6s  (AVX2/NEON, 32 bytes/iter)
  ├─ ankerl::unordered_dense map insert:   -0.8s  (Swiss-table vs Robin Hood)
  ├─ Pre-reserved thread-local maps:       -0.1s  (eliminated rehashing)
  └─ normalizeFastaArena() trailing \n fix: -0.05s
```

**Analysis:**
- **Load time:** 86% improvement — far exceeding the original 20% target
- **SIMD wins:** The biggest win from `simd_find_char()` comes in `normalizeFastaArena()` where multi-line FASTA sequence compaction previously scanned byte-by-byte
- **Hash map dominance:** Insertion into the hash map was the single largest component of parse time (48% in v1.1.0). The ankerl Swiss-table implementation cuts this dramatically
- **Lookup speed:** Unchanged (zero-copy invariant preserved)
- **Memory peak:** Reduced 30% (from 263MB to 185MB) via ankerl's compact table layout
- Gap with SeqKit now reduced to ~9% (0.245s vs 0.22s)

**Conclusion:** v1.2.0 delivers a paradigm shift in load performance. The combination of SIMD boundary scanning and Swiss-table hashing addresses the two dominant bottlenecks identified in v1.1.0 profiling (FASTA parsing at 46%, hash map build at 48%).

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

## Binary Cache Compression: LZ4 Integration (v1.3.0)

### Problem Context

Binary cache files (`.traceon`) created by v1.2.0 are uncompressed (v1 format, magic `"TRO\x01"`). A 100MB FASTA file expands to ~105MB binary cache — no compression benefit despite high data redundancy:

```
Binary cache size (100MB FASTA):
  v1 uncompressed: ~105MB
  Typical repetition: 75-80% nucleotide ('A', 'C', 'G', 'T') + frequent newlines
  Expected LZ4 ratio: 3-4x
```

### Decision: LZ4 Format v2

**Chosen:** LZ4 compression for binary cache payloads with format version bump to `0x02`.

**Rationale:**
- **Decompression speed:** 1–2 GB/s (2-4x faster than GZIP's ~400 MB/s)
- **Compression ratio:** 3–4x (bioinformatics data is highly repetitive)
- **Backward compatibility:** `loadBinary()` detects version byte and routes to appropriate decompressor
- **Low CPU overhead:** LZ4 faster to decompress than GZIP read time (~35MB compressed → 105MB uncompressed in ~0.02s)
- **Already vendored:** LZ4 v1.10.0 at `third_party/lz4/` with CMake integration

### Implementation

#### Format v2 Header Layout

```
[0-3]   "TRO\x02"                   ← Version bumped to 0x02
[4]     mode (0=GenomeIndex, 1=NGSIndex)
[5-12]  original_size (uint64_t)    ← Uncompressed payload size
[13-20] compressed_size (uint64_t)  ← Compressed payload size
[21+]   LZ4-compressed payload      ← count + all records (identical to v1 wire format)
```

**Key design:**
- Decompressed payload is identical to v1 format (count + records)
- Parsing logic shared between v1 and v2 (only data source differs: mmap vs text_arena_)
- `safe_advance()` lambda reusable after decompression into `text_arena_`

#### Serialization Path (saveBinary)

```cpp
1. Serialize payload → uncompressed std::vector<char> buf
   (contains count + all records, identical to v1 format)

2. Compress with LZ4
   compressed_size = LZ4_compress_default(buf.data(), compressed.data(), 
                                          buf.size(), LZ4_compressBound(buf.size()))

3. Write header
   out.write("TRO\x02", 4)
   out.write(&mode, 1)
   out.write(&original_size, 8)       // buf.size()
   out.write(&compressed_size, 8)     // result from compress
   out.write(compressed.data(), compressed_size)
```

#### Decompression Path (loadBinary)

```cpp
1. mmap() compressed file

2. Read header
   magic_bytes = read_advance(4)  // "TRO\x02"
   mode = read_advance(1)
   original_size = read_advance(8)
   compressed_size = read_advance(8)
   compressed_ptr = read_advance(compressed_size)

3. Decompress into text_arena_
   text_arena_.resize(original_size)
   LZ4_decompress_safe(compressed_ptr, text_arena_.data(),
                       compressed_size, original_size)

4. Parse from text_arena_ (same parsing logic as v1, new data source)
   ptr = text_arena_.data()
   end = ptr + text_arena_.size()
   arena_safe_advance() lambda → parse count + records

5. Release mmap (data now owned by text_arena_)
   mmap_handle_.reset()
```

**Key invariant:** After decompression, string_views point into `text_arena_` instead of mmap memory. This is identical to the GZIP path — consistent data ownership model.

#### Backward Compatibility (v1 Support)

`loadBinary()` checks format version byte:

```cpp
if (format_version == 0x01) {
    // v1 path: use mmap'ed data directly
    ptr = static_cast<const char*>(mmap_handle_->data);
    end = ptr + mmap_handle_->size;
    // existing parsing logic
} else if (format_version == 0x02) {
    // v2 path: decompress into text_arena_, then parse
    // LZ4 decompression, parse from text_arena_, release mmap
} else {
    throw "Unknown format version"
}
```

**Result:** Old v1 caches load transparently; new caches use v2 compression.

#### Serialization Helper: serializePayload()

New private method to eliminate code duplication between v1 and v2:

```cpp
void SmartStrategy::serializePayload(std::vector<char>& buf) const {
    // Write count + all records
    // Identical to v1 format (mode-dependent: GenomeIndex vs NGSIndex)
    // Reused by both saveBinary() (after compression) and tests
}
```

### Performance Characteristics

#### Size Reduction

| Dataset | v1 (uncompressed) | v2 (LZ4) | Ratio | Relative Size |
|---------|-------------------|----------|-------|---------------|
| 100MB FASTA | ~105MB | ~35MB | 3.0x | 33% |
| 100MB FASTQ | ~140MB | ~42MB | 3.3x | 30% |

#### Decompression Overhead

```
LZ4 decompression timing (100MB → 35MB compressed):
  - Read mmap header: ~0.00s
  - LZ4_decompress_safe(): ~0.02s (1.75 GB/s throughput)
  - Parse records: ~0.001s (same as v1)
  - Total v2 restore: ~0.021s
  
v1 restore (mmap baseline):
  - Read mmap header: ~0.00s
  - Parse records: ~0.001s
  - Total: ~0.001s

v2 overhead: +0.02s (acceptable tradeoff for 3x file size reduction)
```

#### Lookup Performance

**Zero impact:** Once decompressed and parsed, string_views point into `text_arena_` (identical to GZIP load). Lookup throughput unchanged: **12–18M OPS/s (WGS)**, **28–81M OPS/s (long-read)**.

### Trade-offs

#### Positive Consequences ✅

1. **Storage efficiency:** 3x smaller binary cache files
2. **I/O efficiency:** Faster file transfer (35MB vs 105MB over network)
3. **Decompression speed:** LZ4 much faster than any other algorithm
4. **Backward compatible:** Old v1 caches still work (auto-detection)
5. **No API changes:** `save()` / `restore()` work transparently

#### Negative Consequences ❌

1. **Decompression latency:** +0.02s on restore (rarely significant, cached by downstream tools)
2. **CPU usage:** Minimal (LZ4 highly optimized)

#### Mitigation

The 0.02s decompression overhead is negligible relative to:
- GZIP load: 0.245s (decompression is small fraction)
- Typical bioinformatics workflow: downstream analysis dominates
- Binary cache goal: amortize load cost (load once, restore many times)

### Testing

New test: `TEST_CASE("LZ4 binary cache compression and round-trip integrity", "[cache]")`

```cpp
1. Create FASTA with 3200 nucleotides (highly repetitive)
2. Load into cache (GenomeIndex)
3. Save to v2 binary (LZ4-compressed)
4. Restore from v2 binary
5. Verify all sequences match byte-for-byte
6. Verify compressed_size < uncompressed_size
7. Check round-trip integrity (write → read → compare)
```

---

## Future Work

### ✓ v1.1.0 "Bakuya" (Q1 2026) — Completed

**zlib-ng Integration + Pre-Size + Direct-Write** *(shipped)*
- zlib-ng v2.2.2 via FetchContent, zlib-compat mode
- Pre-size with 3x heuristic + OOM guard (25% available memory)
- Direct-write to text_arena_ (eliminated temp_buffer)
- Result: **0.251s** GZIP decompress (10% improvement over v1.0.0)

### ✓ v1.2.0 "Caladbolg" (Q2 2026) — Completed

**SIMD-Accelerated Parsing** *(shipped)*
- `simd_find_char()` with AVX2 (32 bytes/iter), NEON (16 bytes/iter), runtime dispatch
- Integrated into FASTA/FASTQ single-threaded and multithreaded parsers
- Biggest win in `normalizeFastaArena()` where multi-line FASTA compaction was previously byte-by-byte

**ankerl::unordered_dense Hash Map** *(shipped)*
- Replaced `robin_hood::unordered_flat_map` with Swiss-table design
- Pre-reserved thread-local maps in multithreaded parsers

**Result:** **0.245s** load+parse for 100MB GZIP (**86% reduction** from v1.1.0)

### v1.3.0 "Hrunting" (Q3 2026) — In Progress

**Binary Cache Compression** *(in progress)*
- ✓ LZ4 compression for `.traceon` files (v2 format, `"TRO\x02"`)
- ✓ Format version detection (backward compat with v1)
- ✓ 3x size reduction achieved (105MB → 35MB for 100MB FASTA)
- ✓ Decompression overhead: +0.02s (negligible)
- Status: Implemented and tested

**Parallel GZIP Decompression** *(planned)*
- Integrate Intel ISA-L or custom parallel decompressor
- Target: 4x speedup (0.245s → ~0.06s for 100MB)
- Close remaining ~9% gap with SeqKit (0.245s vs 0.22s)

**Adaptive Chunking** *(future)*
- Dynamic chunk size based on compression ratio
- Optimize for different file characteristics

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

**Status:** ✅ Accepted & Evolving (v1.2.0 base, v1.3.0 binary cache compression added)  
**Last Updated:** 2026-06-23  
**Version:** 1.3.0 "Hrunting" (Binary Cache Compression)

*"Trace On" - Decompressing genomic data at the speed of thought.*