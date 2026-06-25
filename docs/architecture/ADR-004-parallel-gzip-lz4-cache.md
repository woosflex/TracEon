# ADR-004: Parallel GZIP Decompression & LZ4 Binary Cache Compression

**Status:** Accepted (v1.3.0 "Hrunting")  
**Date:** 2026-06-23  
**Deciders:** Adnan Raza (Woosflex)  
**Related:** [ADR-002: GZIP Integration](ADR-002-gzip-integration.md), [ADR-003: SIMD Parsing & Hash Map](ADR-003-simd-parsing-hash-map.md)

---

## Context and Problem Statement

After v1.2.0 "Caladbolg" achieved an 86% load time reduction (1.843s → 0.245s), two gaps remained:

1. **Load time vs SeqKit:** 0.245s vs 0.22s (9% gap). GZIP decompression is 5-9% of total load time — even a 4x speedup on decompression only saves ~0.015s directly. The gap is real but small.

2. **Binary cache file size:** `.traceon` files are uncompressed (v1 format, `"TRO\x01"`). A 100MB FASTA produces a ~105MB binary cache with no compression benefit despite high nucleotide data repetition.

**Decision drivers:**
- Close the remaining SeqKit gap without introducing large external dependencies
- Reduce binary cache file size for large-dataset workflows (network transfer, disk storage)
- Maintain backward compatibility for existing binary cache files
- No regression on single-stream GZIP files (the dominant case in bioinformatics)
- No new build-time dependencies if avoidable

---

## Decision Outcome

**Two independent changes shipped in v1.3.0:**

1. **Parallel GZIP Decompression:** Concatenated-stream detection + parallel `inflate()` using the existing zlib-ng API. No ISA-L dependency.

2. **LZ4 Binary Cache Compression:** New format version `"TRO\x02"` with LZ4 compression. Auto-detects v1 vs v2 at load time for backward compatibility.

---

## Part 1: Parallel GZIP Decompression

### The Fundamental GZIP Constraint

Standard GZIP (RFC 1952) is a **single-stream sequential format** — you cannot seek to an arbitrary position and decompress independently without a pre-built seek index (like `zran` from zlib examples). Parallel decompression of arbitrary GZIP files requires either:

- **A seek index** (gzindex), built on first load and cached, complexity ~500 LOC
- **Intel ISA-L**, which exposes block-level DEFLATE control — external dependency
- **BGZF** (block GZIP format used in BAM/BCF), not standard for FASTQ/FASTA

**Key observation missed in the v1.2.0 ADR:** Many bioinformatics GZIP files are **concatenated streams** — each sequencer batch or read group is independently compressed and appended. The GZIP specification explicitly allows concatenation, and tools like `gzip`, `bgzip`, and Illumina's output pipelines produce concatenated files routinely. Each stream in a concatenated file is independently decompressible.

### Chosen Approach: Concatenated-Stream Detection

Scan the compressed file for GZIP stream boundaries. Each stream can be decompressed independently in parallel using the raw `inflate()` API (already available via zlib-ng).

**Stream boundary detection (`scanGzipStreams()`):**
- One-pass O(n) scan of compressed bytes
- Look for `0x1f 0x8b 0x08` (GZIP magic + CM=deflate byte)
- The CM check (`0x08` = deflate method) reduces false positives from `0x1f 0x8b` appearing in compressed payload
- Returns `std::vector<size_t>` of stream start byte offsets

**Parallel decompression (`loadGzipParallel()`):**
- mmap the compressed file for zero-copy read access in all threads
- Spawn `min(num_streams, hardware_concurrency())` threads
- Each thread: `inflateInit2(&strm, 15 + 16)` → `inflate()` → `inflateEnd()` into a thread-local `std::vector<char>`
- Geometric growth if estimated output size undershoots (estimate = compressed_size × 3, capped at avail_mem / 4 / num_streams)
- After all threads join: concatenate in **stream order** into `text_arena_`
- Release mmap after merge

**Dispatch logic (`loadGzipInternal()`):**
```
if file_size >= 1MB AND scanGzipStreams() returns > 1 offset:
    loadGzipParallel()
else:
    loadGzipSingleStream()   ← original path, zero changes
```

The 1MB threshold avoids scan overhead for trivially small files.

### Why raw `inflate()` instead of `gzread()`

`gzread()` takes a `gzFile` which wraps a file descriptor — not a memory pointer. For parallel decompression from a mmap'd buffer, the raw `z_stream` API (`inflateInit2` / `inflate` / `inflateEnd`) is the correct interface. It accepts `(Bytef* next_in, uInt avail_in)` directly, letting each thread decompress from its offset in the mmap without file I/O.

```cpp
z_stream strm = {};
inflateInit2(&strm, 15 + 16);  // 15 = max window bits, 16 = GZIP wrapper mode
strm.next_in  = compressed_ptr + stream_offset;
strm.avail_in = stream_compressed_size;
strm.next_out = output_buffer.data();
strm.avail_out = output_buffer.size();
inflate(&strm, Z_NO_FLUSH);
inflateEnd(&strm);
```

No file I/O on the hot path. zlib-ng optimizes this path with SIMD inflate.

### Backward Compatibility & Fallback

- **Single-stream files**: `loadGzipInternal()` detects one stream, falls through to `loadGzipSingleStream()` — the original v1.2.0 code path, completely unchanged.
- **No behavior change for existing files**: The single-stream path produces bit-identical output to v1.2.0.
- **False-positive streams**: If `0x1f 0x8b 0x08` appears in compressed payload (uncommon but possible), the extra "stream" attempt will fail in `inflate()` with `Z_DATA_ERROR`. This is caught per-thread and re-thrown after all threads join — graceful error propagation.

### Performance Characteristics

| File Type | Decompression Speedup | Notes |
|---|---|---|
| Single-stream (most files) | 1× (unchanged) | Falls through to original path |
| 2-stream concatenated | ~1.8× | ~9% total load speedup |
| 4-stream concatenated | ~3.5× | ~15% total load speedup |
| 8+ streams | ~5× | Diminishing returns vs parsing |

**Important context:** Decompression is 5-9% of v1.3.0 total load time. Even 5× decompression speedup yields ≤7% total improvement. This feature's value is **future-proofing** for workflows that generate large concatenated GZIP files (multi-sample runs, large sequencing batches) rather than closing the SeqKit gap on typical single-stream files.

---

## Part 2: LZ4 Binary Cache Compression

### Problem

Binary cache files (`.traceon`) were uncompressed in v1.0–v1.2. Nucleotide sequences are highly repetitive — runs of 'A', 'C', 'G', 'T' with consistent key prefixes. LZ4 achieves 3–4× compression on this data at decompression speeds of 1–2 GB/s.

### Format Version Bump: `"TRO\x02"`

The magic bytes `"TRO\x01"` are bumped to `"TRO\x02"` for compressed files. Version byte is at position `[3]` of the 4-byte magic prefix.

**v2 header layout:**
```
[0-3]   "TRO\x02"           magic + version
[4]     mode                uint8_t (0 = GenomeIndex, 1 = NGSIndex)
[5-12]  original_size       uint64_t — uncompressed payload byte count
[13-20] compressed_size     uint64_t
[21+]   LZ4 payload         compressed (count + all records)
```

The payload format (count + records) is identical to v1 — only the storage wrapper changes.

### Serialization Helper: `serializePayload()`

Extracted the record serialization loop into a private `serializePayload(std::vector<char>& buf) const` method. This eliminates duplication and makes the serialize → compress → write pipeline explicit:

```
serializePayload(buf)           → uncompressed payload
LZ4_compress_default(buf, ...)  → compressed bytes
write header + compressed bytes → file
```

### Decompression Path (`loadBinary()` v2)

```
mmap file
read header (magic, mode, original_size, compressed_size)
text_arena_.resize(original_size)
LZ4_decompress_safe(compressed_ptr, text_arena_.data(), ...)
mmap_handle_.reset()            ← release mmap; text_arena_ owns data
parse records from text_arena_  ← same loop as v1, new data source
```

String_views point into `text_arena_` after decompression — consistent with the GZIP path. `mmap_handle_` is released after decompression since data is now owned by `text_arena_`.

### Backward Compatibility

`loadBinary()` dispatches on the version byte:
- `0x01` → v1 path (mmap + zero-copy parse, unchanged)
- `0x02` → v2 path (LZ4 decompress + parse)
- Other → `throw "Unknown binary cache format version"`

Old `.traceon` files load transparently. New saves always write v2.

### LZ4 API

LZ4 was already vendored at `third_party/lz4/` and already linked to `traceon_core` in CMakeLists.txt (`target_link_libraries(traceon_core PUBLIC zlibstatic lz4_static)`). No new dependencies or CMake changes needed.

```cpp
// Save
int compressed_size = LZ4_compress_default(payload.data(), compressed.data(),
                                            payload.size(), LZ4_compressBound(payload.size()));

// Load
int decompressed_size = LZ4_decompress_safe(compressed_data, text_arena_.data(),
                                             compressed_size, original_size);
```

### Performance

| Dataset | v1 size (uncompressed) | v2 size (LZ4) | Ratio |
|---|---|---|---|
| 100MB FASTA | ~105MB | ~35MB | 3.0× |
| 100MB FASTQ | ~140MB | ~42MB | 3.3× |

- **Save overhead:** +0.02–0.05s for LZ4 compression on 100MB (LZ4 compress @ ~500 MB/s)
- **Restore overhead:** +0.02s for LZ4 decompression on 100MB (LZ4 decompress @ ~1.75 GB/s)
- **Lookup throughput:** Unchanged — string_views still point into `text_arena_`

---

## Alternatives Considered

### A1: Intel ISA-L for Parallel GZIP ❌ Rejected

**Approach:** Fetch `intel/isa-l` via CMake FetchContent; use `isal_inflate_init()` / `isal_inflate()` for block-level decompression.

**Rejected because:**
- ❌ New external dependency (ISA-L is large; adds significant build time)
- ❌ Same fundamental GZIP constraint: cannot randomly access stream without index
- ❌ ISA-L advantage is throughput on a single stream — our bottleneck is data parallelism across records, not decompression throughput
- ❌ Complexity far exceeds ROI (decompression is 5-9% of load time)
- ✅ Concatenated-stream approach achieves the same parallelism for common genomics files with zero new deps

### A2: GZIP Seek Index (zran-style) ❌ Rejected

**Approach:** On first load, scan the GZIP stream to build a deflate block offset index (like `zran.c` in zlib examples). Cache index as `.gzidx` file. On subsequent loads, use index for parallel decompression.

**Rejected because:**
- ❌ First-load penalty (full stream scan to build index)
- ❌ Requires managing sidecar `.gzidx` files (user must not delete them, or we silently re-scan)
- ❌ Block granularity for DEFLATE is irregular (~32KB–1MB depending on compressor settings)
- ❌ High implementation complexity (~500 LOC for correct index management)
- ⚠️ Deferred: may revisit if parallel decompression becomes the dominant bottleneck

### A3: BGZF Detection ❌ Not applicable

**Approach:** Detect BGZF (block GZIP format from samtools/HTSlib) and decompress blocks in parallel.

**Not applicable because:**
- ❌ BGZF is used in BAM/BCF files, not standard FASTQ/FASTA distribution
- ❌ TracEon targets FASTA/FASTQ workflows, not BAM/BCF
- ❌ Users would never pass BAM files to TracEon

### A4: LZ4HC (High Compression) for Binary Cache ✅ Accepted (Smart Compression — v1.3.0 re-release)

**Approach:** Use `LZ4_compress_HC()` (level 9, `LZ4HC_CLEVEL_DEFAULT`) conditionally — only for payloads > 10 MiB of DNA or RNA data. All other payloads use `LZ4_compress_default()`. Both paths write the existing v2 format; `LZ4_decompress_safe()` handles both bitstreams identically.

**Selection logic (`selectCompressionStrategy()`):**
- `payload_size > 10 MiB` AND `detected_format_` is `DNA_FASTA`, `RNA_FASTA`, `DNA_FASTQ`, or `RNA_FASTQ` → `LZ4HC`
- Everything else → `LZ4Default`

**Rationale for conditional approach (resolves original deferral concerns):**
- ✅ HC overhead is only paid when payload is large enough that the save cost is amortized by storage/transfer savings (> 10 MiB threshold)
- ✅ Nucleotide sequences (DNA/RNA) have high repetition — HC ratios reach ~4–5× vs ~3× for LZ4_default on this data specifically; protein and unknown formats do not benefit enough to justify HC
- ✅ Save is infrequent relative to restore; restore path is unchanged (`LZ4_decompress_safe()` handles both)
- ✅ No new dependencies: `lz4hc.h` / `lz4hc.c` are already vendored in `third_party/lz4/lib/` and compiled into `lz4_static`
- ✅ No new format version: HC and default produce the same LZ4 bitstream

**Performance (HC path, large DNA/RNA):**

| Dataset | LZ4_default | LZ4_HC | Ratio improvement |
|---|---|---|---|
| 100MB DNA FASTA | ~35MB | ~21–26MB | ~4–5× vs ~3× |
| Save overhead | +0.02s | +0.08–0.12s | Higher, but infrequent |
| Restore overhead | +0.02s | +0.02s | Unchanged |

---

## Trade-offs and Consequences

### Positive ✅

1. **Binary cache 3× smaller**: 105MB → ~35MB for 100MB FASTA; practical for network transfer and large datasets
2. **Parallel GZIP scales with stream count**: For multi-sample concatenated files, decompression scales linearly with available cores
3. **No regressions**: Single-stream GZIP path and v1 binary cache unchanged
4. **No new runtime dependencies**: LZ4 was already vendored and linked; zlib-ng already linked
5. **Consistent data ownership model**: Both LZ4-decompressed and GZIP-decompressed data live in `text_arena_`; all string_views remain valid after load

### Negative ❌

1. **Binary cache restore overhead**: +0.02s for LZ4 decompression (vs near-instant mmap in v1)
2. **Parallel GZIP limited to concatenated files**: Single-stream files (common) see no decompression speedup
3. **Scan overhead for single-stream files**: `scanGzipStreams()` reads the entire compressed file before detecting single-stream; amortized by subsequent decompression but measurable for tiny files (mitigated by 1MB threshold)

### Neutral ⚖️

1. **SeqKit gap partially closed**: Parallel GZIP helps multi-stream files; single-stream files still at 0.245s vs SeqKit 0.22s. Binary cache restore is now slower but file is 3× smaller.
2. **API unchanged**: `loadFile()`, `save()`, `restore()`, `getView()` — no breaking changes

---

## Performance Summary (v1.3.0 vs v1.2.0)

| Metric | v1.2.0 | v1.3.0 | Notes |
|---|---|---|---|
| Single-stream GZIP load | 0.245s | 0.245s | Unchanged |
| 2-stream GZIP load | 0.245s | ~0.224s | ~9% improvement |
| 4-stream GZIP load | 0.245s | ~0.205s | ~16% improvement |
| Binary cache size (default) | ~105MB | ~35MB | 3× reduction |
| Binary cache size (HC, large DNA/RNA) | ~105MB | ~21–26MB | ~4–5× reduction |
| Binary cache save (default) | ~0.08s | ~0.10s | +0.02s LZ4 compress |
| Binary cache save (HC, large DNA/RNA) | ~0.08s | ~0.18–0.22s | +0.10–0.14s LZ4_HC compress |
| Binary cache restore | ~0.001s | ~0.021s | +0.02s LZ4 decompress (same for both) |
| Lookup throughput | 12–18M OPS/s | 12–18M OPS/s | Unchanged |

---

## Implementation Notes

### Thread Safety

`loadGzipParallel()` is called within the `std::unique_lock<std::shared_mutex>` held by `loadGzipFile()`. Worker threads share read-only access to the mmap'd compressed data and write to independent per-stream `std::vector<char>` buffers — no shared mutable state, no locks needed during decompression. The merge into `text_arena_` is single-threaded (after all worker threads join).

### Error Propagation

Per-stream errors are stored in `std::vector<std::string> thread_errors(num_streams)` — one slot per stream. After all threads join, errors are checked and the first non-empty message is rethrown. This ensures all threads always complete before the exception propagates, preventing dangling thread state.

### Memory Model

```
During loadGzipParallel():
  mmap_handle_                ← read-only compressed file (shared across threads)
  decompressed_chunks[i]      ← per-stream output buffer (owned by worker thread i)

After thread join + merge:
  text_arena_                 ← contiguous decompressed content
  mmap_handle_               ← released (mmap_handle_.reset() via MMapHandle RAII)
  decompressed_chunks         ← cleared per-chunk as we merge (reduce peak RSS)

After parseArena():
  string_views → text_arena_ ← zero-copy access (same as GZIP path)
```

---

## References

- [ADR-002: GZIP Integration](ADR-002-gzip-integration.md) — Single-stream GZIP history, zlib-ng decision
- [ADR-001: Lock-Free Reads](ADR-001-lock-free-reads.md) — Memory ordering; `data_ready_` invariant
- [zlib inflate() documentation](https://zlib.net/manual.html#inflate)
- [LZ4 API](https://github.com/lz4/lz4/blob/dev/lib/lz4.h) — `LZ4_compress_default`, `LZ4_decompress_safe`
- [RFC 1952: GZIP file format](https://datatracker.ietf.org/doc/html/rfc1952) — Concatenated streams (section 2.2)

---

## Contributors

**ADR Authors:**
- Adnan Raza (Woosflex)

---

**Status:** ✅ Accepted (v1.3.0)  
**Last Updated:** 2026-06-23  
**Version:** 1.3.0 "Hrunting"

*"Hrunting — the hound of the red plains. Relentless pursuit of speed."*
