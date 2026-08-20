# Changelog

All notable changes to TracEon will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.3.0] - 2026-08-20

**Codename: Harpe** — remote read access + multi-core read scaling.

### Added

- **Remote read access over a minimal TCP protocol** (Slice 1): the v4 in-memory
  cache can now be served over the network and queried by key with zero-copy,
  lock-free reads preserved on the server. New POSIX-only components:
  - `TraceonProto`: length-prefixed framing (64 MiB cap), `HELLO`/`GET`/`HAS`/`STATS`/`BYE`
    messages, and CRC32C-verified `OK`/`OK_STATS` responses (the same Castagnoli
    checksum family as the v4 cache format).
  - `TraceonServer`: thread-per-connection server on an immutable cache
    (graceful shutdown, ephemeral port, `alignas(64)`-isolated read atomics).
  - `TraceonClient`: CRC-verified `getView()`/`has()`/`stats()`, thread-safe.
  - `remote_bench`: `local` / `serve` / `remote` modes; a `docker/` testbed
    (Dockerfile + compose + run.sh) demonstrates a server/client pair over a
    Docker bridge network.
  - ADR-006 documents the protocol and design; the protocol is trusted-network-only
    (no authentication — do not expose to untrusted networks).
- **Multi-core read-scaling benchmark** (Slice 2): `benchmarks/read_scaling.cpp`
  measures lock-free `getView()` throughput at 1/2/4/8/14 threads. On the
  14-core Ultra 5 125H this scales near-linearly through ~8 threads then plateaus
  (~6–7× at 8T, ~7.6–9.4× at 14T, NOT 14×) due to the P/E-core hybrid and shared
  L3 — reported honestly, no manufactured 64-core claim.
- **False-sharing audit**: the hot read-path atomics (`data_ready_`,
  `data_loaded_`, and the `TRACEON_DEBUG_LIFECYCLE` `active_readers_`) are now
  `alignas(64)`-pinned to their own cache lines, separated from `cache_mutex_`.
  Preventive on the default immutable path (RMW-free reads); fixes real
  contention (1.2–1.36× at 4/8/14T) in the opt-in debug diagnostics.
- **kmer C API `kmerindex_reserve()` now rejects absurd sizes deterministically**: ankerl
  `unordered_dense` 4.4.0 `reserve()` clamps to `max_size()` instead of throwing, so
  `reserve(2^62)` silently succeeded (returned 1) on machines where the huge allocation
  landed. The C boundary now pre-checks `n > map.max_size()` (~2^32 entries) and rejects
  such requests with return 0 + diagnostic before any allocation; the no-throw/return-0
  contract is now platform-independent.

### Added

- **Public `Cache` facade completed to match the documented API surface**: the
  README documents `cache.loadGzipFile(...)` and the lifecycle contract names
  `clearCache()` as a reload path, but the high-level `Cache` class did not
  expose them (only `SmartStrategy` did) — the documented examples did not
  compile. `Cache` now forwards `loadGzipFile()`, `clearCache()`, `getAllKeys()`,
  `getQuality()`, and `getDetectedFormat()` to the strategy, with unit tests
  (gzip round-trip, non-gzip rejection + failure atomicity, clear/reload cycle,
  key enumeration, quality round-trip, format detection).

### Fixed

- **`SmartStrategy::loadGzipFile()` now enforces its documented "Throws if file
  is not valid GZIP" contract.** zlib's `gzread()` transparently reads non-gzip
  files as raw data (no error), so an invalid file previously loaded as plain
  text. The loader now validates the gzip magic bytes (`0x1f 0x8b`) up front —
  the same check `loadFile()` already performs — and throws
  `std::runtime_error` otherwise. Failure atomicity is now an explicit,
  tested contract (see `SmartStrategy::loadGzipFile()` docs): pre-flight
  rejections (unreadable / not GZIP) throw BEFORE teardown and preserve the
  loaded snapshot; mid-stream failures (truncated/corrupt stream, trailing
  garbage, OOM guard) throw after teardown has begun and leave the cache
  empty — never partial (pre-loaded-cache tests added for both paths).

### Research

- Surveyed khash/khashl and structurally similar minimizer indexes across Winnowmap,
  minigraph/minigraph-cactus, mm2-fast, BLEND, StrobeMap, strobealign,
  HiFiMapper, LRA, and meryl. The canonical report is
  `outputs/aligner-khash-traceon-survey.md`; it contains no benchmark results.

- Investigated a TracEon `traceon_kmer` C API drop-in for minimap2's per-bucket
  `khash_t(idx)` table. Evidence and benchmark recommendations are recorded in
  `outputs/minimap2-traceon-integration-research.md`; no source files were
  modified. Current experimental patch remains blocked on iterator ABI and
  error-status handling fixes.

---

## [2.3.0] — v2.3.0 "Harpe" (unreleased)

*Harpe slice 1 — remote read access.* Release stamp/date land at release.

### Added

- **Remote read access over TCP** (protocol, server, client):
  - `include/TraceonProto.h` — a minimal length-prefixed wire protocol
    (8-byte LE payload_len prefix, 64 MiB cap, message-type bytes) with pure,
    testable encode/decode codecs and CRC32C-verified OK responses.
  - `TraceonServer` (`include/TraceonServer.h` + `src/TraceonServer.cpp`) —
    thread-per-connection server over an immutable loaded cache, ephemeral
    port support, idempotent `stop()`/restart lifecycle.
  - `TraceonClient` (`include/TraceonClient.h` + `src/TraceonClient.cpp`) —
    CRC32C-verified, thread-safe client (`getView`/`has`/`stats`/`close`),
    hostname resolution, fast-failing connect.
  - `benchmarks/remote_bench.cpp` — three-mode benchmark (`local`/`remote`/
    `serve`), one connection per thread, median/p95/p99 latency + ops/s.
  - Protocol/architecture decision record:
    `docs/architecture/ADR-006-traceon-remote-access.md`.
- **Multi-core read scaling (slice 2)**:
  - `benchmarks/read_scaling.cpp` — lock-free `getView()` read throughput vs
    thread count (1/2/4/8/14/18) against one in-memory immutable cache built
    via the `addEntry()` path; standalone benchmark target (not in the default
    ctest suite), per-thread shuffled key traversal + start barrier, aggregate
    ops/s + scaling factor + median per-thread latency + hit-rate/checksum
    verification. Honest single-socket framing: near-linear to 8 threads
    (~7×), plateau to ~9.5× at 14 threads on the Ultra 5 125H (P/E-core hybrid,
    shared 18 MiB L3, all-core frequency drop) — no 64-core story.
  - **False-sharing audit** (`include/SmartStrategy.h`): pinned the hot
    read-path atomics (`data_ready_`, `data_loaded_`, and opt-in
    `active_readers_` under `TRACEON_DEBUG_LIFECYCLE`) to distinct cache lines
    with `alignas(64)` and separated them from `cache_mutex_`. The default read
    path was already RMW-free (pure acquire loads), so this is preventive there
    and a real fix for the writer-flip/teardown and the debug diagnostic (whose
    per-lookup `active_readers_` RMW previously shared a line with
    `data_ready_`, measured ~1.2–1.36× cost at 4/8/14 threads).
- **Docker testbed** (`docker/`): multi-stage image building the server +
  tools (Release, `TRACEON_BUILD_SERVER=ON`), compose services
  (make-cache → server → client) on a bridge network, and `run.sh` for
  one-shot verification. Added `docker/data/sample.fasta` +
  `.dockerignore`.
- **`tools/traceon_make_cache`**: converts FASTA/FASTQ into a v4 `.traceon`
  cache for server deployment.

### Fixed

- **`TraceonServer::stop()` deadlock**: closing the wake-pipe/listen fds while
  the accept thread is blocked in `poll()` is undefined behaviour on Linux and
  did not reliably wake `poll()`, so `stop()` could join-hang forever. The
  server now wakes the accept loop via the pipe write, joins it first, then
  closes the fds; the accept loop also drops (rather than spawns) a
  connection accepted mid-shutdown.
- **`TraceonClient` hostname support**: only IPv4 literals were accepted
  (`inet_pton`), which broke containerized deployments reaching a peer by
  service name. The client now resolves hosts via `getaddrinfo`.
- **`remote_bench` CLI**: documented `--opt=value` form was rejected by the
  exact-string parser; the parser now accepts both `--opt value` and
  `--opt=value`.



## [2.2.0] — v2.2.0 "Gáe Bolg" (2026-08-14)

### Added
- **Aligner validation matrix COMPLETE** — 5 real tools with byte-identical outputs:
  `woosflex/minimap2` (−13.6% full chr1), `woosflex/Winnowmap` (−80% chr1-50Mb),
  `woosflex/mm2-fast` (−39%/−33% RSS), `woosflex/minigraph` (first khashl adapter),
  `woosflex/BLEND` (fuzzy-seed mapper, 7-config identity). Every fork: stock
  byte-identical + opt-in `TRACEON=1` table backend + `tcache` mmap'd prebuilt index.
  See `RELEASE_v2.2.0.md` and `TracEon - Aligner Integration Results.md` (vault).
- Documented upstream bugs found during integration: Winnowmap applyWeight crash +
  uninitialized rep_len; minigraph GCC-16 crash on its own test data (memcpy-load
  fix in port) + gi->g_own leak; mm2-fast safestringlib GCC-16 CFLAGS patch.

### Changed
- (none — no core/API changes in v2.2.0; additive release)

## [2.1.0] — v2.1.0 "Excalibur" (2026-08-12)

### Added
- **Fuzzing gate**: 7 libFuzzer targets (v4 loader, gzip loader, FASTQ/FASTA parsers, kmer encode, kmer C API, TRKI mmap loader), 83-file seed corpus, `fuzz.yml` CI workflow (per-push + weekly, crash artifacts). Closes the last never-run verification layer.
- **Aligner validation**: `woosflex/minimap2` and `woosflex/Winnowmap` `traceon-backend` forks — `TRACEON=1` table backend + `tcache v2` mmap'd open-addressing prebuilt index. Recurring load+map: minimap2 −13.6% (full chr1), Winnowmap −80% (chr1-50Mb); PAF byte-identical to stock in every mode. See `RELEASE_v2.1.0.md` and `outputs/*-integration-*.md`.

### Changed
- kmer C API: `frozen` flag is now `std::atomic<bool>` (release/acquire) — thread-safe concurrent lookup (the minimap2-family integration requirement). Error path is allocation-free (fixed-size thread-local buffer) — the no-throw C boundary is total under OOM.

### Fixed
- `assert(inserted)` audit in kmer tests/examples — explicit `== 1` checks (the `-1` exception status was a silent trap).
- CI fuzz integration: `extern "C"` linkage on all fuzz targets (libFuzzer link), OSS-Fuzz `fuzzer-no-link` pattern (duplicate ASan runtime on ubuntu clang), 64MiB declared-length harness cap + `allocator_may_return_null` (4GB alloc abort under ASan).

## [2.0.0] — v2.0.0 "Durandal" (2026)

*Integrity, v4 binary format, and the lifecycle contract. Verified:
**128 test cases / 4195 assertions** (was 115 / 4076).*

### 💥 Breaking Changes / Migration

1. **Binary cache:** `.traceon` is **v4 only** (`"TRO\x04"`). v1/v2/v3 caches
   are rejected with "unsupported cache version; regenerate with v2.0.0" —
   there are no legacy readers. Regenerate caches after upgrading.
2. **Input integrity:** truncated/corrupt GZIP and trailing garbage after the
   last GZIP member now throw `std::runtime_error`; partial data is never
   published.
3. **Mutation/lifecycle:** `set()`/`addEntry()` after a load throws
   `std::logic_error`; call `clearCache()` only after all readers stop using
   views.
4. **View/storage:** `getView()` remains zero-copy and non-owning; its result
   must not be used after `clearCache()`, reload, or destruction (reader
   quiescence contract — see ADR-001).
5. **Source/ABI:** `std::string_view`-backed index keys require recompilation;
   **no ABI stability is promised** across 1.x → 2.x.
6. **k-mer surface:** the k-mer C API ships as an opt-in separate target (`traceon_kmer`), hardened: noexcept C boundary with status codes + `kmerindex_last_error()`, caller-owned iterators, freeze semantics, validated TRKI loading. Not part of the core umbrella.

### ✨ Added

- **`.traceon` v4 binary format**: `saveBinary()`/`Cache::save()` now write
  `"TRO\x04"` — a 26-byte header (magic, codec flags, index mode, logical
  payload length u64-LE, compressed frame length u64-LE, CRC32C u32-LE)
  followed by a streamed LZ4 Frame. The CRC32C checksum covers the ENTIRE
  uncompressed logical payload plus the canonical header fields, computed
  **incrementally** as chunks pass through the LZ4F compressor/decompressor
  (no second full-payload pass, no extra allocation).
- **CRC32C with runtime dispatch**: x86 SSE4.2 `crc32` instruction (guarded by
  `TRACEON_HAS_AVX2`, selected via `__builtin_cpu_supports("sse4.2")`),
  AArch64 `crc32cx/crc32cw` intrinsics (guarded by `TRACEON_HAS_NEON` /
  `__ARM_FEATURE_CRC32`), and a portable table fallback — with a tested
  fallback-vs-hardware equivalence guarantee (new `include/Crc32c.h`).
- **Lifecycle contract documentation**: reader quiescence is now prominent in
  `include/Cache.h`, `include/SmartStrategy.h`, `README.md`, and ADR-001;
  `TRACEON_DEBUG_LIFECYCLE` adds a debug-only reader counter/scope that warns
  in `clearInternal()` when a `getView()` lookup overlaps a teardown
  (diagnostic only — not synchronization).

### 🔧 Changed

- **v4 reader hardening (semantics preserved from v2/v3)**: implausible sizes
  rejected before allocation (undersized logical length, 32-bit `size_t`
  guard, OOM guard + `std::bad_alloc` catch), subtraction-form bounds checks
  (`n > end - ptr`), exact-length + exact-frame-termination requirements,
  count-bounded reserve (`remaining_bytes / min-record-bytes`), and
  **failure atomicity** — an invalid load never publishes `data_ready_` or
  partial map state (cache resets to the pristine empty state).
- **Wrong-field rejection**: unknown magic/version, unsupported codec flags,
  invalid index mode, wrong logical length, wrong record count, and
  checksum mismatches each throw a specific `std::runtime_error`.

### 🐛 Fixed

- **FASTQ empty-line handling**: Empty quality lines now parse as `qual=''`
  (previously they corrupted record pairing and dropped every subsequent
  record); empty sequence lines no longer break the 4-line cycle. Both parsers
  consume a strict positional 4-line cycle (header / seq / `+` / qual)
  advancing **exactly one line terminator per boundary** at all 8 sites.
- **FASTA header-only / abutting records**: `normalizeFastaArena()` writes the
  record terminator only when the sequence actually produced bytes.
  Header-only records (`>a\n>b\n...`) and records whose sequence abuts the
  next `>` without a newline (`>a\nACGT>b\n...`) are preserved — previously
  the unconditional terminator clobbered the next record's `>` marker and
  silently truncated the whole file.
- **Multithreaded chunk-boundary classifier**: the backward lookbehind was
  replaced with a **forward** check — a candidate line-start `@` is a genuine
  header iff the line two lines later begins with `+`. Quality lines that start
  with `+` (Phred+33 Q10, ASCII 43) no longer cause every worker after chunk 0
  to drop its whole chunk (repro: 80k reads → 4,468 loaded, 94.4% loss).
- **GZIP integrity validation**: `loadFile()` / `loadGzipFile()` / `restore()`
  now throw `std::runtime_error` on truncated GZIP streams, trailing garbage
  after the last member, and truncated tail members in concatenated GZIP —
  corrupt/truncated `.gz` files that previously loaded silently now fail.
- **Binary-cache hardening**: v2 `original_size` / `compressed_size` >
  `INT_MAX` rejected before allocation; v2 resize carries the v3-style OOM
  guard + `std::bad_alloc` catch; all pointer bounds checks use the
  overflow-safe subtraction form; record counts bounded by remaining payload
  bytes (≥ 12 B/record GENOME, ≥ 20 B NGS) in v1/v2/v3 — crafted caches can
  no longer force multi-GB allocations.
- **Duplicate-key `set()` save/restore consistency**: `addEntry()` incremented
  `serialized_size_estimate_` on a no-op `emplace()`, over-declaring the
  binary header so `restore()` rejected the cache; now first-wins with a count
  precheck.
- **Immutable-after-load enforcement**: `set()`/`addEntry()` after a load
  throws `std::logic_error` (see ADR-001); `clearCache()` reopens the build
  phase.

### 🧪 Testing

- **128 test cases / 4195 assertions** — v1/v2/v3 compatibility tests replaced
  with v4 equivalents: GENOME + NGS round-trips, CRC32C known-answer
  (`0xE3069283`), single-byte mutations at header/middle/tail → checksum
  error, truncation at frame boundaries, wrong logical length, wrong record
  count, wrong magic/codec/mode → rejection, legacy v1/v2/v3 rejection, and
  fallback-vs-hardware CRC equivalence. The v2 `2^32+N` truncation probe,
  count-bomb, and implausible-size hardening tests were ported to v4
  semantics. All parser/GZIP/integrity tests remain green.

## [1.4.0] — v1.4.0 "Caliburn" (2026-06-30)

*The sword of selection — choosing the right index mode.*

### ✨ Added

- **Streaming v3 Binary Cache Format**: `saveBinary()` now writes `.traceon` v3 (`"TRO\x03"`) — a streamed LZ4 Frame instead of a single LZ4 block over the full serialized payload. `serializePayload()` pushes bytes through a 1 MB streaming window, so peak save/load memory stays bounded by `STREAM_CHUNK_SIZE` regardless of dataset size (the old v2 path held text_arena_ + payload + compressed copies all at once). `loadBinary()` still reads v1 (uncompressed) and v2 (LZ4 block) for backward compatibility.
- **IndexMode Selection**: `SmartStrategy(IndexMode mode)` and `Cache(IndexMode mode)` constructors allow choosing between `IndexMode::GENOME` (default, string-keyed `GenomeIndex`) and `IndexMode::NGS` (hash-keyed `NGSIndex`). Previously, `NGSIndex` was compiled in but permanently unreachable.
- **`Cache::getIndexMode()`**: New method exposes the active index mode from the public `Cache` API.
- **`Cache::set()` persists via `save()`**: `set()` now routes through `SmartStrategy::addEntry()` instead of a separate `m_manual_store`. Manually-added entries are now serialized by `save()` and restored by `restore()`. The separate `m_manual_store` field in `Cache` has been removed.

### 🔧 Changed

- **`SmartStrategy` constructor signature**: `SmartStrategy()` → `SmartStrategy(IndexMode mode = IndexMode::GENOME)`. Default behavior is identical.
- **`Cache` constructor signature**: `Cache()` → `Cache(IndexMode mode = IndexMode::GENOME)`. Default behavior is identical.
- **`Cache::get()` / `getView()` / `hasSequence()` / `size()`**: Removed the `m_manual_store` lookup path; all data now lives in a single `SmartStrategy`-managed store.
- **`FastqTests.cpp`**: Replaced fragile relative path (`../test_data/simple.fastq`) with a self-generated temp file. Tests now pass regardless of working directory.

### 🧪 Testing

- **+12 new test cases** (57 → 69), **+56 new assertions** (3776 → 3832)
- **NGSIndex correctness**: FASTA load/lookup, FASTQ load/lookup, save/restore round-trip via NGSIndex
- **Real parallel GZIP coverage**: Two concatenated streams totalling >1 MB (1200 records × 500-char sequences each), verifying `loadGzipParallel` is actually exercised (previous test used tiny data below the 1 MB threshold and fell through to single-stream path)
- **RNA FASTA format detection**: `FileFormat::RNA_FASTA` triggered by U-containing sequences
- **Protein FASTA format detection**: `FileFormat::PROTEIN_FASTA` triggered by protein amino acid sequences
- **`clearCache()` + reload**: Verifies no dangling string_view into old `text_arena_` after clear and re-load
- **v1 binary format round-trip**: In-memory construction of a `TRO\x01` blob, written to disk and loaded via `loadBinary()`, verifying backward-compatibility of the v1 format
- **`Cache::getView()` coverage**: Basic zero-copy view test via public Cache API
- **`Cache::set()` + `save()` + `restore()`**: Verifies manually-added entries survive a save/restore cycle

---

## [1.3.0] — v1.3.0 "Hrunting" (2026-06-23)

### ✨ Added

- **LZ4 Binary Cache Compression**: `.traceon` binary cache files now use LZ4 compression (format version 0x02). Decompresses at ~1–2 GB/s with ~3x size reduction (~105MB → ~35MB for 100MB FASTA).
- **Format Version Detection**: `loadBinary()` now detects format version byte and supports backward compatibility with v1 uncompressed format (0x01). Automatic routing based on magic bytes `"TRO\x01"` (v1) vs `"TRO\x02"` (v2).
- **Serialization Helper**: New private method `serializePayload()` abstracts payload serialization, eliminating code duplication between v1 and v2 paths.
- **Parallel GZIP Decompression**: Concatenated GZIP streams are now decompressed in parallel. `scanGzipStreams()` scans for GZIP stream boundaries; `loadGzipParallel()` spawns one thread per stream using the zlib-ng `inflate()` API. Single-stream files fall through to the existing single-threaded path unchanged.
- **Smart Compression**: `saveBinary()` now auto-selects the LZ4 compression mode based on payload size and detected file format. Payloads > 10 MiB of DNA or RNA data use `LZ4_compress_HC()` (level 9) for ~4–5× size reduction; all other data uses `LZ4_compress_default()`. Both modes wrote the v2 format — `LZ4_decompress_safe()` handles both bitstreams transparently. **Correction:** the "no format version bump" note below is wrong in hindsight — v1.4.0 replaced the v2 single-block payload with the streaming v3 LZ4-Frame format (see the v1.4.0 entry), and `saveBinary()` has written v3 since.

### 🔧 Changed

- **Binary Cache Format**: Version byte bumped from `\x01` to `\x02`. Header now includes `original_size` and `compressed_size` fields for decompression.
- **Binary Cache I/O**: v2 format decompresses LZ4-compressed payload into `text_arena_` before parsing. String views still point into arena (consistent with GZIP path).
- **GZIP Loader**: `loadGzipInternal()` is now a dispatcher — detects stream count and routes to `loadGzipParallel()` (multi-stream) or `loadGzipSingleStream()` (single-stream, unchanged).

### ⚡ Performance

- **Binary Cache Size**: 3x reduction (e.g., 105MB → ~35MB for 100MB FASTA sequences)
- **Restore Latency**: Unchanged or faster (LZ4 decompression ~1–2 GB/s, smaller file size reduces I/O)
- **Parallel GZIP Decompression**: For concatenated GZIP files (common in bioinformatics sequencer output), decompression scales with stream count (~1.8x for 2 streams, ~3.5x for 4 streams)
- **Binary Cache Size (HC path)**: ~4–5× reduction for large DNA/RNA datasets vs ~3× with LZ4 default (e.g., 100MB DNA FASTA → ~21–26MB with HC vs ~35MB with default)
- **Lookup Throughput**: Unchanged (12–18M OPS/s WGS range, 28–81M OPS/s long-read scenarios)

---

## [1.2.0] — v1.2.0 "Caladbolg" (2026-06-09)

### ✨ Added

- **SIMD Record Boundary Scanning**: `simd_find_char()` with AVX2/NEON/scalar dispatch for accelerated record boundary scanning in FASTA/FASTQ parsers
- **ankerl::unordered_dense Hash Map**: Integrated via CMake FetchContent (header-only, MIT license)
- **Pre-reserved Thread-Local Maps**: Multithreaded parsers now pre-reserve thread-local maps to reduce rehashing

### 🔧 Changed

- **Hash Map Replacement**: Replaced `robin_hood::unordered_flat_map` with `ankerl::unordered_dense::map` for better insert performance
- **TransparentStringHash**: Updated to use `std::hash` instead of `robin_hood::hash`
- **Pre-reservation Heuristic**: Conservative formula — `chunk_size/500` (FASTA), `chunk_size/600` (FASTQ) with 1.25x safety margin

### 🐛 Fixed

- **normalizeFastaArena() trailing newline**: Fixed bug where `resize()` zero-initialized trailing `\n` when input lacked trailing newline

### ⚡ Performance

- **100MB GZIP Load Time**: 0.245s (86% reduction from 1.843s v1.1.0 baseline)
- **Peak RSS**: 185MB (30% reduction from 263MB v1.1.0 baseline)
- **Lookup Throughput**: 12M OPS/s (WGS), 28–81M OPS/s (long-read scenarios)

## [1.1.0] — v1.1.0 "Bakuya" (2026-06-07)

### ✨ Added

- **zlib-ng Integration**: Replaced system zlib with zlib-ng v2.2.2 via CMake FetchContent in zlib-compat mode. Provides ~10-15% faster inflate through SIMD-optimized CRC and decompression routines.
- **Pre-Size with OOM Guard**: Pre-allocate `text_arena_` to `min(compressed_size × 3, available_memory × 0.25)` before decompression. Uses `std::filesystem::file_size()` for accurate compressed size and platform-specific memory queries (`/proc/meminfo`, `sysctl`, `GlobalMemoryStatusEx`). On OOM, throws `std::bad_alloc` with diagnostic.
- **Direct-Write Decompression**: Eliminated `temp_buffer` entirely. Decompressed data written directly into `text_arena_` with geometric growth fallback. Retains 1MB chunk buffer for `gzread` source reads.
- **shrink_to_fit**: Release excess arena capacity after decompression to minimize steady-state memory.

### 🔧 Changed

- **GZIP Load Time**: Reduced from 0.28s to 0.251s for 100MB GZIP (-29ms, 10% improvement)
  - zlib-ng SIMD inflate: -25ms
  - Pre-size (no realloc): -3ms
  - Direct-write (no move): -1ms
- **GZIP Memory Peak**: Reduced from ~366MB to ~266MB for 100MB GZIP (27% reduction) by eliminating `temp_buffer` co-allocation
- **Removed `third_party/zlib/`**: zlib-ng is now fetched at build time via FetchContent
- **Update `third_party/` notes**: zlib-ng no longer vendored in-tree

### 🐛 Fixed

- No bug fixes in this phase — pure performance optimization

### 🧪 Testing

- Validated GZIP load times match v1.1.0 targets (0.251s ± 0.01s)
- Memory profiling confirmed 266MB peak RSS for 132MB compressed file
- All existing GZIP test cases pass (auto-detect, explicit load, magic bytes)

### 📊 Benchmarks (Intel Core Ultra 5 125H, 16GB RAM)

**Load Time Component Breakdown (100MB GZIP, total 0.251s):**
- Decompression (zlib-ng inflate): 5-9% of load time
- FASTA parsing (arena construction): 46% of load time
- Hash map building (insertion): 48% of load time

**Throughput:**
- 294-378 MB/s decompressed throughput across file sizes (10MB–1GB)
- Consistent performance regardless of compression ratio due to direct-write design

---

## [1.0.1] — Patch Release

### 🐛 Fixed (Critical)

- **`addEntry()` dangling string_views (UB)**: `SequenceView` members were views into
  function-parameter strings that went out of scope immediately. Replaced with
  `std::deque<std::string> manual_store_` whose elements have reference stability
  (push_back never invalidates existing references per the C++ standard).

- **`loadBinary()` mmap with no error handling**: `open()` / `fstat()` / `mmap()` failures
  on Linux, and `CreateFileA()` / `MapViewOfFile()` failures on Windows, all went
  unchecked — leading to reads from invalid/NULL pointers (crash/UB). Full error
  propagation added via `std::runtime_error`.

- **`loadBinary()` no bounds checking**: Reading from the mmap buffer advanced `ptr`
  without verifying it stayed within `[data, data+size)`. A corrupt or truncated binary
  cache file would cause out-of-bounds reads. A `safe_advance(n)` lambda now throws
  on any overrun, and a record-count sanity cap (1B) prevents runaway loops.

- **Unaligned reads in `loadBinary()`**: `*reinterpret_cast<const uint32_t*>(ptr)` is
  undefined behaviour on strict-alignment architectures (ARM). All integer reads now use
  `std::memcpy` which the compiler optimises away on x86/x64 while remaining portable.

- **Multi-line FASTA sequences contained embedded newlines**: NCBI-style wrapped FASTA
  (60 chars/line) produced sequence `string_view`s with internal `\n` chars. A new
  `normalizeFastaArena()` pass compacts the arena in-place before any parser runs, so
  every sequence is guaranteed to be a single contiguous string. Works for both
  single-threaded and multithreaded parse paths.

### ✨ Added

- **Heterogeneous lookup on hot read path**: `GenomeIndex` now uses `StringHashMap`
  (from `MapDefs.h`) with a `TransparentStringHash` and `std::equal_to<>`. `getView()`
  can now call `map.find(sequence_id)` with a `string_view` directly — no `std::string`
  heap allocation per lookup. Requires C++20; falls back to `std::unordered_map` without
  robin_hood.

- **Binary cache format versioning**: Magic changed from `"MMAP"` to `"TRO\x01"`.
  Opening an old v1.0 binary cache now throws a descriptive `std::runtime_error` instead
  of silently producing garbage data.

### 🔧 Changed

- **GZIP buffer growth**: Removed the single upfront `3×compressed_size` reservation
  (inaccurate for highly compressible or incompressible data). The decompression loop
  now doubles `temp_buffer` capacity geometrically when needed, matching standard
  container growth behaviour.

- **Thread count guard**: Both multithreaded parser templates now use
  `std::max(1u, std::thread::hardware_concurrency())` to avoid division-by-zero on
  platforms where `hardware_concurrency()` returns 0.

- **`clearInternal()`**: Now also clears `manual_store_`, freeing strings added via
  `addEntry()` when the cache is reloaded or destroyed.

### 🧪 Testing

- 15 new assertions across 7 new test cases: multi-line FASTA normalization,
  `addEntry` stable-view correctness, `loadFile` / `loadBinary` on missing files,
  truncated binary, old-format rejection, initial state validation.
- Total: **43 assertions in 12 test cases** (was 28 in 6).

---

## [1.0.0] - 2025-12-16 "Avalon"

### 🎯 Codename
**Avalon** - Like EMIYA's absolute defense from Fate/stay night, this release provides performance isolation through lock-free concurrency and zero-copy architecture.

### ✨ Added

#### Core Features
- **Native GZIP Support** ([ADR-002](docs/architecture/ADR-002-gzip-integration.md))
  - Automatic format detection via `.gz` extension or magic bytes (`0x1f 0x8b`)
  - Streaming decompression with zlib (1MB chunks)
  - Minimal performance impact: < 7% overhead on lookups after decompression
  - Public APIs: `loadFile()` (auto-detect), `loadGzipFile()` (explicit)

- **Lock-Free Concurrent Reads** ([ADR-001](docs/architecture/ADR-001-lock-free-reads.md))
  - C++20 atomics with `memory_order_acquire`/`release` semantics
  - Immutable-after-load design eliminates all read locks
  - `std::atomic<bool> data_ready_` replaces `std::shared_mutex` on hot path
  - 2x potential throughput improvement on concurrent workloads

- **Robin Hood Flat Hash Maps**
  - Switched from `std::unordered_map` to `robin_hood::unordered_flat_map`
  - Open addressing provides ~15% better cache locality
  - Pre-allocation strategy: Reserve 125% of estimated capacity to prevent rehashing

- **Memory-Mapped Binary Cache**
  - Custom `.traceon` format with magic bytes verification (`MMAP`)
  - Instant restoration via `mmap()` / `MapViewOfFile()` (sub-second for GB files)
  - Cross-platform support (Linux, macOS, Windows)

#### Architecture
- **Zero-Copy Data Access**
  - Single arena allocation (`std::vector<char>` text_arena_)
  - All sequences returned as `std::string_view` (pointer + length)
  - No heap fragmentation

- **Hybrid Key Storage**
  - `std::string` keys (leverages SSO for short IDs < 15 chars)
  - `std::string_view` values (zero-copy)
  - Trade-off: 24 bytes overhead per key for improved cache locality

- **Forced GenomeIndex Mode**
  - Disabled NGSIndex (hash-only) despite being implemented
  - Rationale: SSO cache locality outweighs hash speed for typical genomic IDs

#### Performance Optimizations
- **Multithreaded Parsing**
  - Automatic parallelization for files > 10MB
  - Scales linearly up to 8 cores
  - Thread-local caches merged without lock contention

- **Capacity Pre-Allocation**
  - Parsers reserve 125% of estimated record count
  - Prevents rehashing during insertion
  - Heuristics: ~100 bytes/record (FASTA), ~150 bytes/record (FASTQ)

### 📊 Performance

#### Benchmarks (Intel Core Ultra 5 125H, 16GB RAM)

**Small Datasets (10MB, L3 Cache Resident):**
- WGS Short Reads: 40-55M OPS/s
- PacBio Long Reads: 50-60M OPS/s
- Reference Genomes: 40-50M OPS/s

**Medium Datasets (100MB, RAM Latency):**
- WGS Short Reads: 12-18M OPS/s (plain), 11-17M OPS/s (GZIP)
- PacBio Long Reads: 25-35M OPS/s
- Reference Genomes: 15-25M OPS/s

**Large Datasets (500MB-1GB):**
- WGS Short Reads: 10-15M OPS/s
- PacBio Long Reads: 20-30M OPS/s
- Reference Genomes: 12-20M OPS/s

**Memory Efficiency:**
- 10MB dataset: ~25 MB RSS
- 100MB dataset: ~180 MB RSS (5x more efficient than BioPython)
- 500MB dataset: ~900 MB RSS

**Binary Cache:**
- Save: 0.08s (100MB), 0.25s (1GB)
- Restore: ~0.00s (instant via mmap)

### 🔧 Changed

#### API
- `loadFile()` now auto-detects GZIP format (extension + magic bytes)
- Renamed internal `clear()` → `clearInternal()` to fix deadlock
- Made `getView()` lock-free (removed `std::shared_lock`)

#### Internal
- Replaced `std::shared_mutex` with `std::atomic<bool>` on read path
- GZIP decompression uses streaming with geometric buffer growth
- Optimized chunk size: 1MB for GZIP decompression

### 🐛 Fixed

#### Critical
- **Lock Overhead**: Eliminated 25-35% of execution time on read-heavy workloads
- **Performance Regression**: Fixed debug build masking as release (6M → 40M OPS/s)
- **HashMap Load Factor**: Pre-allocation prevents rehashing (was causing 20% degradation)
- **Deadlock in clearCache()**: Refactored to use internal helper

#### Minor
- Missing error handling for `tellg()` failures in GZIP size estimation
- Race condition in format auto-detection (moved detection outside lock)
- Test path dependencies (tests require running from build directory)

### 📚 Documentation

#### Added
- **ADR-001**: Lock-Free Reads (memory ordering justification)
- **ADR-002**: GZIP Integration (design decisions, rejected alternatives)
- **Performance Profile** with regression thresholds and troubleshooting
- **Comprehensive Benchmark Guide** (`benchmarks/README.md`)
- Single-include header (`TracEon.h`) with Doxygen comments
- Working example code (`examples/simple_usage.cpp`)

#### Updated
- README with GZIP feature highlights and accurate system specs
- Architecture diagrams showing lock-free data flow
- API documentation for `loadFile()` / `loadGzipFile()`

### 🧪 Testing

#### Added
- Lock-free concurrency tests (4-thread validation)
- GZIP test suite (8 assertions):
  - Explicit load (`loadGzipFile()`)
  - Auto-detect via extension (`.gz`)
  - Auto-detect via magic bytes (`0x1f 0x8b`)
- Zero-copy memory semantics verification
- Architecture validation (GenomeIndex usage)

#### Coverage
- 100% of public API methods
- Edge cases: empty files, malformed records, truncated streams
- Concurrency: Multi-threaded read validation

### 🔄 Migrations

**No breaking changes** - v1.0.0 is the initial public release.

### 🏗️ Build System

#### Dependencies
- **Required**: CMake 3.20+, C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- **Bundled**: 
  - robin_hood.h v3.11.5 (MIT License)
  - zlib 1.2.11 (zlib License)
  - lz4 1.9.4 (BSD License)
- **Test**: Catch2 v3.6.0 (fetched via CMake FetchContent)

#### Platforms
- ✅ Linux x86_64 (primary development platform)
- ✅ macOS ARM64 (Apple Silicon M1/M2/M3)
- ⚠️ Windows x64 (experimental, MSVC 2019+ required)
- 🔄 ARM64 Linux (untested but should work)

### ⚠️ Known Limitations

1. **RAM Latency Wall**: 3-4x degradation when dataset exceeds L3 cache (~16-32MB)
   - This is **hardware physics** (100ns RAM vs 10ns cache), not a software bug
   - Mitigation: Binary cache reduces load time, not lookup time

2. **Immutable After Load**: Cannot incrementally update cache
   - Required for lock-free reads
   - Alternative: Reload entire cache (fast with mmap)

3. **Single-Node Only**: No distributed caching
   - Target: Local workstation/server workloads
   - Future: v2.1.0 distributed architecture

4. **Test Path Dependency**: Tests must run from build directory
   - Workaround: `cd build && ./unit_tests`

5. **RecordTypes.h Structure**: `SequenceView` defined in `SmartStrategy.h`
   - Architectural decision for internal implementation details
   - Public API abstracts this via `Cache.h`

---

## Version History

| Version | Date | Codename | Key Features |
|---------|------|----------|--------------|
| 1.0.0 | 2025-12-16 | Avalon | Zero-copy, Lock-free, GZIP support, Memory mapping |
| 1.1.0 | 2026-06-07 | Bakuya | zlib-ng integration, Pre-size + Direct-write, Optimized GZIP |
| 1.2.0 | 2026-06-09 | Caladbolg | SIMD parsing, ankerl::unordered_dense hash map, 86% load time reduction |
| 2.0.0 | 2026-Q3 | Durandal | C API, Streaming |
| 2.1.0 | 2026-Q4 | Excalibur | Distributed caching |

---

## Comparison with State-of-the-Art

### vs BioPython
- **Performance**: 26-54x faster on long reads
- **Memory**: 5x more efficient (180MB vs >900MB for 100MB dataset)
- **Concurrency**: Lock-free vs single-threaded

### vs PyFastX
- **Architecture**: Memory-resident vs disk-cached
- **Long Reads**: 26-54x faster (TracEon excels with sparse indexing)
- **Short Reads**: Competitive (disk caching can be faster for specific workloads)
- **Concurrency**: Native multi-threading vs limited parallelism

---

## References

### External Resources
- [C++ Memory Model](https://en.cppreference.com/w/cpp/atomic/memory_order)
- Preshing, *Acquire and Release Semantics* (2012)
- Boehm & Adve, *Foundations of the C++ Concurrency Memory Model* (2008)
- [Robin Hood Hashing](https://github.com/martinus/robin-hood-hashing)
- [What Every Programmer Should Know About Memory](https://people.freebsd.org/~lstewart/articles/cpumemory.pdf) (Drepper, 2007)

### Comparison Tools
- BioPython: https://biopython.org/
- PyFastX: https://github.com/lmdu/pyfastx
- SeqKit: https://bioinf.shenwei.me/seqkit/

---

**Status:** ✅ Production Release  
**Next Version:** v1.3.0 "Hrunting" (Q3 2026)

*"Trace On" - Projecting legendary performance from genomic data across eons.* ⚔️