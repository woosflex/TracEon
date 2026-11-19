# TracEon v3.0.0-dev — Pre-Release Notes

**Date:** 2026-11-19  
**Release type:** Pre-release / Development Candidate (v3.0.0-dev)

## Highlights

### 1. Memory-Mapped Cache Storage (`CacheStorage::Mapped`)
- Portable v5 binary cache format introducing codec 0 (uncompressed, memory-mapped payloads).
- Independent OS reader processes can share payload pages directly through the OS page cache without duplicating memory in private heaps or running a background daemon.
- Zero-copy sequence views (`getView()`) borrow directly from mapped memory.
- Backward compatibility: cache restoration safely reads both v4 and v5 cache formats.

### 2. Atomic POSIX Persistence Guarantees
- Cache saves use synchronized temporary files and atomic POSIX replacement (`rename`).
- Parent-directory synchronization (`fsync`) ensures directory metadata durability across system crashes or power interruptions.
- Pre-commit failures preserve the original cache file without corruption. Symlink and device destinations are explicitly rejected.

### 3. Remote Protocol v2 & Batch Retrieval
- **Snapshot Identity**: Canonical SHA-256 snapshot verification pins content between client and server, preventing version drift.
- **Batch Retrieval**: Added opt-in `getBatch()` (and CLI `--batch-size`) bounded to 256 IDs, preserving order, duplicates, misses, and operation deadlines.
- **Full FASTQ Support**: Wire protocol natively transmits both sequence and quality strings with Castagnoli CRC32C verification.
- **Deadline-Bounded I/O**: Per-frame absolute deadlines cover both read and write operations; client admission cap (default 64) prevents connection exhaustion.
- **Cross-Platform Socket Safety**: Enforced `O_NONBLOCK` descriptor modes, `SO_NOSIGPIPE`, and global `SIGPIPE` suppression for clean execution on Linux and macOS under ThreadSanitizer.

### 4. Unified CLI Tool (`traceon`)
- Added unified multi-command CLI:
  - `traceon build`: Build portable binary caches from FASTA/FASTQ inputs with optional `--mapped` mode.
  - `traceon serve`: Launch resident zero-copy TCP lookup service with configured admission caps.
  - `traceon query`: Perform sub-millisecond single and batch queries from stdin or command-line arguments against a resident server.
  - `traceon info`: Inspect served snapshot ID, capabilities, and record counts.

### 5. Relocatable Modern CMake Export
- Full `find_package(TracEon 3 REQUIRED COMPONENTS core)` support with relocatable targets (`TracEon::core`, `TracEon::remote`, `TracEon::kmer`).
- Decoupled testing and dependency options (`BUILD_TESTING=OFF`, `TRACEON_BUILD_SERVER=OFF`) allowing lightweight embedded or core-only distributions.

## Verification

- **Comprehensive Test Suite**: 23/23 CTest targets passing across Linux and macOS.
- **ThreadSanitizer Clean**: Verified under Apple Clang and GCC ThreadSanitizer (`-fsanitize=thread`).
- **Installed Consumer Validation**: Automated downstream CMake build and consumption check verified in CI.
- **Fuzz Hardened**: 7 libFuzzer targets covering parser, cache, k-mer, and protocol decoding boundaries.

See [CHANGELOG.md](CHANGELOG.md) for the complete change ledger.
