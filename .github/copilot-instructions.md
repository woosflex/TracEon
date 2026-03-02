# Copilot Instructions for TracEon

## What This Project Is

TracEon is a zero-copy, lock-free genomic data caching library for FASTA/FASTQ files (including `.gz`). It indexes sequences into a single arena buffer and returns `std::string_view` handles — never copying data on reads. Requires C++20.

## Build Commands

```bash
# Configure and build (Release is required for performance-sensitive work)
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Debug build (for development/testing)
cmake -B build -S . -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
```

## Test Commands

```bash
# Run all tests
cd build && ctest --output-on-failure
# or directly
./build/unit_tests

# Run a single test suite (Catch2 tag syntax)
./build/unit_tests "[cache]"
./build/unit_tests "[gzip]"
./build/unit_tests "[strategy]"

# Verbose output
./build/unit_tests -d yes
```

Tests use **Catch2 v3** (auto-fetched via CMake `FetchContent` — no manual install needed).

## Architecture

The public API lives in `include/Cache.h` (with `include/TracEon.h` as a convenience include). Callers interact only with `Cache`.

Internally, `Cache` routes to one of two backends:
- **`SmartStrategy`** — used for file-backed loading. Owns a single `text_arena_` (arena allocator) and a `robin_hood::flat_map` (`file_cache_`) mapping sequence IDs to `std::string_view` slices into the arena. Handles GZIP streaming via zlib (1 MB chunks).
- **`robin_hood::unordered_map`** — used for manually inserted key/value pairs.

Lock-free read safety is enforced via `std::atomic<bool> data_ready_` with acquire/release semantics. Reads that arrive before loading completes spin-wait on this flag.

**Index modes** are selected automatically based on input:
- `IndexMode::GENOME` → `GenomeIndex` (string keys, e.g. chromosome names)
- `IndexMode::NGS` → `NGSIndex` (uint64_t hash keys, for short-read datasets)

Format detection happens by file extension first, then magic bytes (GZIP: `\x1f\x8b`; FASTQ: leading `@`).

## Key Conventions

**Naming:**
- Classes: `PascalCase` (`Cache`, `SmartStrategy`, `IEncodingStrategy`)
- Methods: `camelCase` (`loadFile()`, `getView()`, `hasSequence()`)
- Private/member variables: `snake_case_` with trailing underscore (`text_arena_`, `data_ready_`)
- Enums: `PascalCase::UPPER_SNAKE` (`FileFormat::DNA_FASTA`, `IndexMode::GENOME`)

**Return types on hot paths always use `std::string_view`** — never `std::string`. Nullable returns use `std::optional`.

**Third-party libraries** are in `third_party/` (vendored: `robin_hood.h`, zlib, lz4). Do not replace these with system libraries.

**`lz4`** is integrated but not yet active — it is reserved for future binary cache serialization.

**Architecture Decision Records** are in `docs/` — check there before making structural changes to understand prior design decisions.

## CI Workflows

Three workflows run on PRs (`.github/workflows/`):
- `cmake-multi-platform.yml` — builds on Linux, macOS, Windows
- `unit-tests.yml` — runs the Catch2 test suite
- `performance-check.yml` — runs benchmark comparisons
