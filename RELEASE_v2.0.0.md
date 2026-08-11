# TracEon v2.0.0 — "Durandal"

**Release date:** 2026-08-11 (target)
**Codename:** Durandal — per the naming convention (alphabetical order of Emiya Shirou's Noble Phantasms: Caladbolg → Caliburn → Durandal → ...)

**Tag:** `v2.0.0`

---

## Summary

TracEon v2.0.0 is a **major release** focused on data integrity, a new binary format, and API safety. All prior hardening work (parser/gzip/binary fixes, string_view keys, wyhash+reserve performance) plus the new v4 format, k-mer C API, and lifecycle contract land together. **No backward compatibility is provided** — the library has no external adopters, so a clean break is intentional (SemVer major).

## Breaking changes (migration required)

1. **Binary cache format v4 only**: `.traceon` files now use the v4 format (`"TRO\x04"`, LZ4 Frame + whole-payload CRC32C, streamed). v1/v2/v3 caches are **rejected** with a clear error — regenerate caches with `saveBinary()`.
2. **GZIP integrity enforcement**: corrupt/truncated `.gz` files now **throw** `std::runtime_error` instead of loading silently.
3. **Malformed binary caches rejected**: crafted/truncated caches fail loudly and never publish partial state (failure atomicity).
4. **`set()`/`addEntry()` after a load throws** `std::logic_error` (immutable-after-load).
5. **`std::string_view` keys**: index keys are now string_view-backed — recompile consumers; view lifetime is tied to the loaded snapshot (reader-quiescence contract, ADR-001).
6. **No ABI stability promised** across 1.x → 2.x.

## New in v2.0.0

- **v4 binary format + CRC32C**: whole-logical-payload checksum streamed through LZ4F (no second pass), SSE4.2/AArch64/table dispatch, exact-length + exact-frame-termination validation.
- **k-mer C API** (opt-in `traceon_kmer` target): hardened C boundary — noexcept + status codes + `kmerindex_last_error()`, caller-owned iterators, freeze semantics, validated TRKI loading (overflow-safe, bounds-checked).
- **Reader-quiescence lifecycle contract**: reads safe only while the same snapshot is installed; returned views invalid after `clearCache()`/`reload`/destruction. Debug diagnostics via `TRACEON_DEBUG_LIFECYCLE`.
- **Performance**: string_view keys + wyhash + pre-reserved maps — +33% isolated insert throughput, −14..−28% parse-dominated loads.

## Verification (all run on this release)

- **Tests**: 128 cases / 4,195 assertions pass (normal build); suite split into 6 domain executables + kmer_tests (15/490) — all CTest-registered.
- **Sanitizers**: ASan+UBSan — full suite passes (127/4,193 with the single `[oom]` test excluded — documented RLIMIT_AS/shadow incompatibility, passes in normal build); kmer under ASan passes; 80k/80k repro passes under ASan.
- **Repro**: 80k-record FASTQ regression (`@`-in-quality) — 80k/80k both control and `+`-quality files.
- **Benchmarks** (clean machine, median of 3, 1000MB): WGS_GZ **3.06s** · RefGenome_GZ **1.37s** · WGS plain **1.00s**. CRC32C overhead negligible (v4 ≈ pre-v4 within noise). Absolute numbers are host-dependent; relative A/B on identical fixtures is the regression gate.

## Repository

- `outputs/traceon-v2-design-review.md` — Feynman design review (checksum choice, thread-safety, perf freeze, versioning, git flow)
- `docs/architecture/ADR-005-traceon-v4-binary-format.md` — v4 wire spec
- `docs/architecture/ADR-001-lock-free-reads.md` — lifecycle contract
- CHANGELOG.md — full history + migration section
