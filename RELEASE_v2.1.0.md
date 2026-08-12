# TracEon v2.1.0 "Excalibur" — Release Notes

**Date:** 2026-08-12
**Release type:** Minor — verification + tooling hardening (no format change; v4 `.traceon` caches remain compatible)

## What's in this release

### 1. Fuzzing gate (libFuzzer) — the last never-run verification layer
- 7 fuzz targets covering the full attack surface: v4 binary loader, gzip loader, FASTQ parser, FASTA parser, kmer encoder, kmer C API, TRKI mmap loader
- 83-file curated seed corpus (including byte-exact vuln PoCs for vuln-0001/0003/0006)
- CI workflow (`fuzz.yml`): runs on every push to `main` + PRs + weekly; 120s budget per target; crash reproducers uploaded as artifacts even on failure
- OSS-Fuzz canonical pattern (`-fsanitize=fuzzer-no-link` on libs, single runtime on targets) — verified green on clang 22 (Fedora) and CI (ubuntu clang)
- All 7 targets validated with real libFuzzer runs: zero findings

### 2. kmer C API hardening (thread-safety + total no-throw)
- `frozen` flag is now `std::atomic<bool>` (release/acquire) — concurrent `kmerindex_get()` from mapping worker threads no longer races (verified with an 8-thread TSan harness + positive control)
- Error path is now allocation-free: fixed-size thread-local buffer replaces `std::string` — the advertised no-throw C boundary is now total, even under OOM
- `assert(inserted)` audit — callers check `== 1` explicitly (the `-1` exception status was a silent trap)
- New test: 100k keys → freeze → 8 threads × 25k gets + 50k insert attempts (all must fail frozen). kmer_tests now 150,532 assertions in 16 cases

### 3. Real-world aligner validation (the kmer C API's raison d'être)
- **`woosflex/minimap2` `traceon-backend`** — minimap2 v2.31 fork with:
  - `TRACEON=1`: per-bucket `traceon_kmer` table (ankerl::unordered_dense) replacing `khash_t(idx)`
  - `tcache v2`: mmap'd open-addressing prebuilt index — zero-rebuild load, O(1) lookup, whole-file CRC32C
  - Recurring load+map: **−13.6% wall at full chr1 (249Mb)**, PAF byte-identical; the delta is pure load phase (mmap+CRC ~134ms vs khash rebuild)
- **`woosflex/Winnowmap` `traceon-backend`** — same integration on the second mapper:
  - Recurring load+map: **−80% wall / −31% RSS at chr1-50Mb** (stock Winnowmap has no working prebuilt-index path — `-d` disabled upstream + an `applyWeight` bloom-deref crash; every stock run rebuilds)
  - Found + documented two upstream bugs (applyWeight crash; uninitialized `rep_len` in the SV-aware MCAS path) — not fixed in the stock path, byte-identity preserved
- Every integration verified: **PAF byte-identical to stock** at all scales, `.mmi` round-trips both directions, corruption detected via CRC

## Verification summary
- Unit suite: **128 test cases / 4,195 assertions** (plus kmer 16 / 150,532)
- Sanitizers: ASan + UBSan (split suite, `~[oom]`), TSan (C API thread test), libFuzzer (7 targets)
- Release builds: CI on 6 platforms (ubuntu/windows/macOS × x86_64/arm64)
- Performance check: no regression vs v2.0.0 baselines

## Migration notes
- None — v2.1.0 is additive. `.traceon` v4 caches, the C++ API, and the kmer C API are unchanged.

## Codename
"Excalibur" — 4th in the alphabetical Noble Phantasm sequence (Caladbolg → Caliburn → Durandal → **Excalibur** → Gáe Bolg → …)
