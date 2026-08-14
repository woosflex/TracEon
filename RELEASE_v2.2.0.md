# TracEon v2.2.0 "Gáe Bolg" — Release Notes

**Date:** 2026-08-14
**Release type:** Minor — aligner validation matrix completion (no format/API change; v4 `.traceon` caches remain compatible)

## What's in this release

### The full aligner validation matrix — 5 real tools, all PAF byte-identical

The kmer C API's real-world proof is now complete across the exact-semantic-fit candidates identified by the Feynman survey (`outputs/aligner-khash-traceon-survey.md`). Every fork keeps the stock build byte-identical and adds an opt-in `TRACEON=1` mode (per-bucket `traceon_kmer` table) + `tcache` (mmap'd open-addressing prebuilt index — zero-rebuild load, whole-file CRC32C).

| Tool | Fork branch | Recurring load+map result |
|---|---|---|
| **minimap2** (v2.31) | `woosflex/minimap2` `traceon-backend` | **−13.6%** wall at full chr1 (249Mb); RSS ≈ stock |
| **Winnowmap** (72403b3) | `woosflex/Winnowmap` `traceon-backend` | **−80%** wall / −31% RSS (chr1-50Mb); stock has NO prebuilt-index path |
| **mm2-fast** (14fe36c) | `woosflex/mm2-fast` `traceon-backend` | **−39%** wall / −33% RSS (chr1-50Mb, Hermes-verified) |
| **minigraph** (2f569eb) | `woosflex/minigraph` `traceon-backend` | First **khashl** (linear-probing) adapter; GAF/PAF identity; fixtures too small to time |
| **BLEND** (6f19e37) | `woosflex/BLEND` `traceon-backend` | PAF identity across 7 configs; load **~3.3×** faster (tcache vs khash rebuild) |

Key findings across the matrix:
1. **The table swap alone is memory-positive, time-neutral** — mapping wall time is dominated by chaining/alignment, not minimizer lookups.
2. **The recurring-run win is the zero-rebuild load** — mmap + CRC in ~ms vs khash rebuild in hundreds of ms–minutes (grows with reference size; whole-human `.mmi` rebuild is minutes).
3. **Winnowmap's stock has no working prebuilt-index path** (upstream `-d` disabled + `applyWeight` crash) — every stock run rebuilds; tcache fixes it.
4. **LISA_HASH** (mm2-fast's learned-index mode) kept separate — a stock-khash comparator, not part of the swap.

### Upstream bugs found + documented (none fixed in stock paths — byte-identity preserved)
- **Winnowmap**: `applyWeight` bloom-filter deref crashes any prebuilt-index load; uninitialized `rep_len` read in the SV-aware MCAS path (`rl:`/`mapq` binary-dependent garbage)
- **minigraph**: pristine r613 crashes on its own `test/MT.gfa` at `-O3 -msse4` under GCC 16 (unaligned u64 casts vectorized into faulting SSE loads) — fixed in the port with `memcpy` loads, output byte-identical to `-O0`; also found + fixed a `gi->g_own` leak (ASan/LSan clean)
- **mm2-fast**: safestringlib submodule needs a local GCC-16 CFLAGS patch (documented)

## Verification summary
- Unit suite: **128 test cases / 4,195 assertions** (unchanged — no core code changes in v2.2.0)
- kmer: 16 / 150,532 assertions
- Sanitizers: ASan + UBSan + TSan + libFuzzer (all green from v2.1.0, unchanged)
- Every fork gate: stock machine-code byte-identical · PAF/GAF byte-identity · `.mmi`/tcache round-trips · corruption detection (CRC, no crash)

## Migration notes
- None — v2.2.0 is additive. v4 caches, C++ API, and kmer C API unchanged.

## Codename
"Gáe Bolg" — 5th in the alphabetical Noble Phantasm sequence (Caladbolg → Caliburn → Durandal → Excalibur → **Gáe Bolg** → Harpe → Hrunting)
