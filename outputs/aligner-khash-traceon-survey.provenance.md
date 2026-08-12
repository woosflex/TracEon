# Provenance: aligner khash/TracEon survey

## Canonical artifact

`outputs/aligner-khash-traceon-survey.md`

## Scope

Surveyed the explicitly requested tools: Winnowmap, minigraph, minigraph-cactus, strobemap/StrobeMap, strobealign, BLEND, HiFiMapper, LRA, meryl, mm2-fast, and minimap2 forks/stock modes. Added MoleMap and ParBLiSS/kmerhash as non-target structural controls.

## Evidence collection

- Current official repository snapshots were shallow-cloned under `/tmp/aligner-research` and inspected with `grep`, `sed`, `nl`, and `git rev-parse`.
- Official GitHub blob/raw pages were fetched for exact source claims, especially the `idx` declarations, key/value types, lookup/build paths, and Cactus external `minigraph` boundary.
- Web search was used for current repository/paper discovery and primary paper identification.
- Primary papers/docs were consulted for algorithm context: Winnowmap2, Minigraph-Cactus, BLEND, Strobemers, Strobealign, LRA, and the minimap2 cookbook/paper.
- Local current TracEon files were inspected directly. The report uses current checkout commit `9c54742a7699e0701c0d836ba68c812ed510a6fe`, not stale observations from earlier minimap2 integration notes.

## Key verified source observations

- Winnowmap, mm2-fast stock, BLEND, and canonical minimap2 instantiate a classic `uint64_t -> uint64_t` khash index with `hash(a)=a>>1`, equality on `a>>1`, bit-0 singleton flag, and packed multi-hit offset/count.
- Minigraph core instantiates the same logical type/semantics with `khashl.h`; the khashl probing implementation is linear open addressing.
- Minigraph-Cactus invokes an external `minigraph` executable; Cactus does not contain the minigraph table implementation in the inspected tree.
- StrobeMap uses robin-hood open addressing with a `uint64_t -> tuple<uint64_t,uint32_t>` map, not the TracEon packed value contract.
- Current strobealign production lookup uses sorted vectors and bucket ranges, not a hash table.
- HiFiMapper uses `std::unordered_map<uint64_t,uint64_t>` with `>>1` hash/equality and a low-bit convention, but the source does not establish open addressing.
- LRA's primary index is sorted vectors; its visible `unordered_map<Tuple,int>` is statistics-only.
- meryl uses packed prefix/suffix arrays and is not an aligner or open-addressing replacement target.

## Local commands/observations

1. `git clone --depth 1` official repositories for the named projects; source grep identified exact declarations and index paths.
2. `git rev-parse HEAD` recorded pinned snapshot commits and dates used in the report.
3. Current TracEon source inspection confirmed `std::atomic<bool> frozen`, fixed-size thread-local diagnostics, custom `ankerl::unordered_dense` hash/equality, and opt-in `traceon_kmer` target.
4. No aligner fork was changed.
5. No end-to-end A/B benchmark, index build, mapping run, performance counter run, or output-equivalence test was executed for this survey.

## Verification status

- Independent reviewer and verifier passes completed. The verifier found the local TracEon commit was one commit ahead of the public remote and therefore not reachable as a GitHub blob; the report now cites the public repository plus local file paths/commit and explicitly labels the local evidence. MoleMap references were changed from dead pinned URLs to live moving-branch comparator links with a pinning caveat.
- Direct source observations are cited with pinned URLs in the canonical report.
- Quantitative impact is intentionally marked as hypothesis/measurement plan only.
- Maintenance claims are snapshot-level observations, not support guarantees.
- AI-policy claims are bounded negative findings from visible repository docs; absence of a policy is not permission or prohibition evidence.
- StrobeMap root licensing remains unresolved because the inspected source tree lacked a root license while Bioconda metadata labels the package GPL-3.0.

## Known residual risks

- A real adapter still requires tool-specific build/link changes, caller-side status handling, explicit freeze before mapping workers, and logical output tests.
- The C API adds an out-of-line boundary and may not beat khash's cheap inline shifted hash.
- Per-bucket map overhead, occurrence distribution, cache state, compiler/LTO, allocator, and ISA can dominate results.
- `std::unordered_map` storage in HiFiMapper is implementation-dependent.
- No claim is made that arbitrary minimap2 forks share the inspected table path.
