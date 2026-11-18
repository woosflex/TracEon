# Internal minimap2 lookup pilot — 2026-09-22

Follow-up: the [resident multi-process experiment](shared-minimap-index-pilot.md)
tests cache-once index reuse with a native-minimap2 residency control.

TracEon's existing `kmerindex_c_api` successfully replaced the lookup performed
by minimap2's `mm_idx_get`. All 75 mapping runs produced matching canonical SAM
results. This demonstrates internal lookup compatibility, not shared indexing
across processes or a general performance improvement.

The experiment used minimap2 commit
`79c9cc186b95f50bd899f69b48eba995ced810c6`, compiled twice with the same Makefile
defaults. The modified copy links the portable Release TracEon kmer archive.
Every seed lookup uses TracEon's C API, with no khash lookup fallback. Singleton
flags and posting-list offsets retain minimap2's semantics; posting lists,
sketching, chaining and alignment remain minimap2's implementation.

For this compatibility pilot, original khash tables remain available for index
statistics and serialization. Each bucket is copied into TracEon after building
or loading, and every copied key/value is checked before mapping. The extra
table and exhaustive verification add startup time and memory. The measured
whole-process times therefore cannot isolate hash lookup performance or predict
the performance of a finished shared index.

## Results

Six repetitions per arm, alternating execution order, at one and four mapper
threads. Local, nonexclusive host; no CPU affinity or explicit cold-cache
control. Times include loading, conversion/verification where applicable,
mapping and SAM output. No statistical speedup claim is made.

| Input | Threads | Original median seconds | TracEon median seconds |
|---|---:|---:|---:|
| E. coli, 10,000 simulated single-end short reads | 1 | 0.1208 | 0.1635 |
| Same | 4 | 0.0705 | 0.1166 |
| E. coli, first 1,000 corrected PacBio reads | 1 | 2.0727 | 2.0460 |
| Same | 4 | 0.7232 | 0.7563 |
| SIRV, 53 Iso-Seq reads, spliced alignment | 1 | 0.1921 | 0.2117 |
| Same | 4 | 0.0914 | 0.0956 |

The 72 indexed mapping runs plus three direct-reference construction/mapping
runs matched all SAM alignment fields and tag values, ignoring only record/tag
ordering and `@PG` provenance. Primary mapped counts were 9,988/10,000,
982/1,000 and 53/53 respectively. The three indexes written by the patched
mapper were byte-identical to the original mapper's indexes. This establishes
equivalence on these inputs, not independent biological correctness.

The existing C ABI check passed, and the kmer suite passed 153,638 assertions
across 12 cases. The pilot itself retains executable assertions for copied
entries, SAM equivalence and serialized-index equality.

## Reproduction and evidence

```sh
python benchmarks/minimap_lookup_pilot.py \
  --source build/research-inputs/minimap2-79c9cc186b95f50bd899f69b48eba995ced810c6 \
  --library /tmp/traceon-release-check-20260922/libtraceon_kmer.a \
  --cookbook build/research-inputs/cookbook \
  --output build/minimap-lookup-pilot-new
```

The output directory must be new. Local raw SAM, indexes and binaries are in
`build/minimap-lookup-pilot-2026-09-22/`. The durable
[evidence archive](results/minimap-lookup-pilot-2026-09-22.tar.gz) contains the
exact source patch, build logs, commands, timing/RSS observations, SAM hashes,
input/binary hashes, provenance and harness snapshot. Large SAM files and
binaries remain local and are regenerable.

The next architecture experiment is an immutable minimizer/postings index that
independent mapping processes can attach to without rebuilding private tables.
That needs its own shared-memory or batched-service integration and aggregate
memory/throughput tests. The current resident sequence-record service does not
implement that interface. This result supports pursuing that experiment; it
does not establish the cache-once/many-processes benefit inside minimap2 yet.
