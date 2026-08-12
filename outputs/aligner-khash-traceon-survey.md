# DNA/RNA aligners using khash-style minimizer/k-mer tables: TracEon validation survey

**Research status:** source-backed survey and benchmark design. No aligner fork was modified and no end-to-end A/B performance result was generated. All performance statements below are hypotheses or measurements to collect, never claimed results.

**Question:** beyond minimap2, which DNA/RNA aligners use `khash`/`klib` or structurally similar open-addressing minimizer/k-mer indexes, and where could TracEon's `kmerindex` C API replace the table?

## Executive recommendation

Among the explicitly named and inspected projects, the strongest first-wave validation targets are:

1. **Winnowmap** — exact minimap2-derived `khash_t(idx)` semantics, high value on repetitive long-read mapping, and a relatively recent source snapshot.
2. **mm2-fast, stock seeding mode** — exact `uint64_t -> uint64_t` minimap2 index, but benchmark its optional `LISA_HASH` mode separately because it can bypass the khash lookup path.
3. **minigraph** — exact key/value and bit-flag semantics implemented with lh3's newer `khashl` linear-probing table; this validates TracEon beyond a literal minimap2 fork and adds graph-mapping relevance.
4. **BLEND** — minimap2-derived khash index with the same table contract, but a modified seed/hash-generation algorithm; useful for testing whether the adapter remains correct when seed generation changes.
5. **HiFiMapper, optional** — C++ `std::unordered_map<uint64_t,uint64_t>` with the same `>>1` hash/equality and a low-bit singleton convention. It is a semantic control, not proven open addressing; its inspected snapshot is older and its visible repository scale is smaller than the main candidates. This is a heuristic prioritization, not a field-wide maintenance/popularity ranking.

**Do not present StrobeMap, strobealign, LRA, or meryl as direct drop-in candidates.** Their production indexes are, respectively, a robin-hood map with a tuple value, a sorted vector/range index, a sorted vector/range index, and packed prefix/suffix arrays. They are useful negative or architecture controls.

The common semantic contract is stricter than “a 64-bit hash map”:

```c
#define idx_hash(a) ((a) >> 1)
#define idx_eq(a,b) (((a) >> 1) == ((b) >> 1))
// key bit 0: singleton marker
// singleton: stored key = base | 1, value = one reference position
// repeated:  stored key = base, value = (position-array-offset << 32) | count
```

TracEon's current `traceon_kmer` implementation mirrors this with an `ankerl::unordered_dense::map<uint64_t,uint64_t>` and custom hash/equality, a C ABI, explicit freeze, atomic phase state, and a fixed non-allocating diagnostic buffer. Those properties make it a plausible **source-level experimental backend**, not an ABI-compatible replacement for khash's raw layout or iteration order.

---

## 1. Scope and method

I inspected pinned official repository snapshots and relevant primary papers/docs. For each project I checked:

- whether the production seeding/index path uses classic `khash.h`, `khashl.h`, or another open-addressing table;
- exact key/value types and hash/equality semantics;
- whether bit 0 is a flag and whether values are packed;
- implementation language and visible license;
- snapshot-level maintenance evidence;
- visible `CONTRIBUTING`, `README`, `SECURITY`, or code-of-conduct AI-policy text.

“Open addressing” means the production table probes an in-memory slot array. A `std::unordered_map` does **not** count as proven open addressing because the C++ standard does not specify its storage strategy. An auxiliary statistics map does not count as the aligner's seeding index.

The repository snapshots were pinned for source inspection. Activity labels below mean “activity visible at the inspected snapshot,” not a guarantee of future maintenance.

---

## 2. TracEon compatibility target

### 2.1 Current TracEon API

The current local checkout at commit `9c54742a7699e0701c0d836ba68c812ed510a6fe` (one commit ahead of the public `origin/main` at inspection) exposes an opaque handle and C functions for; the local evidence is `include/kmerindex_c_api.h` and `src/kmerindex_c_api.cpp`:

- `uint64_t` key/value insertion and lookup;
- reserve, size, freeze, and destruction;
- caller-owned iteration;
- status-code error handling;
- stored-key return so callers can inspect the singleton bit.

The implementation uses:

```cpp
using MinimizerMap = ankerl::unordered_dense::map<
    uint64_t, uint64_t, MinimizerKeyHash, MinimizerKeyEq>;
```

where the custom hash and equality operate on `key >> 1`. The current implementation uses `std::atomic<bool>` for the freeze transition and a fixed-size thread-local error buffer. These are important because mapping workers perform concurrent lookups and `get()` returns an interior value pointer.

Local source: `include/kmerindex_c_api.h`, `src/kmerindex_c_api.cpp`, and `CMakeLists.txt` at the commit above. Public repository context: [TracEon](https://github.com/woosflex/TracEon). Dependency: [ankerl::unordered_dense](https://github.com/martinus/unordered_dense).

### 2.2 What “drop-in” can mean here

There are three distinct claims; only the first is plausible for the exact candidates:

1. **Logical table replacement:** same key matching, singleton/multi-hit interpretation, positions, and mapping output. This is the immediate goal.
2. **Source/API replacement:** existing C code calls a narrow adapter rather than `kh_get`/`kh_put`. This is feasible with a C shim and a C++ final link.
3. **Physical replacement:** same slot order, `khint_t` behavior, pointer/iterator invalidation, serialized byte order, or raw khash ABI. This is not a goal and should not be claimed.

The adapter must preserve the caller's grouping invariant: minimizers are sorted/grouped before insertion, so every logical base key should be inserted once. It must check every create/reserve/insert status and freeze every index bucket before mapping workers begin. The current `traceon_kmer` API exposes the needed primitives, but the atomic flag is not a substitute for a build-to-read barrier: callers must guarantee that no insert/reserve operation can overlap any worker lookup, iteration, or returned-value pointer. Each integration still needs caller-side correctness and lifecycle tests.

---

# Q1 — Survey of aligners and related tools

## 3. Exact khash-lineage candidates

| Tool | Repository / inspected snapshot | Production table and exact type | Hash/equality and flag | Language / license | Maintenance evidence | Candidate status |
|---|---|---|---|---|---|---|
| **Winnowmap** | [marbl/Winnowmap](https://github.com/marbl/Winnowmap), `72403b3` | Vendors `src/khash.h`; `KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)` in [`src/index.c`](https://github.com/marbl/Winnowmap/blob/72403b350fafd3aec98d1a6c67123d0d0d1dff7c/src/index.c#L13-L38) | `idx_hash(a)=a>>1`; `idx_eq(a,b)=((a>>1)==(b>>1))`; key bit 0 marks singleton; repeated values pack offset/count into `uint64_t` | C/C++ build (`.c` compiled through C++ compiler); joint-work license with minimap2 MIT/public-domain/source-specific notices | Pinned snapshot commit dated 2025-12-16; recent source activity observed | **Yes — exact semantic fit** |
| **minigraph** | [lh3/minigraph](https://github.com/lh3/minigraph), `2f569eb` | Vendors `khashl.h`; `KHASHL_MAP_INIT(KH_LOCAL, idxhash_t, mg_hidx, uint64_t, uint64_t, idx_hash, idx_eq)` in [`index.c`](https://github.com/lh3/minigraph/blob/2f569ebe3071fcb242f1bf0eabbb917941b47239/index.c#L1-L17) | Same `>>1` hash/equality and singleton/multi-hit encoding; `khashl` uses linear probing | C/C99, MIT | Pinned snapshot dated 2025-08-11; NEWS has a 0.21-r606 release entry | **Yes — requires khashl adapter mapping** |
| **mm2-fast** | [bwa-mem2/mm2-fast](https://github.com/bwa-mem2/mm2-fast), `14fe36c` | Vendors `khash.h`; stock [`index.c`](https://github.com/bwa-mem2/mm2-fast/blob/14fe36c100f6c2aab224d000f3903ca5909640cd/index.c#L17-L33) has the exact minimap2 `uint64_t -> uint64_t` map | Same shifted hash/equality, low-bit singleton marker, packed repeated payload | C/C++; MIT text for minimap2-derived code | Pinned snapshot dated 2024-06-04; less recent than Winnowmap/minigraph | **Yes in stock mode; separate LISA mode** |
| **BLEND** | [CMU-SAFARI/BLEND](https://github.com/CMU-SAFARI/BLEND), `6f19e37` | Vendors `khash.h`; exact [`KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)`](https://github.com/CMU-SAFARI/BLEND/blob/6f19e37f5120fc8d3b389978e7e4056c51810d0b/src/index.c#L16-L33) | Same `>>1` semantics and singleton bit; BLEND changes seed/hash generation and some shifts, so source-local tests are required | C, MIT | Pinned snapshot dated 2023-05-10; research code with older visible activity | **Yes — exact table shape, specialized algorithm** |
| **Canonical minimap2** | [lh3/minimap2](https://github.com/lh3/minimap2), v2.31/r1302 source | Baseline classic khash `uint64_t -> uint64_t`, plus a separate string-name map | Exact source contract in [`index.c`](https://github.com/lh3/minimap2/blob/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/index.c#L11-L33) and construction/lookup at [L93-L110](https://github.com/lh3/minimap2/blob/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/index.c#L93-L110) | C-oriented, MIT | Latest baseline for the primary validation | **Target #1, not “beyond” list** |

### Winnowmap details

Winnowmap is a minimap2-derived long-read mapper with weighted minimizer sampling and additional repetitive-region handling. Its `src/index.c` retains the exact `idx` khash shape. Lookup selects the bucket from low minimizer bits, forms the bit-0-clear query key, then uses the stored key's bit to choose singleton versus multi-hit interpretation. The separate `str` map is unrelated to seeding.

The critical source evidence is the [`idx` declaration](https://github.com/marbl/Winnowmap/blob/72403b350fafd3aec98d1a6c67123d0d0d1dff7c/src/index.c#L24-L38) and the [`mm_idx_get`/construction path](https://github.com/marbl/Winnowmap/blob/72403b350fafd3aec98d1a6c67123d0d0d1dff7c/src/index.c#L88-L105). The vendored header identifies klib khash 0.2.8 and uses open addressing with power-of-two capacity and probing: [`khash.h`](https://github.com/marbl/Winnowmap/blob/72403b350fafd3aec98d1a6c67123d0d0d1dff7c/src/khash.h#L117-L126).

The algorithmic paper is Jain et al., **“Long-read mapping to repetitive reference sequences using Winnowmap2,” Nature Methods (2022)**, DOI [10.1038/s41592-022-01457-8](https://doi.org/10.1038/s41592-022-01457-8). The table swap would not test Winnowmap's weighted-minimizer algorithm itself; it would test the storage backend under that workload.

### minigraph details

Minigraph uses `khashl.h`, not classic `khash.h`. Its core minimizer index is still the exact `uint64_t -> uint64_t` contract. `khashl`'s implementation visibly advances by one slot under collision, establishing open addressing with linear probing: [`khashl.h`](https://github.com/lh3/minigraph/blob/2f569ebe3071fcb242f1bf0eabbb917941b47239/khashl.h#L189-L213). The core graph-edit code contains other `khashl` maps with different types, including `uint32_t -> int32_t` and string maps; those are not direct TracEon targets.

The minigraph core index declaration and lookup are at [`index.c`](https://github.com/lh3/minigraph/blob/2f569ebe3071fcb242f1bf0eabbb917941b47239/index.c#L1-L17) and [lookup/build](https://github.com/lh3/minigraph/blob/2f569ebe3071fcb242f1bf0eabbb917941b47239/index.c#L50-L72). The project is MIT-licensed and C-oriented: [Makefile/license](https://github.com/lh3/minigraph/blob/2f569ebe3071fcb242f1bf0eabbb917941b47239/Makefile#L1-L10).

**minigraph-cactus is not a separate table implementation.** The Cactus pipeline constructs and executes an external `minigraph` command in [`cactus_minigraph.py`](https://github.com/ComparativeGenomicsToolkit/cactus/blob/68252c2a9406bbf82b62434e7060403d44bed5b8/src/cactus/refmap/cactus_minigraph.py#L506-L546). Cactus itself does not vendor minigraph's `index.c`/`khashl.h` in the inspected tree. A TracEon experiment therefore targets the separately rebuilt minigraph binary, then plugs that binary into Cactus; editing the Python wrapper alone cannot replace the internal table.

The relevant primary paper is **“Pangenome graph construction from genome alignments with Minigraph-Cactus,” Nature Biotechnology (2023)**, DOI [10.1038/s41587-023-01793-w](https://doi.org/10.1038/s41587-023-01793-w).

### mm2-fast details

mm2-fast retains the exact stock minimap2 table in [`index.c`](https://github.com/bwa-mem2/mm2-fast/blob/14fe36c100f6c2aab224d000f3903ca5909640cd/index.c#L67-L72), with the same lookup and packed payload at [L93-L110](https://github.com/bwa-mem2/mm2-fast/blob/14fe36c100f6c2aab224d000f3903ca5909640cd/index.c#L93-L110).

However, mm2-fast documents an optional `LISA_HASH` learned-index path. In [`seed.c`](https://github.com/bwa-mem2/mm2-fast/blob/14fe36c100f6c2aab224d000f3903ca5909640cd/seed.c#L69-L111), stock builds call the normal lookup while `LISA_HASH` builds call an external learned-index batch lookup. The README says learned-index seeding is disabled by default: [build modes](https://github.com/bwa-mem2/mm2-fast/blob/14fe36c100f6c2aab224d000f3903ca5909640cd/README.md#L46-L61).

Therefore report two separate experiments:

- `mm2-fast-stock-khash` vs `mm2-fast-traceon`;
- `mm2-fast-LISA` as a separate comparator, not as evidence about the khash swap.

The repository is MIT-licensed for the minimap2-derived code and its README describes AVX2/AVX-512 optimizations: [README](https://github.com/bwa-mem2/mm2-fast/blob/14fe36c100f6c2aab224d000f3903ca5909640cd/README.md).

### BLEND details

BLEND integrates fuzzy seed generation into a minimap2-derived codebase. Its storage table is still the classic minimap2-style khash map, while the seed values and sketching path differ. The exact table declaration is at [`src/index.c`](https://github.com/CMU-SAFARI/BLEND/blob/6f19e37f5120fc8d3b389978e7e4056c51810d0b/src/index.c#L16-L33); lookup and build are at [L83-L110](https://github.com/CMU-SAFARI/BLEND/blob/6f19e37f5120fc8d3b389978e7e4056c51810d0b/src/index.c#L83-L110) and [L230-L294](https://github.com/CMU-SAFARI/BLEND/blob/6f19e37f5120fc8d3b389978e7e4056c51810d0b/src/index.c#L230-L294). The bundled khash header documents quadratic probing: [`khash.h`](https://github.com/CMU-SAFARI/BLEND/blob/6f19e37f5120fc8d3b389978e7e4056c51810d0b/src/khash.h#L51-L58).

The primary paper is Firtina et al., **“BLEND: A Fast, Memory-Efficient, and Accurate Mechanism to Find Fuzzy Seed Matches in Genome Analysis,” NAR Genomics and Bioinformatics (2023)**, [PMC9853099](https://pmc.ncbi.nlm.nih.gov/articles/PMC9853099/). The repository is MIT-licensed: [LICENSE](https://github.com/CMU-SAFARI/BLEND/blob/6f19e37f5120fc8d3b389978e7e4056c51810d0b/LICENSE).

BLEND is a good validation target because it distinguishes two possible benefits:

- storage-table effects from TracEon;
- seed-generation effects from BLEND.

The benchmark must avoid attributing a result from the latter to the former.

---

## 4. Structurally similar but non-exact candidates

### StrobeMap / strobemers

The C++ StrobeMap implementation in [ksahlin/strobemers](https://github.com/ksahlin/strobemers) uses bundled `robin_hood::unordered_map`, not khash. Its active index is:

```cpp
robin_hood::unordered_map<uint64_t,
                          std::tuple<uint64_t, unsigned int>>
```

The key is a 64-bit strobemer hash; the value is `(offset in flat vector, occurrence count)`, not a single packed `uint64_t`. See [`index.hpp`](https://github.com/ksahlin/strobemers/blob/be7fa7a3bed408063dc5b8ff0c9744ea72850e9f/strobemers_cpp/index.hpp#L16-L57) and [`main.cpp`](https://github.com/ksahlin/strobemers/blob/be7fa7a3bed408063dc5b8ff0c9744ea72850e9f/strobemers_cpp/main.cpp#L660-L677). The bundled table is robin-hood open addressing: [`robin_hood.h`](https://github.com/ksahlin/strobemers/blob/be7fa7a3bed408063dc5b8ff0c9744ea72850e9f/strobemers_cpp/robin_hood.h#L8-L12).

This is **not** a direct TracEon C API drop-in:

- value is a tuple containing at least 64+32 bits of information;
- there is no documented low-bit singleton flag matching minimap2;
- the table and flat occurrence vector are coupled;
- the root source license was not found in the inspected snapshot, while the Bioconda recipe labels `strobemap` GPL-3.0. Treat licensing as unresolved until clarified.

The primary paper is **“Effective sequence similarity detection with strobemers,” Genome Research (2021)**, [10.1101/gr.219477.116](https://genome.cshlp.org/content/31/11/2080).

### strobealign

Current `ksahlin/strobealign` is primarily Rust and its production index is **not a hash table**. It stores a sorted `Vec<RefRandstrobe>` plus `bucket_starts`; lookup performs a short linear scan or binary search in a hash-prefix bucket. See [`src/index.rs`](https://github.com/ksahlin/strobealign/blob/8ba69020e86905b249b51ede94ca6f992c555bc5/src/index.rs#L12-L87) and [lookup](https://github.com/ksahlin/strobealign/blob/8ba69020e86905b249b51ede94ca6f992c555bc5/src/index.rs#L126-L195).

The repository also retains legacy C++ sources, but the current production seeding index is the Rust sorted-vector/range design cited above. A TracEon table swap would therefore be a redesign, not a drop-in validation.

Strobealign is MIT-licensed and its [`CONTRIBUTING.md`](https://github.com/ksahlin/strobealign/blob/main/CONTRIBUTING.md) explicitly welcomes normal issue and PR contributions. It does not contain an explicit AI-use rule. The primary paper is **“Flexible seed size enables ultra-fast and accurate read alignment,” Genome Biology (2022)**, DOI [10.1186/s13059-022-02831-7](https://doi.org/10.1186/s13059-022-02831-7).

### HiFiMapper

HiFiMapper uses a C++ `std::unordered_map` with a custom hash/equality:

```cpp
std::unordered_map<uint64_t, uint64_t, Hash, KeyEqual> locator;
hash(key) = std::hash<uint64_t>()(key >> 1);
eq(a,b) = ((a >> 1) == (b >> 1));
```

See [`minimizer_engine.hpp`](https://github.com/staverm/hifimap/blob/7de7deae4ff0c9fb025c40b4d69c166450309c5/include/hifimap/minimizer_engine.hpp#L111-L132). Its construction stores `value<<1|1` for singleton entries and packs repeated-region offset/counts into a 64-bit value: [`minimizer_engine.cpp`](https://github.com/staverm/hifimap/blob/7de7deae4ff0c9fb025c40b4d69c166450309c5/src/minimizer_engine.cpp#L182-L195) and [multi-hit path](https://github.com/staverm/hifimap/blob/7de7deae4ff0c9fb025c40b4d69c166450309c5/src/minimizer_engine.cpp#L251-L264).

Semantically, this is close to TracEon's key contract. But `std::unordered_map` does not prove open addressing; the target libstdc++/libc++ implementation must be measured or inspected separately. It is a high-risk/low-priority C++ internal refactor because the project snapshot is old and the tool has low adoption relative to the minimap2 lineage.

The repository is MIT-licensed: [LICENCE](https://github.com/staverm/hifimap/blob/main/LICENCE). The official repository includes a project thesis, **“Rapid Alignment of High-Fidelity Sequencing Data”**, rather than a separate peer-reviewed HiFiMapper paper: [PDF](https://raw.githubusercontent.com/staverm/hifimap/main/hifimap_paper.pdf).

### LRA

LRA's production minimizer index is a sorted vector/range design, not khash or open addressing. `GenomeTuple` contains a 64-bit 2-bit k-mer and a 32-bit position; mapping uses sorted vectors and boundaries. The visible `unordered_map<Tuple,int>` is used in `CalculateMinimizerStats()` for statistics, not primary seed lookup. See [`MMIndex.h`](https://github.com/ChaissonLab/LRA/blob/6221610ddef555af76985a91160a9e336e3d9035/MMIndex.h#L46-L110) and [`TupleOps.h`](https://github.com/ChaissonLab/LRA/blob/6221610ddef555af76985a91160a9e336e3d9035/TupleOps.h#L35-L92).

LRA is C++ but its USC-RL v1.0 license restricts use to academic/non-commercial use unless separately licensed: [LICENSE.txt](https://github.com/ChaissonLab/LRA/blob/master/LICENSE.txt). The primary paper is **“lra: A long read aligner for sequences and contigs,” PLOS Computational Biology (2021)**, DOI [10.1371/journal.pcbi.1009078](https://doi.org/10.1371/journal.pcbi.1009078).

LRA is a valuable **negative architecture control**: if TracEon wins only on khash-like indexes but not on sorted range indexes, that is a more honest conclusion than claiming broad aligner acceleration.

### meryl

meryl is a k-mer counter/database, not an aligner. The inspected counting path uses a prefix-partitioned, packed suffix/value/label structure rather than an open-addressed `uint64 -> uint64` map. The core counting object explicitly collects and sorts suffixes in segments: [`merylCountArray.H`](https://github.com/marbl/meryl/blob/3b3866ca23001b6c1e5edcaa7ea8a465fb89a5e6/src/meryl2/merylCountArray.H#L27-L41). Its utility types include 128-bit k-mer data and 32-bit values: [`kmers-tiny.H`](https://github.com/marbl/meryl-utility/blob/f800fc4dada365a701122b5ac1ecb946eb402cb4/src/kmers-v2/kmers-tiny.H#L30-L38).

Do not integrate TracEon into meryl as a table replacement. Use meryl only as:

- a k-mer-counting comparator;
- a source of repetitive-k-mer lists for Winnowmap;
- a separate packed-array versus hash-table memory/throughput control.

Its visible license notice is public-domain/U.S.-Government-work oriented with file-specific exceptions and citation requirements: [README.licenses](https://github.com/marbl/meryl/blob/3b3866ca23001b6c1e5edcaa7ea8a465fb89a5e6/README.licenses).

---

## 5. Additional structural controls

Two non-requested projects are useful for interpretation, not first-wave aligner integration:

- **MoleMap** uses a genuine custom open-addressing k-mer index, but its table stores a 32-bit bucket representation and external `dir/pos/ref` arrays, not TracEon's packed `uint64 -> uint64` contract. It is GPLv3. See [paper](https://www.biorxiv.org/content/10.1101/2022.06.20.496811v2.full.pdf), current [index.h](https://github.com/kehrlab/molemap/blob/master/src/index.h#L1-L86), and current [GPLv3 LICENSE](https://github.com/kehrlab/molemap/blob/master/LICENSE). These are moving-branch links, so MoleMap is a comparator rather than a pinned integration target.
- **ParBLiSS/kmerhash** is a generic C++ Robin Hood/open-addressing table library with templated key/value widths. It can instantiate `uint64_t,uint64_t`, but does not carry minimap2's low-bit semantics and is not an aligner. It is a useful microbenchmark control. See [README](https://github.com/ParBLiSS/kmerhash/blob/327d0652c922235cc4bdf99ddb05bbbc30eec370/README.md), [implementation](https://github.com/ParBLiSS/kmerhash/blob/327d0652c922235cc4bdf99ddb05bbbc30eec370/include/kmerhash/hashmap_robinhood.hpp#L41-L57), and [Apache-2.0 license](https://github.com/ParBLiSS/kmerhash/blob/327d0652c922235cc4bdf99ddb05bbbc30eec370/LICENSE).

---

# Q2 — Drop-in ranking and key-type differences

## 6. Compatibility ranking

| Rank | Tool | Semantic similarity | Maintenance snapshot | Benchmark value | Integration judgment |
|---:|---|---|---|---|---|
| 1 | **Winnowmap** | Exact classic khash contract | Recent | High: repetitive long reads and a recognized minimap2 derivative | Lowest-risk high-value fork |
| 2 | **mm2-fast stock** | Exact classic khash contract | Moderate/older | High: optimized seeding/chaining exposes whether table gains survive an optimized pipeline | Low source risk; mode confounder |
| 3 | **minigraph** | Exact key/value/flag contract; different `khashl` API/probing | Recent | High: pangenome graph mapping | Medium source risk; high validation breadth |
| 4 | **BLEND** | Exact table contract; different seed-generation path | Older research snapshot | Medium-high: fuzzy seeds/accurate reads | Low table risk; algorithmic confounder |
| 5 | **HiFiMapper** | Same shifted key/equality and similar low-bit encoding | Old | Low-medium: independent C++ implementation | Medium-high internal refactor; not proven open addressing |
| 6 | **StrobeMap** | 64-bit key, but tuple value and no same flag | Older | Medium algorithmic interest, low reuse certainty | Not a C API drop-in; license unresolved |
| 7 | **strobealign** | No map; sorted vector/range lookup | Recent | High as a negative/control workload | No table replacement; Rust redesign required |
| 8 | **LRA** | Production index is sorted vectors; stats map is incidental | Older | Medium as negative architecture control | Not a candidate |
| 9 | **meryl** | Packed prefix/suffix arrays; 128-bit data | Recent | High as k-mer comparator, not aligner | Not a candidate |
| 10 | **minigraph-cactus** | Inherits minigraph externally | Recent wrapper | High deployment-level value | Target rebuilt minigraph binary, not Cactus table |

## 7. Key/value compatibility classes

### Class A — direct/adaptable with current API

- Winnowmap: `uint64 -> uint64`, exact `>>1` hash/equality, exact singleton flag.
- mm2-fast stock: same.
- minigraph core: `uint64 -> uint64`, same semantics; adapter must translate classic khash calls to khashl equivalents or a narrow TracEon wrapper.
- BLEND: same table shape, but validate its source-specific shifted minimizer keys.

### Class B — semantically similar but not direct

- HiFiMapper: `uint64 -> uint64`, same custom hash/equality and low-bit convention, but C++ class internals, region payload, sharding, and standard-library storage differ.
- StrobeMap: `uint64` key but tuple value; packing the tuple into one `uint64` is not generally lossless.
- MoleMap: custom open addressing but 32-bit bucket keys and external occurrence arrays.

### Class C — different production index

- strobealign: sorted `Vec`/`vector` plus bucket boundaries and binary search.
- LRA: sorted minimizer vectors/ranges; `unordered_map` only for statistics.
- meryl: prefix/suffix bit-packed arrays and sorted blocks.

---

# Q3 — Expected impact and measurement

## 8. Impact hypotheses

No direction or magnitude should be assumed before measurement.

| Metric | What may improve | What may regress / confound | How to interpret |
|---|---|---|---|
| **Index build wall time** | Contiguous flat-map storage and different probing may reduce table-build work | Both implementations reserve before insertion; many small per-bucket maps may make allocation/object overhead dominate. The current API has one C call per inserted key, so C-ABI/status overhead is part of the real build comparison unless a batch API is added | Attribute only the post-sort table-build phase, and report API-call overhead separately if instrumented |
| **Peak RSS** | Different metadata and capacity behavior may reduce flags/unused slots | `ankerl`'s indexing metadata and per-map allocations can be larger; transient reallocation can dominate | Measure `/usr/bin/time -v` peak RSS and logical entries separately |
| **Index bytes / load time** | Logical serialization remains compact and independent of physical table | Iteration order differs; reloading rebuilds the map; physical layout is not serialized | Compare logical `.mmi` records and load time, not byte-identical table layout |
| **Seeding throughput** | Better cache locality and flat value storage may reduce lookup latency | `kh_get` is generated static-inline C; TracEon adds an out-of-line C call, custom hash, freeze check, and C++ map lookup; `idx_hash` itself is only a shift | Measure exact hit/miss traces and API call rate; include no-LTO and LTO variants |
| **Total mapping time** | Could improve if seeding dominates the chosen workload | Chaining, DP, I/O, and SIMD base alignment may dominate; a lookup win can disappear end-to-end | Report mapping wall/CPU time and phase breakdown; do not infer total speed from microbenchmarks |
| **Thread scaling** | A map that is fully built and then frozen should support concurrent lookups | An atomic flag alone does not make an in-progress transition safe; any mutation/destruction overlap invalidates results | Enforce a single-writer build phase, explicitly freeze all buckets, then create workers; run TSan and a stress test separately |
| **Correctness** | Should be unchanged if logical records and flags are preserved | Wrong `>>1` equality, lost bit 0, incorrect packed offset/count, or unordered iteration can alter results | Require logical-index and output-equivalence gates before performance claims |

The most important negative hypothesis is that TracEon can be **faster as a generic flat hash map but slower in this application** because minimap2's hash is unusually cheap and its per-bucket maps are pre-sized. Generic hash-table benchmarks are not sufficient; the [udb3 benchmark](https://github.com/attractivechaos/udb3/blob/master/README.md) itself warns that results depend on workload and CPU.

## 9. Required A/B methodology

For each realistic candidate:

1. Pin the exact upstream commit, compiler, C/C++ standard, ankerl version, optimization flags, ISA flags, allocator, and OS.
2. Maintain two build targets from the same tree:
   - `stock`: original khash/khashl or original C++ map;
   - `traceon`: only the table backend changed, with the same seed generation, chaining, alignment, output, and CLI parameters.
3. For mm2-fast, make three explicit modes: stock khash, TracEon stock-seeding, and LISA_HASH. Do not combine them.
4. Use the same reference/read files, checksums, preset, `-k/-w/-H/-I`, thread count, output sink, and index mode.
5. Randomize A/B order. Run exactly three measured repetitions for the primary median-of-3 result and report all three raw values. Before looking at A/B deltas, predeclare an extension rule: if either build has MAD/median >5%, repeat both variants to five runs and report median/MAD of all five. Treat a median improvement below 5% as “no engineering win claimed” unless a later powered experiment supports it. Use a warm-cache and separately labelled cold-cache series, not an unlabelled mixture.
6. Capture:
   - wall time and CPU time;
   - peak RSS;
   - output/index file sizes;
   - cycles, instructions, cache misses, branch misses with `perf stat` where available;
   - number of indexed entries, singleton fraction, occurrence histogram;
   - seed lookups, hits, candidate positions, and seeds/read if instrumentation is added;
   - API calls and failed status counts;
   - mapping output and truth metrics.
7. Report deltas as `100 * (traceon - stock) / stock`, with units and median/MAD. A result below run-to-run noise is “no detectable effect under this protocol,” not “zero overhead.”

### Minimum correctness gates

- exact singleton query: stored `base|1`, queried `base`, must hit and return the singleton position;
- repeated key: stored `base`, packed `(offset,count)`, must return all positions exactly;
- base-key collision: `base` and `base|1` must be treated as the same logical key;
- empty and absent buckets;
- max legal packed offset/count values;
- index iteration over occupied entries only;
- fresh build versus `.mmi` load;
- stock and TracEon logical entry-set equality;
- canonicalized PAF/SAM equivalence, allowing only explicitly documented record-order differences;
- simulated truth evaluation where available.

Tool-specific output oracles should be used rather than one universal equality rule: minimap2-family tools compare logical bucket entries, seed-hit traces, and canonicalized PAF/SAM; minigraph additionally compares GAF/graph-chain invariants; BLEND retains its fuzzy-seed invariants; strobealign/LRA remain native sorted-index controls rather than forced common-output targets.

## 10. Benchmark matrix

### Dataset tiers

| Tier | Reference / reads | Tools | Purpose |
|---|---|---|---|
| Smoke | Minimap2/minigraph small test fixtures, including MT-human/MT-orang and a tiny GFA | Winnowmap, mm2-fast, minigraph, BLEND | Build, load, output, and adapter debugging at low cost |
| Bacterial | **E. coli K-12** reference, accession [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3), plus PBSIM2-simulated ONT/PacBio reads | All four exact candidates; HiFiMapper if built | Fast repeated A/B runs; mapping truth; low-RSS environment |
| Human subset | Proposed fixed GRCh38 `chr1:1-50,000,000` half-open interval, or a separately labelled CHM13 equivalent; extract with `samtools faidx`, then record FASTA/checksum; simulate reads with a pinned PBSIM2 commit/configuration | Winnowmap, mm2-fast, BLEND; minigraph if represented as graph | Table occupancy and cache behavior at a larger scale without whole-genome cost |
| Human accurate reads | HG001/NA12878 or HG002 PacBio HiFi subset, exact files/checksums from [GIAB](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/) or an equivalent public release | Winnowmap, mm2-fast, BLEND, optional HiFiMapper | Accurate long-read workload and output/truth assessment |
| Repetitive stress | CHM13/T2T repeat-rich interval or segmental-duplication subset; keep reference and read subset fixed | Winnowmap first; mm2-fast/BLEND second | Tests occurrence distribution, singleton fraction, and high-hit lookup behavior |
| Graph | Minigraph's small graph fixtures first; then a pinned HPRC/minigraph-cactus graph/query subset | minigraph, then Cactus deployment wrapper | Tests khashl adapter and external-binary integration; report graph construction and graph mapping separately |
| Sorted/control | strobealign's documented baseline dataset and LRA/HiFiMapper only if their native builds are retained | No TracEon swap initially | Demonstrates whether observed effects are specific to hash-table backends |

Use minimap2's [cookbook](https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/cookbook.md), which documents E. coli/PBSIM2 examples, prebuilt-index workflows, and `paftools.js mapeval`, as the starting reproducibility template. For simulated datasets, fix simulator version, random seed, read count, read-length distribution, and error model. For real reads, record accession, release, platform/basecaller, subset rule, and SHA-256.

### Tool-by-metric matrix

| Tool | Index build | Peak RSS | Seed lookup rate | Total mapping | Correctness/truth |
|---|---:|---:|---:|---:|---:|
| Winnowmap | Required | Required | Required | Required | PAF/SAM + repeat-rich truth |
| mm2-fast stock/TracEon | Required | Required | Required | Required | Exact output and separate LISA comparison |
| minigraph | Required | Required | Required | Required | GAF/PAF logical graph mapping |
| BLEND | Required | Required | Required | Required | Fuzzy-seed outputs + accurate-read truth |
| HiFiMapper | Optional exploratory | Required if run | Required | Required | Native output; layout-dependent result |
| StrobeMap | No TracEon drop-in | Optional | Native map lookup only | Optional | Control only; license first |
| strobealign/LRA/meryl | Native baseline only | Optional | Native lookup/index rate | Optional | Negative architecture controls |

### Instrumentation boundary

Measure both:

- a **table/API microbenchmark** driven by recorded minimap2-style `(query_key, expected hit/miss)` traces; and
- the real aligner's **full mapping path**.

This is necessary because khash lookup is generated inline C while a C-to-C++ shim may be out-of-line. A table microbenchmark can show the backend's behavior; only the full aligner run shows whether that behavior matters to users.

---

# Q4 — AI-generated code and PR policies

## 11. Findings

In a bounded inspection on 2026-08-12, I located no explicit AI rule in the named repositories' pinned README, LICENSE, CONTRIBUTING/DEVELOPMENT, SECURITY, and code-of-conduct files where present. This means “not located in those files,” not repository-wide permission or hostility and not permission to submit code.

| Project | Visible policy evidence | Conclusion |
|---|---|---|
| minimap2 | README/developer guidance and code-of-conduct surface; no explicit AI/Copilot/LLM rule located | No explicit policy found; ask maintainer before collaboration |
| Winnowmap | README, license, source tree; no visible `CONTRIBUTING.md` or AI rule | No explicit policy found |
| minigraph | README, code of conduct, source tree; no explicit AI rule | No explicit policy found |
| minigraph-cactus | README, DEVELOPMENT, license/tree; no explicit AI rule located | No explicit policy found |
| mm2-fast | README, code of conduct, source tree; no explicit AI rule | No explicit policy found |
| BLEND | README, code of conduct, license; no explicit AI rule | No explicit policy found |
| strobemers/StrobeMap | README/source tree; root license itself unresolved; no explicit AI rule | No explicit AI policy; licensing is the larger blocker |
| strobealign | [`CONTRIBUTING.md`](https://github.com/ksahlin/strobealign/blob/main/CONTRIBUTING.md) explicitly welcomes ordinary issues/PRs; no AI rule | Most clearly open to normal collaboration, but AI provenance still should be disclosed if asked |
| HiFiMapper | README/license/source; no explicit AI rule | No explicit policy found |
| LRA | README/license/source; no explicit AI rule | No explicit policy found; license restrictions matter |
| meryl | README/license/source; no explicit AI rule | No explicit policy found |

A strobealign issue contains an individual contributor's Claude Code disclosure ([issue #568](https://github.com/ksahlin/strobealign/issues/568)); that is not repository policy. Silence in the other repositories is not evidence of either acceptance or hostility. Since TracEon will use forks/validation branches rather than upstream PRs, the practical action is to preserve source provenance, keep changes isolated, and ask before proposing collaboration.

---

# Q5 — ROI and integration sequence

## 12. Recommended 3–5 tools beyond minimap2

| Order | Tool | Effort | Expected impact to measure | Risk | What success would prove |
|---:|---|---|---|---|---|
| 1 | **Winnowmap** | Low–medium | Index RSS/build and seed lookup under weighted minimizers; end-to-end repetitive long-read mapping | Medium: C/C++ build and extra Winnowmap logic; license is source-specific/joint-work | TracEon's table behavior generalizes from canonical minimap2 to an independently maintained, accuracy-focused derivative |
| 2 | **mm2-fast stock** | Low | Whether an optimized minimap2 pipeline still benefits after table replacement; separate build and mapping phases | Medium: C++/AVX build, explicit linker, LISA_HASH confounder | Any gain is not merely an artifact of unoptimized minimap2; any regression exposes C-ABI overhead under a high-performance pipeline |
| 3 | **minigraph** | Medium | Khashl-vs-TracEon table behavior in sequence-to-graph indexing and mapping | Medium–high: khashl API differs; other khashl maps must remain untouched; graph workload is different | Generalization to lh3's graph mapper and a linear-probing khash variant, not just classic khash forks |
| 4 | **BLEND** | Low–medium | Table effects under fuzzy/BLEND seed values and accurate long reads | Medium: older research fork and changed seed generation; output semantics differ | TracEon preserves the table contract despite a changed seeding algorithm; separates storage from hash-generation effects |
| 5 | **HiFiMapper (optional)** | Medium–high | Direct C++ replacement of an independent `std::unordered_map` backend; memory/layout sensitivity | High relative to value: old project, low adoption, implementation-dependent STL layout, dependencies | A semantic-control result: whether TracEon can replace a custom `>>1` C++ map beyond the C minimap2 family |

**Prioritized integration sequence:**

1. **Winnowmap first.** It is the best balance of exact compatibility, real mapper value, and a workload where occurrence distributions are interesting.
2. **mm2-fast second.** It is mechanically close and tests whether the result survives a performance-optimized minimap2 derivative. Keep stock and LISA modes separate.
3. **minigraph third.** It broadens the story to graph mapping and khashl linear probing, with manageable but real adapter work.
4. **BLEND fourth.** It tests specialized fuzzy-seed generation after the adapter is stable.
5. **HiFiMapper only if an independent C++ control is worth the integration cost.**

StrobeMap, strobealign, LRA, and meryl should remain controls/read-only audits until the exact candidates produce trustworthy results.

## 13. What each successful result would and would not establish

- **Winnowmap success:** evidence that the logical `uint64 -> uint64` table swap is not specific to minimap2's default minimizer distribution. It does not establish universal long-read mapper acceleration.
- **mm2-fast success:** evidence that a seeding-table change can survive a more optimized pipeline. It does not compare against LISA_HASH unless that mode is run separately.
- **minigraph success:** evidence that TracEon's storage backend can support an lh3 graph index with khashl-like semantics. It does not establish that all graph indexes are compatible.
- **BLEND success:** evidence that storage replacement remains correct under changed seed values and specialized fuzzy matching. It does not isolate storage from BLEND's algorithmic effects without a stock-BLEND control.
- **HiFiMapper success:** evidence about an independent C++ map with similar equality semantics. It does not prove HiFiMapper's original map was open addressing or that its total mapper is popular/maintained.

---

# Honest risks and caveats

1. **Same logical type is not same workload.** Two tools may both say `uint64_t -> uint64_t` while differing in bucket partitioning, occurrence distribution, singleton fraction, and lookup hit/miss ratio.
2. **`idx_hash` is exceptionally cheap.** Replacing a one-shift hash with a stronger generic hash may hurt even if the new table has better locality.
3. **The C ABI may dominate.** khash operations are macro/static-inline C; an out-of-line C-to-C++ call and freeze check can erase backend gains, especially on short reads or small buckets.
4. **Many small maps matter.** Minimap2-family indexes use one map per bucket. Whole-genome RSS may be driven by map-object and allocator overhead rather than slot payload alone.
5. **Freeze is a lifecycle contract, not a magic lock.** Build/reserve/insert must complete before read workers; destruction must wait for all readers and returned pointers.
6. **Serialization is logical.** A TracEon-backed table will not preserve khash slot order. Compare logical records and mapping outputs, not physical table bytes.
7. **Different aligners have different bottlenecks.** A seed lookup improvement may not change total mapping time when chaining, DP, SIMD alignment, or I/O dominates.
8. **Fork drift is real.** A patch against one minimap2-derived repository cannot be mechanically generalized to all forks. Pin every source revision and inspect the actual index path.
9. **Licenses differ.** Winnowmap has source-specific joint-work notices; LRA is USC-RL academic/non-commercial; StrobeMap's root license is unresolved; external dependencies need a separate license inventory.
10. **Maintenance labels are snapshots.** Commit dates show activity at inspection, not support commitments.
11. **No performance numbers exist yet.** Until the repaired forks run the matrix, all impact claims are hypotheses.

---

# Recommended aligner validation matrix

| Phase | Tool(s) | Reference/read tier | Primary metrics | Gate |
|---|---|---|---|---|
| 0 | minimap2 baseline + TracEon | Tiny fixtures and E. coli | Build/load time, RSS, exact logical entries, mapping output | API and lifecycle correctness |
| 1 | Winnowmap | E. coli simulated; human chr1 subset; repeat-rich CHM13 subset; HiFi subset | Index build/RSS, lookup rate, mapping wall/CPU, output/truth | First generalization result |
| 2 | mm2-fast stock + LISA comparator | Same bacterial/human datasets | Same metrics, plus ISA/build mode and C++ link overhead | Optimized-pipeline sensitivity |
| 3 | minigraph | Tiny GFA, then pinned graph/query subset | Graph index build/RSS, lookup rate, GAF/PAF wall time, output equivalence | khashl/graph validation |
| 4 | BLEND | Accurate-read E. coli and human subset | Build/RSS, lookup rate, total time, fuzzy-seed correctness | Specialized-seed validation |
| Control | strobealign, LRA, meryl; optional StrobeMap/HiFiMapper | Native recommended datasets | Native index/mapping profile only | Avoid overgeneralization |

For every A/B cell: same machine, same thread count, same compiler/flags/allocator, same input checksums, randomized order, exactly three primary repetitions with the predeclared MAD extension rule, **median-of-3 plus all raw values**, and explicit warm/cold-cache status. Capture `uname -a`, `lscpu`, compiler versions, git SHAs, allocator/runtime, CPU affinity/NUMA binding, filesystem, and command lines. Define lookup rate as instrumented lookup calls/sec, hits/sec, and cycles/lookup; collect `perf stat` cycles, instructions, cache-misses, and branch-misses when available.

---

# Sources

## Core table/index sources

- [minimap2 `index.c`, v2.31](https://github.com/lh3/minimap2/blob/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/index.c)
- [minimap2 `khash.h`, v2.31](https://github.com/lh3/minimap2/blob/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/khash.h)
- [Winnowmap `index.c`](https://github.com/marbl/Winnowmap/blob/72403b350fafd3aec98d1a6c67123d0d0d1dff7c/src/index.c)
- [Winnowmap vendored `khash.h`](https://github.com/marbl/Winnowmap/blob/72403b350fafd3aec98d1a6c67123d0d0d1dff7c/src/khash.h)
- [minigraph `index.c`](https://github.com/lh3/minigraph/blob/2f569ebe3071fcb242f1bf0eabbb917941b47239/index.c)
- [minigraph `khashl.h`](https://github.com/lh3/minigraph/blob/2f569ebe3071fcb242f1bf0eabbb917941b47239/khashl.h)
- [mm2-fast `index.c`](https://github.com/bwa-mem2/mm2-fast/blob/14fe36c100f6c2aab224d000f3903ca5909640cd/index.c)
- [mm2-fast `seed.c`](https://github.com/bwa-mem2/mm2-fast/blob/14fe36c100f6c2aab224d000f3903ca5909640cd/seed.c)
- [BLEND `index.c`](https://github.com/CMU-SAFARI/BLEND/blob/6f19e37f5120fc8d3b389978e7e4056c51810d0b/src/index.c)
- [BLEND vendored `khash.h`](https://github.com/CMU-SAFARI/BLEND/blob/6f19e37f5120fc8d3b389978e7e4056c51810d0b/src/khash.h)
- [StrobeMap `index.hpp`](https://github.com/ksahlin/strobemers/blob/be7fa7a3bed408063dc5b8ff0c9744ea72850e9f/strobemers_cpp/index.hpp)
- [StrobeMap `main.cpp`](https://github.com/ksahlin/strobemers/blob/be7fa7a3bed408063dc5b8ff0c9744ea72850e9f/strobemers_cpp/main.cpp)
- [strobealign `src/index.rs`](https://github.com/ksahlin/strobealign/blob/8ba69020e86905b249b51ede94ca6f992c555bc5/src/index.rs)
- [HiFiMapper `minimizer_engine.hpp`](https://github.com/staverm/hifimap/blob/7de7deae4ff0c9fb025c40b4d69c166450309c5/include/hifimap/minimizer_engine.hpp)
- [HiFiMapper `minimizer_engine.cpp`](https://github.com/staverm/hifimap/blob/7de7deae4ff0c9fb025c40b4d69c166450309c5/src/minimizer_engine.cpp)
- [LRA `MMIndex.h`](https://github.com/ChaissonLab/LRA/blob/6221610ddef555af76985a91160a9e336e3d9035/MMIndex.h)
- [meryl `merylCountArray.H`](https://github.com/marbl/meryl/blob/3b3866ca23001b6c1e5edcaa7ea8a465fb89a5e6/src/meryl2/merylCountArray.H)
- [meryl-utility k-mer types](https://github.com/marbl/meryl-utility/blob/f800fc4dada365a701122b5ac1ecb946eb402cb4/src/kmers-v2/kmers-tiny.H)
- [Cactus minigraph wrapper](https://github.com/ComparativeGenomicsToolkit/cactus/blob/68252c2a9406bbf82b62434e7060403d44bed5b8/src/cactus/refmap/cactus_minigraph.py)
- [TracEon public repository](https://github.com/woosflex/TracEon) — current local API evidence is `include/kmerindex_c_api.h` and `src/kmerindex_c_api.cpp` at unpushed commit `9c54742a7699e0701c0d836ba68c812ed510a6fe`
- [ankerl::unordered_dense](https://github.com/martinus/unordered_dense)

## Papers and benchmark/data sources

- [Winnowmap2, Nature Methods (2022)](https://doi.org/10.1038/s41592-022-01457-8)
- [Minigraph-Cactus, Nature Biotechnology (2023)](https://doi.org/10.1038/s41587-023-01793-w)
- [BLEND, PMC9853099](https://pmc.ncbi.nlm.nih.gov/articles/PMC9853099/)
- [Strobemers, Genome Research (2021)](https://genome.cshlp.org/content/31/11/2080)
- [Strobealign, Genome Biology (2022)](https://doi.org/10.1186/s13059-022-02831-7)
- [HiFiMapper project thesis](https://raw.githubusercontent.com/staverm/hifimap/main/hifimap_paper.pdf)
- [LRA, PLOS Computational Biology (2021)](https://doi.org/10.1371/journal.pcbi.1009078)
- [Minimap2 paper, arXiv:1708.01492](https://arxiv.org/abs/1708.01492)
- [Minimap2 cookbook](https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/cookbook.md)
- [GIAB HG001/NA12878 GRCh38 directory](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/)
- [T2T CHM13 repository](https://github.com/marbl/CHM13)
- [E. coli K-12 NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3)
- [Unordered Dictionary Benchmark 3](https://github.com/attractivechaos/udb3/blob/master/README.md)

## Policy/licensing sources

- [strobealign CONTRIBUTING](https://github.com/ksahlin/strobealign/blob/main/CONTRIBUTING.md)
- [strobealign issue #568](https://github.com/ksahlin/strobealign/issues/568)
- [StrobeMap Bioconda recipe](https://bioconda.github.io/recipes/strobemap/README.html)
- [LRA license](https://github.com/ChaissonLab/LRA/blob/master/LICENSE.txt)
- [MoleMap current GPLv3 license](https://github.com/kehrlab/molemap/blob/master/LICENSE)
