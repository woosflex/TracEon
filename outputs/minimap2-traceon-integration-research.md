# TracEon `traceon_kmer` as a minimap2 `khash_t(idx)` replacement

**Research status:** source-backed design review; no production integration or end-to-end performance result has been produced. Local compile probes are reported explicitly below.

**Observation scope:** GitHub refs, release metadata, and pinned source files were queried during this research run; time-dependent claims below are anchored to the cited commit/tag rather than to a moving `master` file. The local checkout is TracEon commit `9566b24dfee04e7d582fa082744cf37c8d59b99b`; GCC is `16.1.1`. [33]

**Current baseline observed:** the checkout's experimental copy is minimap2 **2.28-r1209** (`benchmarks/minimap2_src/minimap.h`). The latest canonical release observed through the official release/API sources is **v2.31/r1302**, commit `3c28777e7e2dcc90f825de1b9f17a89cca7d4452`. [1–3, 33]

## Executive recommendation

Build the first proof as an independent **`mm2-traceon` fork of minimap2 v2.31**, pinned to the full upstream commit SHA, while retaining stock khash as the default build. Use TracEon's C API as a narrow adapter, compile its C++ implementation with `$(CXX)`, and drive the final executable link with `$(CXX)` (or a separate `CXXLD` variable). Do not start by rewriting the hash table in C or by modifying upstream's public API.

Before benchmarking, fix three integration blockers in the current experimental copy:

1. Replace the three old shared-handle iterator call patterns with the current caller-owned iterator ABI.
2. Use one coherent `kmerindex_c_api.h`/`.cpp` pair; the copied benchmark `.cpp` is stale and conflicts with the tracked header.
3. Check `create`, `reserve`, and `insert` results explicitly. `assert(inserted)` is incorrect for the documented `-1` exception status because `assert(-1)` succeeds.

There are two production-level C-API issues to resolve before mapping workers run. First, the current C++ shim stores `frozen` as a plain mutable `bool` and writes it from `kmerindex_get()`. Minimap2 performs concurrent lookups, so the first concurrent lookup phase can race on that flag. Second, the advertised no-throw error model is not actually total: `set_error()` assigns to a thread-local `std::string` from `noexcept` functions, and that assignment can allocate and terminate the process if it throws. Replace diagnostics with a fixed-size non-allocating thread-local buffer (or catch diagnostic failure with a non-allocating fallback).

For synchronization, finish all inserts/reserves, validate the index, and call `kmerindex_freeze()` on every non-null bucket after generation/load and before creating mapping workers. Make the state atomic with release/acquire operations, and document the stronger lifetime contract: no insert/reserve/destroy or iterator invalidation may overlap any read, and the index must remain alive until all returned pointers/iterators are done. Prefer explicit pre-freeze plus an atomic guard; an atomic flag alone does not make concurrent map mutation safe.

This is better framed initially as a **drop-in experimental fork / optional backend**, not as an upstream minimap2 PR. Upstream's default build is C + GNU make + zlib; adding a C++ static dependency and a new container requires a compelling, independently reproducible end-to-end benefit and a platform-safe default-off path.

---

## Q1. Current minimap2 state

### Canonical repository and release

The canonical project is [**`lh3/minimap2`**](https://github.com/lh3/minimap2). The GitHub repository metadata reports `fork: false`, `archived: false`, default branch `master`, and recent pushes/releases. The latest release observed is **[Minimap2-2.31 (r1302)](https://github.com/lh3/minimap2/releases/tag/v2.31)**, published 19 May 2026, at commit [`3c28777e7e2dcc90f825de1b9f17a89cca7d4452`](https://github.com/lh3/minimap2/commit/3c28777e7e2dcc90f825de1b9f17a89cca7d4452). The upstream README's binary installation command also names v2.31. [1–3, 10]

This is not a maintenance-only or archived project. The release includes recent bug fixes, including supplementary/secondary flag handling, a rare Smith–Waterman out-of-bounds issue, and a mappy use-after-free. The release notes explicitly warn that the fixes can make v2.31 occasionally produce alignments different from the prior release. Therefore, pin both the source SHA and expected correctness behavior; do not treat a version bump as performance-only. [1–3]

### License

The repository's pinned [`LICENSE.txt`](https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/LICENSE.txt) is the **MIT License**. Preserve the copyright and permission notices in a fork and in binary/source redistribution. [14]

### Maintenance and API implications

The v2.31 README's [Developer's Guide](https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/README.md) says minimap2 aims to keep APIs in `minimap.h` stable, while `mmpriv.h` contains private APIs that may change frequently. The `idx` khash instantiation in `index.c` is even more private than the documented C API. That makes a source patch feasible but not upstream-stable. [10]

### AI-generated code / AI PR policy

I found no tracked `CONTRIBUTING.md` and no AI/LLM/Copilot policy in the checked upstream repository tree, README, NEWS, or code of conduct. This is a bounded repository-tree/search observation, not proof of a global policy absence. The absence of a policy is **not** evidence that AI-assisted contributions are either permitted or prohibited. For an AI-assisted PR, ask the maintainer or disclose provenance according to TracEon's own contribution policy before submitting. [1, 10, 33]

### Is there a de-facto maintained fork?

There are maintained or at least active derivatives, but no evidence that one has replaced `lh3/minimap2` as the general long-read mapping standard:

- The [VGP assembly pipeline](https://github.com/VGP/vgp-assembly/tree/master/pipeline/purge_dups) invokes `minimap2`; its [asset documentation](https://github.com/VGP/vgp-assembly/tree/master/pipeline/asset) records a historical validation against v2.17. This supports the command/tool lineage, not a current universal version pin. [15, 16]
- The [nf-core minimap2 module](https://github.com/nf-core/modules/tree/master/modules/nf-core/minimap2/align) is a maintained workflow module that packages/configures minimap2; it is not a minimap2 fork. [17]
- [`bwa-mem2/mm2-fast`](https://github.com/bwa-mem2/mm2-fast) is a real specialized fork with AVX2/AVX-512, chaining/alignment optimizations, and optional learned-index seeding. It is not the de-facto standard and, importantly for Q5, its `index.c` retains minimap2's khash instantiation. [18, 19]
- [`Minimap2onGPU/mm2-gb`](https://github.com/Minimap2onGPU/mm2-gb) is a specialized GPU derivative, not a general CPU replacement. [36]

**Recommendation:** target `lh3/minimap2`; use workflow forks only as benchmark comparators or if their explicit specialization is relevant. [15–19]

---

## Q2. Which minimap2 line to target?

### Exact compatibility target

At upstream v2.31, [`index.c` lines 17–20](https://github.com/lh3/minimap2/blob/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/index.c#L17-L20) contains the same semantics:

```c
#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;
```

The current raw upstream source omits the optional semicolon after the macro invocation; that formatting difference is irrelevant to the API semantics described in the question. [4]

The exact three-line pattern is present in the checked tags:

| Release | Commit | Result |
|---|---|---|
| v2.26 | [`e28a55be`](https://github.com/lh3/minimap2/blob/e28a55be86b298708a2a67c924d665a00b8d829c/index.c#L18-L20) | exact `idx_hash`/`idx_eq`/`KHASH_INIT(idx)` semantics |
| v2.27 | [`b6762503`](https://github.com/lh3/minimap2/blob/b6762503a9c9af34ae02f56784035ff852ed04a5/index.c#L18-L20) | same |
| v2.28 | [`8170693d`](https://github.com/lh3/minimap2/blob/8170693de39b667d11d8931d343c94a23a7690d2/index.c#L18-L20) | same |
| v2.29 | [`1fd85be6`](https://github.com/lh3/minimap2/blob/1fd85be6e2515c9194740e1d2e6a2625be36f508/index.c#L19-L21) | same |
| v2.30 | [`79c9cc18`](https://github.com/lh3/minimap2/blob/79c9cc186b95f50bd899f69b48eba995ced810c6/index.c#L19-L21) | same |
| v2.31 | [`3c28777e`](https://github.com/lh3/minimap2/blob/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/index.c#L17-L20) | same |

The six rows were checked with `git ls-remote --tags` plus `git show <tag>:index.c | grep -E 'idx_hash|idx_eq|KHASH_INIT\\(idx'`; the full SHAs and source links above make the result reproducible. This is empirical stability, not a future compatibility promise. [4–9]

### Recommendation

Fork from **v2.31/r1302**, not v2.28, because the semantics match and v2.31 is the latest stable upstream release. Pin the full SHA and add a source sentinel in CI that checks both:

```text
idx_hash(a) == a >> 1
idx_eq(a,b) == ((a >> 1) == (b >> 1))
```

If the immediate objective is only to bring up the existing local copy, v2.28-r1209 is a valid staging target: the local file identifies that version and the pinned v2.28 source has matching semantics. However, do not use that copied source as the production baseline without porting the patch to v2.31 and re-running correctness tests. The v2.31 source differs beyond the `idx` block, so a patch copied mechanically from v2.28 is not a version upgrade strategy. [7, 33, 4]

---

## Q3. Integration surface and operation mapping

The v2.31 source does **not** contain `mm_idx_prep`. The construction path is: [4]

```text
mm_idx_gen -> mm_idx_post -> worker_post
mm_idx_build -> mm_idx_gen
mm_idx_str -> mm_idx_post
```

The relevant source is [upstream `index.c`](https://github.com/lh3/minimap2/blob/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/index.c); the specific operation forms are also visible in the [v2.30 source](https://github.com/lh3/minimap2/blob/79c9cc186b95f50bd899f69b48eba995ced810c6/index.c), which is structurally close to the local experimental copy.

### Complete `idx` inventory

| Context | Upstream khash operations | TracEon mapping | Required notes |
|---|---|---|---|
| `mm_idx_destroy` | `kh_destroy(idx, h)` per bucket | `kmerindex_destroy(h)` | Null-safe destruction is required. The separate sequence-name `khash_t(str)` remains unchanged. |
| `mm_idx_get` | `kh_get`, `kh_end`, `kh_key`, `kh_val` | `kmerindex_get(h, query, &matched_key)` | Return the stored value pointer and actual stored key. `matched_key & 1` replaces `kh_key(h,k)&1`; `*val` replaces `kh_val(h,k)`. |
| `mm_idx_stat` | `kh_size`; loop `0..kh_end`; `kh_exist`; `kh_key`; `kh_val` | `kmerindex_size`; `kmerindex_iter_begin(h,&it)`; `kmerindex_iter_next(&it,...)` | The iterator must yield occupied entries only. Sum singleton as 1, multi-hit value low 32 bits as count. |
| `mm_idx_cal_max_occ` | `kh_size`; occupied-slot iteration; key/value reads | same size + caller-owned iterator | Build the occurrence-count array from exact singleton/multi encoding. |
| `worker_post` | `kh_init`; `kh_resize`; `kh_put`; post-insert `kh_key |= 1`; `kh_val = ...` | `kmerindex_create`; `kmerindex_reserve`; one-shot `kmerindex_insert(final_key,value)` | Pre-grouped insertion makes one-shot final-key insertion valid. Singleton: `base|1 -> position`; multi: `base -> (start_p<<32)|n`. |
| `mm_idx_dump` | direct logical size; occupied-slot iteration; `kh_key`; `kh_val` | `kmerindex_size`; caller-owned iterator | Serialize logical key/value pairs, including singleton bit. Do not serialize ankerl's raw table layout. |
| `mm_idx_load` | `kh_init`; `kh_resize`; `kh_put`; `kh_val` assignment | create; reserve; insert serialized final key/value | The existing `.mmi` logical format can remain compatible if the patched code writes/reads the same records. |
| `mm_idx_reader_read` | indirect `mm_idx_load`, `mm_idx_dump`, or generation | covered above | Add explicit freeze after load/generation before concurrent mapping. |

The operation inventory above is from the pinned upstream source; the TracEon API signatures and freeze/iterator contracts are from the pinned header and implementation. [4, 25, 26]

### Operations not needed for `idx`

No canonical `index.c` call site uses `kh_del(idx)`, `kh_clear(idx)`, `kh_foreach`, `kh_foreach_value`, or `kh_begin(idx)`. It uses a literal slot loop for inspection. The C API's iterator replaces the whole occupied-slot loop. `kh_n_buckets` is not part of the logical dump contract.

The documented invariant is useful: `worker_post` sorts each bucket by minimizer and groups identical minimizers before inserting once per group. Therefore the replacement may form the final key before insertion and avoid unstable in-place key mutation. This does **not** eliminate duplicate detection or failure handling: `kmerindex_insert()` must still be checked. [4, 25, 26]

### Current local patch blockers

The local experimental `benchmarks/minimap2_src/index.c` is not build-ready against the tracked C API. It calls:

```c
kmerindex_iter_begin(h);
kmerindex_iter_next(h, &kk, &vv);
```

but the header requires:

```c
kmerindex_iter_t it;
kmerindex_iter_begin(h, &it);
kmerindex_iter_next(&it, &kk, &vv);
```

This occurs in `mm_idx_stat`, `mm_idx_cal_max_occ`, and `mm_idx_dump`. **Local probe, not an end-to-end build:** at TracEon commit `9566b24d`, GCC 16.1.1 was invoked as:

```sh
cd benchmarks/minimap2_src && \
 gcc -std=c99 -Wall -Wextra -Werror -DHAVE_KALLOC \
   -I../../include -I../../build/_deps/ankerl_unordered_dense-src/include \
   -c index.c -o /tmp/traceon-index.o
```

It exited nonzero; the primary diagnostics were “too few arguments” for `kmerindex_iter_begin` and an incompatible first argument to `kmerindex_iter_next` at the three call sites. The same command also surfaced unrelated pre-existing `-Werror` warnings in the copied v2.28 source. [33] The copied `benchmarks/minimap2_src/kmerindex_c_api.cpp` is stale: its reserve and iterator declarations conflict with the tracked header and it lacks the tracked exception translation, freeze, and caller-owned iterator behavior. [25, 26, 33]

The patched C source ignores `create()` and `reserve()` failures and treats `insert()` as a Boolean assertion. The API documents `insert()` as `1` inserted, `0` collision/frozen/null, and `-1` translated C++ exception. In particular, `assert(-1)` is true, so allocation failure can be mistaken for success. Use explicit status checks and fail before publishing an incomplete bucket. On failure, destroy the temporary handle, free any position array, leave `b->h` null, and propagate a fatal/index-build status rather than publishing a partial bucket. The load path should similarly destroy the temporary handle on short input, reserve failure, duplicate, or insert failure.

The fixed-size iterator object also needs a small ABI hardening pass: assert both `sizeof` and `alignof` suitability for placement-new storage, define what repeated `iter_begin()` on one object does, and either provide a cleanup operation or ensure the stored iterator state is trivially destructible. The public lifetime rule must remain explicit: the index cannot be destroyed while a caller-owned iterator or returned value pointer is live. [25, 26]

---

## Q4. Build and link feasibility

### Does the header compile as C99?

Yes, in the current checkout. The public header uses C headers, opaque structs, integer/pointer types, and an `extern "C"` block only under `__cplusplus`; `TRACEON_CAPI_NOEXCEPT` expands to nothing for C. I ran:

```sh
gcc -std=c99 -pedantic-errors -Wall -Wextra -Werror \
  -Iinclude -fsyntax-only /tmp/traceon_c99_probe.c
```

where the probe included `kmerindex_c_api.h` and instantiated `kmerindex_iter_t`. Result: `C99_HEADER_PASS`.

The header is therefore a valid C-facing surface. This does not prove the whole minimap2 adapter compiles; the current adapter's iterator calls fail as described above. [25, 33]

### Important correction to the advertised error model

The header documents that no C++ exception crosses the C ABI, but the current implementation does not fully establish that guarantee. `set_error()` writes to a thread-local `std::string` from functions declared `noexcept`; assigning a new diagnostic may allocate, and an allocation failure in that path invokes `std::terminate`. This can occur inside an exception handler that is attempting to translate another allocation failure. Treat the current status model as **not yet exception-safe under memory exhaustion**. Before integration, replace the diagnostic with a fixed-size thread-local character buffer and non-allocating copy/formatting, or catch all diagnostic failures and use a static fallback string. Add an allocation-failure test that exercises `create`, `reserve`, `insert`, `get`, and the error-reporting call itself. [26]

### Option (a): final link with a C++ driver — recommended

Keep all minimap2 `.c` files compiled by `CC`, compile the shim/library with `CXX`, and use a separate final linker variable:

```make
CC    ?= gcc
CXX   ?= g++
CXXLD ?= $(CXX)
CFLAGS   += -std=c99
CXXFLAGS += -std=c++20 -I/path/to/pinned/traceon/include \
            -I/path/to/pinned/ankerl_unordered_dense/include

traceon_kmer.o: /path/to/pinned/traceon/src/kmerindex_c_api.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Either add traceon_kmer.o to libminimap2.a, or link a pinned archive here.
minimap2: main.o libminimap2.a
	$(CXXLD) $(CXXFLAGS) main.o -o $@ -L. -lminimap2 \
	  -L/path/to/pinned/traceon/lib -ltraceon_kmer $(LIBS)
```

A C++17 standard is sufficient for the shown shim (`try_emplace` and structured bindings); C++20 is appropriate when consuming the TracEon target as specified in the question. Apply sanitizer and optimization flags consistently to C, C++, archive, and final-link commands. If the shim object is already inside `libminimap2.a`, omit `-ltraceon_kmer`. If it is a separate static archive, put it after `-lminimap2`; static archive ordering matters. Remove hard-coded `-lstdc++` when the C++ driver supplies the runtime, but record whether the tested runtime is GNU libstdc++ or libc++.

GCC's [official documentation](https://gcc.gnu.org/onlinedocs/gcc/Invoking-G_002b_002b.html) states that `gcc` does not add the C++ library whereas `g++` automatically links it. A concrete [mixed C/C++ static-library example](https://github.com/noahp/c-cpp-lib-example) compiles C with `gcc`, C++ with `g++`, archives the objects, and links the C executable with `gcc ... -lstdc++`. Both patterns are practical on GNU toolchains; the C++ driver is clearer because it selects the C++ runtime through the toolchain's normal driver rules. Cross-toolchain portability still requires explicit testing of the compiler, runtime, and target platform. [23, 24, 33]

The current local Makefile demonstrates the idea but contradicts itself: its comment says to link with g++, while the recipe still invokes `$(CC)` and appends `-lstdc++`. Make the linker choice explicit rather than relying on that accidental GNU-specific behavior.

### Option (b): pure-C reimplementation behind the same API

This would avoid a C++ runtime and preserve the simplest minimap2 build, but it is no longer using the TracEon `ankerl::unordered_dense` backend. It becomes a new C hash-table implementation that must reproduce:

- `key >> 1` hash/equality semantics;
- one-shot final-key insertion;
- stable value pointers after build;
- caller-owned iteration;
- reserve/freeze/error behavior.

This is reasonable if portability and upstreamability outweigh reuse of TracEon's C++ backend. It is not the best first experiment because it confounds the question “is this table faster?” with “did a second table implementation change the result?”

### Option (c): C-compatible subset of the C++ library

A header-only C++ map cannot be made C-compatible by changing the header alone. A C-compatible subset means writing a separate C implementation or a C ABI wrapper, which collapses into option (b) plus the C++ adapter approach. Keep C++ templates and exceptions private to `.cpp`; the current opaque-handle C API is the right architecture.

### Build recommendation

Use option (a) for the first fork:

1. Pin `traceon_kmer` and `ankerl::unordered_dense` versions.
2. Compile the C API implementation once as a C++ object or `libtraceon_kmer.a`.
3. Keep the public header C99-only.
4. Drive the final minimap2 link with `CXXLD`.
5. Preserve a stock-khash build target so default users do not acquire a C++ dependency.
6. Test GCC, Clang, x86-64, and ARM64/NEON before making any portability claim.

---

## Q5. Prior art

### Direct table replacement

I found no checked project that provides a published, end-to-end, drop-in replacement for minimap2's exact `KHASH_INIT(idx, ...)` table using an AVX2/SIMD or C++ hash map. This is a bounded search result, not proof that no private or obscure fork exists; the search covered the repositories and query angles recorded in the provenance sidecar. [20–22, 34]

Relevant prior art is adjacent rather than identical:

1. **mm2-fast:** [`bwa-mem2/mm2-fast`](https://github.com/bwa-mem2/mm2-fast) reports up to 1.8× speedup using AVX-512/AVX2 across seeding, chaining, and alignment, and cites the Nature Computational Science paper [Kalikar et al., 2022](https://doi.org/10.1038/s43588-022-00201-8). Direct inspection of its [`index.c`](https://github.com/bwa-mem2/mm2-fast/blob/master/index.c) shows it still includes `khash.h` and retains the `idx` instantiation. It is useful SIMD and benchmark methodology prior art, but not a khash replacement. [18, 19, 31]
2. **minimap2-pure-rs:** The Rust reimplementation's [`OPTIMIZATION.md`](https://github.com/henriksson-lab/minimap2-pure-rs/blob/main/OPTIMIZATION.md) reports a 6–7% CPU-cycle improvement in its own benchmark and describes a unified index representation that changes singleton/multi-hit storage and downstream position expansion. Do **not** describe this as a two-probe `mm_idx_get()` baseline: pinned upstream v2.31 `mm_idx_get()` performs one `kh_get()` and then interprets the stored key/value. The Rust notes discuss broader lookup/seed-collection structure, not evidence that upstream's `mm_idx_get()` itself probes twice. This remains useful algorithmic prior art: storage representation may matter more than the hash function. It is not a C++ drop-in backend. [20, 4]
3. **minimap2 index modifier:** [`ispras/minimap2_index_modifier`](https://github.com/ispras/minimap2_index_modifier) changes index generation to incorporate VCF variants and has a 2024 publication, but it does not establish a different per-bucket hash container. [21]
4. **Generic hash benchmarks:** The [C/C++ hash-table benchmark](https://github.com/JacksonAllan/c_cpp_hash_tables_benchmark/) and its [published report](https://jacksonallan.github.io/c_cpp_hash_tables_benchmark/) include khash and `ankerl::unordered_dense` and measure insertion, lookup, deletion, and iteration. Those are useful smoke tests, not evidence for minimap2: the minimap2 workload is per-bucket, pre-reserved, bit-0-normalized, and split between build-time mutation and frozen pointer-returning lookup. [22, 4, 25, 26]
5. **Recent mapper comparisons:** [`rammap` performance notes](https://github.com/jwanglab/rammap/blob/main/docs/performance.md) compare a reimplementation with minimap2 v2.31 on public human and long-read data, but the results cover a broad implementation and should not be attributed to a hash-table substitution. [30]

### Prior-art conclusion

The novelty claim should be narrow: **a C ABI adapter allowing a TracEon open-addressing map to replace minimap2's per-bucket `idx` table while preserving singleton-bit semantics and `.mmi` logical serialization**. Do not claim that minimap2 has never had index alternatives or that generic robin-hood tables are faster until the exact workload is measured.

---

## Q6. Validation and adoption path

### Official minimap2 tests and examples

The official v2.31 Makefile has targets `all`, `extra`, `clean`, and `depend`; it has **no `make test` target** ([pinned Makefile](https://raw.githubusercontent.com/lh3/minimap2/3c28777e2dcc90f825de1b9f17a89cca7d4452/Makefile)). [11] The README's canonical smoke test is: [10, 11]

```sh
./minimap2 -a test/MT-human.fa test/MT-orang.fa > test.sam
./minimap2 -x map-ont -d MT-human-ont.mmi test/MT-human.fa
./minimap2 -a MT-human-ont.mmi test/MT-orang.fa > test.sam
```

The pinned [`test/` tree](https://github.com/lh3/minimap2/tree/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/test) includes `MT-human.fa`, `MT-orang.fa`, `q-inv.fa`, `q2.fa`, `t-inv.fa`, `t2.fa`, and other regression fixtures. Upstream CI ([pinned `ci.yaml`](https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/.github/workflows/ci.yaml)) compiles with GCC/Clang and ARM/macOS variants but does not execute a mapping regression suite. A successful upstream `make` is therefore not enough. [7, 11, 13]

The pinned v2.31 [cookbook](https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/cookbook.md) supplies a useful reproducible simulated long-read path: [12]

```sh
minimap2 -ax map-pb -t4 ecoli_ref.fa ecoli_p6_25x_canu.fa > mapped.sam

# PBSIM2 recipe from the cookbook
src/pbsim --depth 1 --length-min 5000 --length-mean 20000 \
  --accuracy-mean 0.95 --hmm_model data/R94.model ../ecoli_ref.fa
paftools.js pbsim2fq ../ecoli_ref.fa.fai sd_0001.maf > ../ecoli_pbsim.fa
```

The read names encode truth positions, enabling `paftools.js mapeval` where the input format supports it.

### Human-scale dataset tiers

Use tiers rather than jumping directly to an expensive whole-genome run:

- **Tier 0: official smoke/regression.** Minimap2's six small fixtures, both direct-reference and `.mmi` index-first mapping.
- **Tier 1: simulated E. coli long reads.** The cookbook/PBSIM2 recipe. This gives truth-bearing reads and isolates mapping correctness.
- **Tier 2: NA12878/HG001.** The [Garvan long-read data documentation](https://github.com/GenTechGp/gtgseq/blob/main/docs/data.md) describes an ONT R10.4.1 NA12878 dataset (about 20×, roughly 11M reads/70.5 Gbases) and records public FASTQ/ENA links and checksums. It also documents minimap2-mapped derivatives. Use a deterministic downsample or chromosome subset first, then a fixed whole-genome subset. [27]
- **Tier 3: CHM13/T2T.** The [official marbl/CHM13 repository](https://github.com/marbl/CHM13) provides the T2T-CHM13v2.0 assembly, analysis-set links, public sequencing data, and CC0/public-domain reuse terms. Use `chm13v2.0_noY.fa.gz` when comparing CHM13-derived assemblies unless a Y chromosome experiment is intentional. [29]
- **Optional HG001 control.** GIAB's [HG001/NA12878 reference directory](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/) and the same Garvan documentation provide a widely used long-read control. For a genuinely separate HG002/NA24385 control, use GIAB's [Ashkenazim Trio/HG002 directory](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/). Keep HG001, HG002, GRCh38, and CHM13 results stratified; they are not interchangeable workloads. [27, 28, 29]

Record accession/release, checksum, reference build, platform/basecaller, read length distribution, coverage, and preset. Do not report candidate dataset numbers until the exact files are pinned.

### Correctness gates

Before timing, prove that stock and patched builds produce equivalent logical indexes and mappings:

1. Hand-constructed C API tests:
   - insert `base|1`, query `base`, recover the singleton bit and value;
   - insert `base`, query `base`, recover multi-hit value;
   - reject duplicate base keys;
   - verify reserve-before-insert and insert-after-freeze behavior;
   - interleave two caller-owned iterators.
2. Index-level tests:
   - compare bucket count, logical entry count, singleton count, multi-hit counts, and sorted position arrays;
   - compare direct-reference and `.mmi`-loaded mapping;
   - run ASan/UBSan and a ThreadSanitizer-oriented frozen-read test after fixing the flag race.
3. Serialization compatibility:
   - stock build writes `.mmi`, patched build loads it;
   - patched build writes `.mmi`, stock build loads it;
   - both patched/stock combinations map the same reads;
   - validate every bucket's serialized `b->n`, entry count, singleton bit, multi-hit offset/count bounds, and exact iterator count; reject short reads and trailing/early EOF rather than publishing a partial bucket;
   - permit bytewise `.mmi` differences caused by iterator order, but require logical and mapping equivalence unless byte identity is explicitly a product requirement.
4. Mapping output:
   - compare mapped/unmapped counts, primary/secondary/supplementary counts, MAPQ, reference/query coordinates, CIGAR/PAF fields, and tags;
   - canonicalize record order if hash iteration changes output order; retain raw output for investigation.

### What to measure

Separate the measurements because the hash table is only one part of minimap2:

| Phase | Measurement |
|---|---|
| Reference indexing | minimizer extraction + sorting/grouping + table insertion; wall, CPU, cycles, RSS |
| Index serialization | `.mmi` write time and size |
| Index loading | `.mmi` read/reconstruction time and RSS |
| Mapping with prebuilt index | lookup-heavy mapping, excluding index build |
| End-to-end | reference + reads to controlled sink; report separately from component timings |
| Microbenchmark | exact minimizer key/value traces, hit/miss mix, singleton/multi mix, per-bucket sizes, reserve and frozen lookup |

Use release optimization with identical compiler/flags. Run 1/2/4/8 thread matrices, at least five measured repetitions after warmups, randomized A/B order, CPU affinity, fixed governor where possible, and the same filesystem/output sink. Record:

```text
source SHA, TracEon SHA, ankerl version, compiler/version,
CPU model/ISA, OS, -k/-w/-H/-I, preset, -t,
reference/read checksums, command line, cache regime
```

Use `/usr/bin/time -v` for wall/user/system/peak RSS and `perf stat` for cycles, instructions, cache references/misses, branches, branch misses, and page faults where available. Measure cold-cache and warm-cache cases separately; do not label a run “cold” unless the cache state was controlled. Keep **index-construction threads** and **mapping threads** as separate factors: benchmark parallel index construction independently, then run prebuilt-index mapping with identical `-t` for stock and patched builds. Report end-to-end timing only alongside these decomposed measurements.

### Expected performance trade-off

No honest speedup percentage can be predicted from the container names alone.

- `ankerl::unordered_dense` is an open-addressing, cache-oriented map with robin-hood-style behavior; it may improve probe locality and lookup throughput. [32]
- khash is also an open-addressing table with quadratic probing and very low implementation overhead; it is already specialized and embedded in a C codebase.
- Minimap2 creates many **per-bucket** tables. Some buckets may be small, so per-map metadata, allocation behavior, load factor, and reserve rounding can dominate. A globally favorable generic benchmark may lose on this workload.
- The replacement's stored value pointer and freeze rule are useful for minimizing lookup overhead, but the first-read freeze transition must be made thread-safe.
- Mapping time may be dominated by chaining, alignment, I/O, or decompression. A faster table can produce a small end-to-end gain even if its lookup microbenchmark is materially faster.

The first credible result is not “X% faster”; it is a table with index-build, load, mapping, RSS, logical-entry, and correctness deltas, with confidence intervals or at least repeat distributions.

---

## Recommended integration plan

### Phase 1 — source and API bring-up

1. Fork `lh3/minimap2` at v2.31 commit `3c28777e7e2dcc90f825de1b9f17a89cca7d4452`.
2. Vendor or pin the TracEon C API implementation and `ankerl::unordered_dense`; do not build against an unpinned working tree.
3. Keep `khash_t(str)` untouched. Replace only `idxhash_t` and its per-bucket operations in `index.c`.
4. Add a small adapter block for:
   - create/destroy/reserve/insert/get/size;
   - caller-owned iteration;
   - explicit freeze;
   - fail-loud status handling.
5. Fix the TracEon shim's frozen state for concurrent reads. Prefer `std::atomic<bool>` and an explicit `freeze-all-buckets` step after index generation/load and before mapping workers begin.
6. Preserve the existing `.mmi` logical stream: dump stored key/value pairs, not container internals.
7. Add a compile-time source sentinel for the upstream `idx_hash`/`idx_eq` semantics and a full-SHA manifest.

### Phase 2 — Makefile

Use `CC` for all C translation units, `CXX` for the shim, and `CXXLD ?= $(CXX)` for final executable links. Either:

- add the shim object to `libminimap2.a`, or
- link a separately pinned `libtraceon_kmer.a` after `-lminimap2`.

Keep the stock build as the default and expose the backend through an explicit option such as `TRACEON_KMER=1`. Do not make every minimap2 user install TracEon merely to compile the default binary.

### Phase 3 — correctness

Run the hand-built API tests, all official fixtures, direct/index-first mapping, four-way stock/patched `.mmi` compatibility, ASan/UBSan, and race testing. Fail on any non-`1` insert status; never allow an incomplete index to become visible.

### Phase 4 — benchmarks

Run Tier 0/1 first, then fixed NA12878 and CHM13 subsets, then whole-genome workloads if the signal survives. Report index build, serialization, load, mapping, end-to-end, peak RSS, and microbenchmark results separately.

### Phase 5 — adoption decision

- If there is no statistically meaningful end-to-end benefit, keep the C API as an independent reusable backend and stop short of an upstream PR.
- If there is a robust benefit with no correctness, memory, or portability regression, publish `mm2-traceon` with reproducible scripts and open an upstream issue first.
- Only propose upstream integration as default behavior if it removes a measurable bottleneck for a broad workload while preserving the C-only build and supported architectures. Otherwise propose an optional backend or maintain the fork.

---

## Prioritized action list

1. **P0:** port the patch to v2.31 and fix the iterator ABI at all three call sites.
2. **P0:** discard the stale copied shim; build exactly the tracked header/implementation pair.
3. **P0:** make the C error path genuinely non-allocating/no-throw; add allocation-failure tests.
4. **P0:** fix status handling and failure-unwind behavior in C callers; never publish partial buckets.
5. **P0:** fix the mutable-`bool` frozen-read race; add explicit pre-freeze before multithreaded mapping.
6. **P1:** implement stock↔patched `.mmi` compatibility tests and logical index comparison.
7. **P1:** add a default-off Makefile backend with `CXXLD`, GCC/Clang, x86-64/ARM64 CI coverage.
8. **P1:** add exact-workload microbenchmarks and official E. coli/MT smoke tests.
9. **P2:** benchmark fixed NA12878/HG001 and CHM13/T2T subsets with reproducible provenance.
10. **P2:** compare against a minimally changed control (same v2.31, same flags, stock khash), not only against unrelated minimap2 forks.
11. **P2:** publish results and ask upstream maintainers about appetite before preparing a PR.

## Risks and caveats

- `index.c` and `khash_t(idx)` are private implementation details; upstream can refactor them without an API promise.
- v2.31 changed behavior relative to v2.30 in rare cases due to bug fixes; correctness baselines must be version-specific.
- The current experimental source is v2.28-r1209, not v2.31.
- The current local patch does not compile against the tracked iterator header; no integration performance result should be inferred from its existing object files.
- The C API's plain mutable `bool frozen` is unsafe for first concurrent lookups; fix before using minimap2's mapping threads.
- The current `thread_local std::string` diagnostic path can allocate inside `noexcept` error handling and terminate; the advertised no-throw C boundary needs implementation hardening.
- The fixed-size placement-new iterator contract needs alignment, repeated-initialization, and lifetime tests.
- `ankerl::unordered_dense` iteration order differs from khash slot order; `.mmi` byte output and alignment record order may differ even with equivalent logical content.
- Per-bucket map overhead may erase or reverse the expected lookup gain.
- The C++ runtime complicates cross-compilation, Windows, macOS, ARM, libc++/libstdc++ selection, packaging, and downstream builds.
- TracEon's C API and C++ map implementation are currently a source/static-link contract, not a stable shared-library ABI.
- Upstream acceptance is uncertain and probably low for a default dependency change: the current documented build is C/GNU-make/zlib-based, while the proposed change touches a private data structure and adds C++/dependency complexity. This is an inference about integration friction, not a claim about the maintainer's personal preference. A focused optional backend with reproducible benefit has a better chance than a default replacement.
- No upstream AI contribution policy was found; do not infer permission from silence.

## Sources

Numbered citations in the body refer to these direct sources.

1. Official minimap2 repository metadata/API: https://api.github.com/repos/lh3/minimap2
2. Minimap2 v2.31/r1302 release: https://github.com/lh3/minimap2/releases/tag/v2.31
3. v2.31 release commit: https://github.com/lh3/minimap2/commit/3c28777e7e2dcc90f825de1b9f17a89cca7d4452
4. v2.31 `index.c`: https://github.com/lh3/minimap2/blob/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/index.c
5. v2.26 `index.c`: https://github.com/lh3/minimap2/blob/e28a55be86b298708a2a67c924d665a00b8d829c/index.c
6. v2.27 `index.c`: https://github.com/lh3/minimap2/blob/b6762503a9c9af34ae02f56784035ff852ed04a5/index.c
7. v2.28 `index.c` and pinned test tree: https://github.com/lh3/minimap2/blob/8170693de39b667d11d8931d343c94a23a7690d2/index.c and https://github.com/lh3/minimap2/tree/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/test
8. v2.29 `index.c`: https://github.com/lh3/minimap2/blob/1fd85be6e2515c9194740e1d2e6a2625be36f508/index.c
9. v2.30 `index.c`: https://github.com/lh3/minimap2/blob/79c9cc186b95f50bd899f69b48eba995ced810c6/index.c
10. v2.31 README/developer guide: https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/README.md
11. v2.31 Makefile: https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/Makefile
12. v2.31 cookbook: https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/cookbook.md
13. v2.31 CI: https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/.github/workflows/ci.yaml
14. v2.31 MIT license: https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/LICENSE.txt
15. VGP assembly pipeline: https://github.com/VGP/vgp-assembly/tree/master/pipeline/purge_dups
16. VGP asset/version documentation: https://github.com/VGP/vgp-assembly/tree/master/pipeline/asset
17. nf-core minimap2 module: https://github.com/nf-core/modules/tree/master/modules/nf-core/minimap2/align
18. mm2-fast fork README: https://github.com/bwa-mem2/mm2-fast
19. mm2-fast `index.c`: https://github.com/bwa-mem2/mm2-fast/blob/master/index.c
20. minimap2-pure-rs optimization notes: https://github.com/henriksson-lab/minimap2-pure-rs/blob/main/OPTIMIZATION.md
21. minimap2 index modifier: https://github.com/ispras/minimap2_index_modifier
22. Generic C/C++ hash benchmark and report: https://github.com/JacksonAllan/c_cpp_hash_tables_benchmark/ and https://jacksonallan.github.io/c_cpp_hash_tables_benchmark/
23. GCC C++ driver documentation: https://gcc.gnu.org/onlinedocs/gcc/Invoking-G_002b_002b.html
24. Mixed C/C++ static-library example: https://github.com/noahp/c-cpp-lib-example
25. TracEon C API header at observed commit: https://github.com/woosflex/TracEon/blob/9566b24dfee04e7d582fa082744cf37c8d59b99b/include/kmerindex_c_api.h
26. TracEon C API implementation at observed commit: https://github.com/woosflex/TracEon/blob/9566b24dfee04e7d582fa082744cf37c8d59b99b/src/kmerindex_c_api.cpp
27. Garvan long-read data documentation: https://github.com/GenTechGp/gtgseq/blob/main/docs/data.md
28. GIAB HG001/NA12878 reference directory: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/
29. T2T CHM13 data/reference repository: https://github.com/marbl/CHM13
30. rammap performance comparison: https://github.com/jwanglab/rammap/blob/main/docs/performance.md
31. mm2-fast paper, Kalikar et al. 2022: https://doi.org/10.1038/s43588-022-00201-8
32. ankerl unordered_dense: https://github.com/martinus/unordered_dense
33. Local verification log (commands, compiler, outputs): `outputs/.plans/minimap2-traceon-verification.log`
34. Research provenance sidecar: `outputs/minimap2-traceon-integration-research.provenance.md`
35. GIAB HG002/NA24385 directory: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/
36. Minimap2-on-GPU derivative: https://github.com/Minimap2onGPU/mm2-gb
