# TracEon direction audit — earliest prototype through current development

Reviewed 2026-09-22. Read-only inspection of all locally available branch histories,
release tags, historical implementations, benchmarks and the current working tree.
No checkout, reset, code change or benchmark rerun of the historical C++ binaries.
One minimal Python probe reproduced a historical argument-parsing defect.

**Finding: there has been a substantial shift in product emphasis and in the
question being benchmarked. The original local-library proposition has not been
fairly settled by the recent remote-service comparisons.** This does not prove
that returning to it will produce a superior tool or that the old performance
claims were sound.

## What the project originally set out to do

The first release README, commit `499d7f8` (27 September 2025), describes an
in-memory genomic cache intended to eliminate repetitive file parsing. Its API
loads FASTA, retrieves by ID, and saves/restores prepared data. The initial
implementation has plaintext and two-bit strategies; the benchmark uses the
constructor's default plaintext strategy. Python/R bindings and concurrent
access are future accessibility goals. There is no mapper-engine replacement
or network service in this prototype.

The v0.0.1 documentation keeps the same purpose, adding FASTQ, gzip, parallel
parsing and deliberately **uncompressed binary persistence prioritizing restore
speed**. Its planned next step is a general C ABI followed by Python/R bindings.

By v1.0.0, the main design is an embedded C++ library: an immutable arena,
borrowed `string_view` results, concurrent reads and a memory-mapped binary
artifact. Initial loading can consume gzip, but lookups access already decoded
records. The README discusses interactive exploration and integration into
local applications. Its roadmap also speculates about compression, distributed
caching and remote access, so later expansion was already contemplated; it was
not invented entirely during our recent conversation.

## Where the direction changed

| Stage | Evidence | Change and implication |
|---|---|---|
| September–October 2025 prototype | `499d7f8`, `388f3b1`, v0.0.1 | Avoid repeated parsing; cache records locally; save/restore across runs. |
| December 2025 v1 | `95afc06`, `3bf2f9e`, v1.0.0 | Borrowed views and mapped payloads make local access cheap; Python/R accessibility remains planned. |
| June 2026 compression | `49e13ed`, later `883a81f` / v1.3.0 | New cache writes become compressed. Restore allocates a decoded arena and decompresses into it. Smaller disk artifacts trade away direct access to mapped payload bytes. v1.3 still accepts old uncompressed files. |
| August 2026 k-mer expansion | `bca8a79` / v2.0.0 | Adds an optional k-mer/minimizer C API, a different problem from record caching. It is not the originally planned general FASTA/FASTQ binding API. |
| August 2026 remote access | `baaf247`, `e0275c9` / v2.3.0 | TCP server/client introduces independent-process service access and owned transferred records. Remote access becomes possible but incurs serialization, copies and transport. |
| Current uncommitted publication work | current README, ROADMAP, manuscripts and studies | Shared service becomes the primary narrative; mapper replacement and competitive service performance dominate evaluation. These changes are not historical release commits. |

Safety work—parser repairs, record integrity, explicit lifetimes, bounds checks,
atomic persistence and tests—strengthens the original product. It should not be
confused with strategic drift or discarded in a rollback.

## The most consequential technical change: mapped payloads versus decoded heaps

In v1.0.0 `SmartStrategy::loadBinary`, the implementation maps the uncompressed
file, walks serialized records, rebuilds a hash table, and constructs sequence
views pointing directly into the mapped bytes. The mapping stays alive.

Under read-only access, the OS can share those file-backed payload pages among
independently started processes. **That is an architectural property, not a new
measurement in this audit.** Each process still reconstructs and owns its hash
index; v1 is not a fully shared serialized index and restore is not O(1).
Mapping also does not eliminate page faults or the cost of eventually reading
the sequences. The README's “instant restore” wording overstates what the code
guarantees.

Commit `49e13ed` introduces the compressed restore path: allocate
`text_arena_`, decompress with LZ4, release the mapping, then index the decoded
bytes. Subsequent formats retain a decoded in-process working set. Independent
restores therefore do not get the same file-backed sharing of decoded payloads.
A single server avoids separate decoded heaps, but serving bytes adds IPC costs.
Compression is a legitimate storage trade-off; the history does not establish
that it was a performance mistake for every workload.

The current local `Cache::getView` API still exists. Thus the embedded-library
idea has not disappeared. What changed is the persistence model and which access
path we emphasize and evaluate.

## Why the prototype and recent results look different

The original local operation can return a pointer and length. The recent
service operation delivers owned sequence/quality bytes across processes and
validates them with a common consumer. Long records increase transfer cost in
the service but do not make a borrowed-view lookup proportionally more expensive.

Therefore a remote long-record loss is **not a measurement of local borrowed-view
lookup performance**. Equally, a very fast view lookup does not establish that
an application can scan, transform or align the whole record equally quickly.
Both lookup-only and realistic payload-consuming workloads are needed.

Recent comparisons correctly challenge universal speed, memory and mapper-engine
claims. They do not alone settle whether an embedded TracEon library improves a
pipeline currently reparsing data or allocating large per-record objects. My
previous interpretation put too much weight on service results when discussing
the whole project's usefulness.

## The old evidence is not a reliable gold standard

1. **First 64,000× headline:** `499d7f8:benchmark_disk.cpp` scans the entire FASTA
   file for each request, even after the requested record has ended. The cache
   benchmark times one already-loaded lookup separately from setup. This tests
   avoiding repeated scans; it does not establish superiority over indexed
   random access. It also lacks an exact-output comparison between arms.
2. **v0.0.1 lookup coverage:** `gcp_benchmark.cpp` constructs `chr`-prefixed IDs
   instead of sampling actual input names, and does not check each hit. On
   inputs with other naming conventions, misses can contribute to throughput.
   Its restore “integrity” check compares record counts, not all payloads.
3. **v1 lookup semantics:** `traceon_driver.cpp` calls `getView()` and consumes
   only `.size()`. Its reported “Data_Retrieved_MB” is summed view lengths, not
   bytes actually read or copied by that loop. Long-read lookup figures must
   not be interpreted as whole-record processing bandwidth.
4. **Broken tagged v1 competitor invocation:** `benchmark_runner.py` passes
   `file, iterations, prefix, count`, but the pyfastx script reads iterations
   from `argv[3]` (the prefix) and count from nonexistent `argv[5]`. The minimal
   reproduction gives `ValueError: invalid literal for int() with base 10: 'r'`.
   The outer reporter still divides requested iterations by elapsed subprocess
   time without rejecting the nonzero return code. That tagged script cannot
   substantiate a valid pyfastx speedup as written. This does not establish which
   script/version produced every historical README number.
5. Even with its arguments repaired, the historical comparison times Python
   interpreter/index initialization and key generation against a C++ inner loop
   whose setup and precomputed keys are excluded. The BioPython parser stores
   a list of rich objects; that is not an exhaustive comparison of BioPython's
   possible access modes or a universal memory baseline.

The relevant distinction is: **a plausible original problem and useful design,
but incomplete and sometimes defective evidence for the scale of its advantage.**

## Recommended correction in direction

Reassess the original proposition explicitly: a reusable local genomic record
store for applications repeatedly accessing the same dataset, with borrowed
views and predictable lifecycle costs. Treat service access as an optional
integration path, not as the sole definition of TracEon.

The next discriminating experiment should compare the current local API with
appropriate indexed readers under identical semantics: lookup-only separately
from reading every sequence/quality byte, one-time preparation and restore,
multiple independent processes, and total memory. Then isolate the historical
storage trade-off with a bounded comparison of safe uncompressed mapped payloads
versus compressed decoded arenas. This is a hypothesis to test, not authorization
to restore the unsafe v1 loader or build a new architecture immediately.

The original Python/R accessibility goal also remains largely unfulfilled in
the inspected tree: there is no shipped general record-access binding in the
tracked/current source files examined. The k-mer C API serves a different use
case. That gap between a fast C++ component and access from researchers' existing
workflows may matter more for adoption than another synthetic throughput record.

No roadmap or production code was changed by this audit. The result supports
re-centering the evaluation, while preserving the modern correctness work and
the counterexamples already measured.

## Evidence and reproduction

Historical source snapshots and the argument-parsing probe are saved locally in
`build/history-direction-audit/`. They are extracted without checking out old
code. Tag objects were resolved to commits: v0.0.1 → `6f3b30c`, v1.0.0 →
`600734f`, v1.3.0 → `883a81f`, v2.0.0 → `98ea2c2`, v2.3.0 → `e0275c9`.
Dates in the narrative follow the relevant implementation commits rather than
assuming a tag's version name identifies its earliest implementation.

Key read-only commands:

```sh
git show 499d7f8:README.md
git show 499d7f8:benchmark_disk.cpp
git show v0.0.1:README.md
git show v0.0.1:gcp_benchmark.cpp
git show v1.0.0:src/SmartStrategy.cpp
git show v1.0.0:traceon_driver.cpp
git show v1.0.0:benchmarks/benchmark_runner.py
git show 49e13ed:src/SmartStrategy.cpp
git show bca8a79 --stat
git show baaf247 --stat
```

Subsequent evidence: the [earliest-claim rerun](earliest-claims-rerun.md) executes
the same-day repaired first prototype. It confirms a local reused-cache benefit
but does not reproduce the original 64,000× magnitude or three-request payback.

The [v1 mapped-cache follow-up](mapped-cache-study.md) now measures the sharing
property described above. On its 100 MB long-record input, four independent v1
workers use 105.75 MiB aggregate PSS versus 390.94 MiB for current local restores.
Both local implementations beat the persistent HTSlib reader on the tested
repeated sequence workloads; HTSlib remains smaller in process PSS. This is
single-host synthetic evidence for the local-cache direction, not mapper
replacement or proof that all historical benchmark claims hold.
