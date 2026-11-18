# Rerunning the earliest reported TracEon numbers — 2026-09-22

**The original cache-versus-rescan advantage is reproduced; the published
64,000× magnitude and three-request payback are not reproduced on this host.**
A stronger native indexed-reader control also finds a local advantage after
sufficient reuse. This is evidence for the original embedded-cache use case,
not a universal advantage or an exact reconstruction of historical hardware.

## What was reconstructed

The first numerical claim is in commit 499d7f8 (27 September 2025): generated
FASTA with 50,000 sequences of 1,000 bases;83.83 ms scan,174.47 ms cache load,
0.0013 ms cached lookup, and cost recovered after 3 requests. No exact input file,
generator seed, compiler flags or complete machine specification accompanies it.

The exact commit fails to compile on the current compiler because of missing
integer-type includes. With a forced cstdint include, it still fails to link:
PlainTextStrategy.cpp contains declarations rather than its method definitions.
Both failures are preserved. Runnable measurements use same-day historical
follow-up fa95f1d. **Cache.cpp and both original benchmark drivers are byte-
identical between these commits**; the follow-up restores the implementation
and fixes destructor linkage. A forced cstdint include remains necessary for
modern compiler compatibility. No current TracEon cache implementation is used.

Two builds cover plausible missing original flags: C++20/O0 and C++20/O3/NDEBUG.
Inputs use Python Random seed 20260922,50,000 unique IDs and 1,000 random ACGT bases
per record. Single-line and 80-column-wrapped files encode the same records.
The source generator, source hashes, compile commands, binaries' hashes, host
metadata and [prospective protocol](earliest-claims-protocol.md) are preserved.
Executables run on CPU 8 with warm OS file caches; the CPU is not exclusively
reserved. Historical “disk” access is a full-file scan, not necessarily physical
storage I/O in these warm runs. No system cache eviction was performed.

## Original unmodified benchmark drivers

There are 72 paired trials: two builds × two layouts × first/middle/last target
IDs × six repetitions. Each target/build/layout has three scan-first and three
cache-first repetitions. Cache load and one lookup are timed separately exactly
as the historical driver does. Table entries are medians across 18 observations
per build/layout. Ratios are medians of paired scan/lookup ratios.

| Configuration | Full-file scan ms | Cache construction ms | Single cached lookup µs | Paired scan/lookup ratio |
|---|---:|---:|---:|---:|
| Original README claim |83.83|174.47|1.30|over 64,000×|
| O0, single-line |24.89|130.07|8.60|2,884×|
| O0,80-column wrapped |39.34|161.02|8.70|4,465×|
| O3, single-line |13.81|81.36|4.82|2,842×|
| O3,80-column wrapped |23.01|95.36|4.92|4,636×|

Single-lookup microsecond measurements are sensitive to timer overhead, first
access and allocation effects. All ranges are retained in summary.json; the
observed O3 single-line lookup range is 3.975–5.958 µs. These measurements do not
prove the old numbers could never have occurred under different conditions.
They do show that the exact published magnitude is not reproducible from the
available specification on this machine.

The scan traverses the entire file even after the target has been found. Its
large disadvantage is consequently evidence against repeated full scans, not
against established random-access indexes.

## Was the claimed three-request payback real here?

A separate driver consumes all returned bytes through an independently checked,
length-framed SHA256 digest. Each cache arm includes fresh parsing/construction;
the other arm performs the original full scan once per request. Six pairs per
request count are balanced 3/3. Total is setup plus verified retrieval; process
launch, reading the query-ID file, and digest initialization are excluded.

| O3 single-line requests | Repeated scans median ms | Cache setup + queries median ms | Cache wins / pairs |
|---|---:|---:|---:|
|1|13.86|79.28|0/6|
|2|27.44|80.55|0/6|
|3|40.91|80.76|0/6|
|4|54.82|81.00|0/6|
|8|110.97|82.71|6/6|
|16|218.93|82.40|6/6|

The O0 experiment has the same win/loss pattern. **Three requests do not pay
back setup here. Eight is the first tested winning count**, not a determination
that the exact crossover is eight; counts 5–7 were not tested. The original
README's three-request conclusion follows arithmetically from its stated inputs,
but those inputs were not reproduced.

## Does the prototype beat a real indexed reader?

Yes, in the measured local repeated-access workload, after setup is amortized.
The stronger control uses existing HTSlib 1.24, an already-built FAI, and one
persistent handle per measured query sequence. The same historical Cache and
common byte-consuming SHA256 verification are used. This is a separate modern
control, not a reproduction of an indexed comparison from the first README.

The initial 64-query control materializes HTSlib output into std::string to match
Cache::get. It measures 0.128 ms cache versus0.946 ms indexed retrieval, excluding
setup. Because the extra indexed string copy is avoidable, a subsequent native
buffer control consumes HTSlib's allocated buffer directly. That stronger control
uses six balanced pairs and verifies identical outputs:

| Requests | Cache queries ms | Native HTSlib queries ms | Cache setup + queries ms | Native HTSlib setup + queries ms | Cache total wins |
|---|---:|---:|---:|---:|---:|
|64|0.135|1.021|88.04|26.70|0/6|
|10,000|15.16|134.84|93.62|153.64|6/6|

Thus the cache is faster in active retrieval, but indexed access is preferable
for the small total workload. At10,000 requests, the cache wins even including
its full parse/construction, by about 1.64× in median total time. At64 requests,
its faster lookups cannot recover setup. These ratios describe the tested
50 MB corpus and this host, not general library superiority.

An exploratory reuse sweep with the string adapter gives 0/6 cache wins at 1,000
requests,4/6 at 5,000, and 6/6 at 10,000 and 50,000. The sweep and native-buffer
control were selected after examining the initial control and are marked as
exploratory follow-ups. Their plans, source variants and all results are saved.
FAI creation took 0.230s separately; it is deliberately not charged to each
indexed arm. The cache starts from the FASTA each time. No universal amortization
threshold or statistical significance claim is made from these six pairs.

## Correctness and scope

Small plain/wrapped mixed-case/ambiguous-base fixtures pass 8 smoke cases across
both builds and four access/roundtrip paths. All 50,000 source records are checked
for cache load and save/restore under both builds and both layouts:400,000 full
record verifications. Across those checks, initial controls, scan-amortization
trials and indexed follow-ups, **1,314,736 returned record instances match the
independent expected digests** (excluding the small smoke checks). Original
single-lookup drivers do not emit records; their correctness is supported by
the separately verified unchanged parser/retrieval/scanner implementations,
not by pretending their timing output contains a digest.

This evaluates the first September 2025 plaintext Cache, not the v1.0 mapped
implementation or v0.0.1's multi-gigabyte GCP claims. Those other historical
figures remain unverified. The first benchmark uses the default plaintext
strategy; it does not validate the two-bit strategy's advertised compression.
No memory-efficiency claim is tested here. Existing-index retrieval, source
layout, reuse level, machine and compiler all affect the practical outcome.

**Interpretation:** the earliest project was addressing a real computational
cost—repeated parsing—and this reconstructed prototype additionally demonstrates
an advantage over native indexed access for a sufficiently reused, memory-resident
local workload. Its original headline overstated what its benchmark established.
This supports evaluating the embedded-library direction on real workflows;
it does not undo the remote-service counterexamples or establish novelty/adoption.

## Artifacts

Full local artifacts: `build/earliest-claims-2026-09-22/`.
[Compact archive](results/earliest-claims-2026-09-22.tar.gz) and
[SHA256 checksum](results/earliest-claims-2026-09-22.tar.gz.sha256) retain historical
source, harnesses, compile failures, commands, seals, raw outputs and summaries.
Large generated inputs/cache files and executables are omitted; deterministic
regeneration and hashes are recorded. The archive includes the generated native
buffer source and its controller so that the last control is reviewable.

Runnable harnesses:

```sh
python benchmarks/earliest_claim_study.py prepare NEW_OUTPUT
python benchmarks/earliest_claim_study.py execute NEW_OUTPUT
python benchmarks/earliest_claim_study.py summarize NEW_OUTPUT
python benchmarks/earliest_indexed_followup.py NEW_OUTPUT
```

Before prepare, extract 499d7f8 and fa95f1d into matching subdirectories of the new
output, as in the archive. Local HTSlib/tool paths are explicit in the harness;
adapt and record them on another machine. The archived native_control.py runs
its supplementary native-buffer check. Sources must remain unchanged after
sealing; use fresh output directories for repetitions.
