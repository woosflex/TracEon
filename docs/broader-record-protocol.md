# Broader record-access validation protocol — 2026-09-22

Prospective extension of the persistent-reader experiment. No universal speed,
memory, novelty, mapper replacement or cross-host claim is presumed true.

## Comparators and inputs

One C++ driver uses native HTSlib 1.24, TracEon, or embedded CPython with pyfastx
2.3.1's public subscription/sequence/quality API. Common C++ SHA consumption,
result destruction and barriers preserve the existing metric. Python interpreter
startup/PSS belongs only to pyfastx. This measures backend/API integration,
not ordinary Python application throughput or language-independent algorithms.
Fasta preserves case; only KeyError represents a miss. Smoke checks FASTA/FASTQ,
plain/gzip, mixed case, punctuation, duplicates and partial batches. A separate
missing-ID probe is retained: preliminary pyfastx2.3.1 FASTQ checks found that
a miss poisons subsequent valid lookups on the same handle, for plain and gzip.
Performance fixtures contain only present IDs; no failed lookup is scored as fast.

Inputs: Drosophila 100k FASTQ, synthetic 1m FASTQ, corrected E. coli long FASTA,
from the earlier hash-recorded experiment. Gzip level6/mtime0 outputs must
round-trip byte-for-byte. Three fresh pyfastx indexes per representation and
three cache builds from gzip are measured. The original plain cache/HTSlib
indexes are reused with hashes checked. Compressed cache builds must expose the
same snapshot and requested records. Gzip generation is reported separately.

## Warm local comparison

Four arms: HTSlib plain, pyfastx plain, pyfastx ordinary gzip, TracEon (decoded
snapshot, identical irrespective of original compression). Six sessions per
input; server restarted each session. Cases: one worker batch64 random, four
workers scalar random, four workers batch64 random (primary), four workers
batch64 source-ordered. Three arm permutations and their reverses provide equal
pairwise precedence across six sessions (not equal position counts); case order shuffled with fixed seed.
Each worker requests 1,024 records (512 long), prefixes of independently
seeded prior session random ID lists; ordered cases reorder the same selection
by source position before the execution seal; warmup/timing/memory passes reuse that list.
Primary time includes common digest consumption. Memory measured separately.
All outputs must match expected framed digests. Every outcome is retained.
Report six session-paired differences and exploratory bootstrap intervals,
not individual requests as independent replicates. Gzip decompression versus
predecoded residency must not be described as equivalent storage behavior.

## Isolated memory limits

Synthetic 1m, four workers, batch64 random, same four arms, six sessions;
limits 256/768/1536 MiB and zero swap, in transient systemd user services.
Controller, server and workers share each cgroup; worker/server CPUs8–11,
controllerCPU16. No host-wide cache eviction or pressure allocation.
Create fresh source/index/cache copies by ordinary read/write inside each
cgroup so their newly instantiated file pages are charged to that group.
Do not benchmark externally warmed original files inside a limit. File copying
is preparation, outside measured retrieval, with cgroup counters before/after.
Record memory.current/peak/stat/events/pressure/swap, failures, OOM and timeouts.
Preserve failed arms as failures, never treat missing timings as fast results.
Three passes match main study. Compare successful paired sessions only while
reporting completion denominator; do not infer universal minimum memory.

## Scope and provenance

Same local host, nonexclusive CPUs, warm-cache main study; limit tests involve
real reclaim where counters show it, not a cold-cache/storage benchmark.
Independent-host reproduction requires a second machine and remains separate.
Freeze source/helper/protocol/tool/input/request hashes after smoke and before
measurements; record environmental limits, errors, all raw data and analysis.
