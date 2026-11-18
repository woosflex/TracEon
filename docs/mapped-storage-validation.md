# Restoring TracEon's local-cache direction

22 September 2026, uncommitted 3.0.0-dev candidate.

**The current library now supports validated mapped persistence.** Prepare a
cache once, open it in independent processes, and query locally. Payload pages
are shared by the OS; each process owns its hash index. Compression and the
optional TCP service remain supported.

## Implemented path

```sh
traceon_make_cache reads.fastq reads.traceon --mapped
# With remote components enabled, this equivalent builder is also available:
traceon build reads.fastq reads.traceon --mapped
```

```cpp
TracEon::Cache cache;
cache.restore("reads.traceon"); // auto-detect storage; validate before publication
if (cache.hasSequence("read1")) {
    auto sequence = cache.getView("read1"); // borrowed, local sequence bytes
    auto record = cache.getFastqRecord("read1"); // owned sequence + quality
}
```

C++ writers select `cache.save(path, TracEon::CacheStorage::Mapped)`.
One-argument `save(path)` and builders without `--mapped` retain LZ4 behavior.
The added v5 codec 0 stores the existing portable payload directly. Both codecs
use the same CRC32C, checked record parser, duplicate rejection, and atomic
writer. Earlier development v5 readers reject codec 0; released v2.3.0 readers
support neither v5 codec.

Restore checks every payload byte and reconstructs the private index. This is
not the old unchecked reader or an instantaneous restore. Mapped files must
remain immutable: never modify or truncate an open file in place. POSIX atomic
replacement preserves existing mappings while new readers open the replacement.
Borrowed views still expire on cache clear, reload, or destruction.

## Controlled comparison

Unlike the [v1/current historical comparison](mapped-cache-study.md), both arms
here use **the same updated library and worker executable**. Only storage differs.
Two deterministic FASTA inputs contain approximately 100 decimal MB of sequence
payload: 666,666 records of 150 bases, and 3,333 records of 30,000 bases. They reuse
the historical study's inputs and independent reference SHA-256 digests.

Six repetitions per input; reproducibly shuffled arm order; query worker on CPU 8;
independently exec'd memory workers on CPUs 8–11; one Intel Core Ultra 5 125H host.
After the same 10,000 untimed query IDs, time 10,000 complete sequence queries
with common length-framed ID/sequence SHA-256 consumption. Restore is measured
separately. For memory, every worker verifies every record before waiting for
PSS measurement. No cache is inherited through fork.

All **40,439,940 complete record deliveries** matched their reference digests:
24 query observations and 48 memory cohorts containing 120 independent workers.

Medians of six observations:

| Metric | 150-base records: mapped / LZ4 | 30-kb records: mapped / LZ4 |
|---|---:|---:|
| Restore, ms | 198.92 / 318.82 | 27.27 / 138.82 |
| 10,000 consuming queries, ms | 11.54 / 12.87 | 170.83 / 163.02 |
| One-reader PSS, MiB | 160.35 / 160.54 | 98.05 / 98.30 |
| Four-reader aggregate PSS, MiB | 311.71 / 639.78 | 104.07 / 390.67 |

Mapped restore was faster in all six pairs for both inputs. Exploratory 95%
bootstrap intervals of paired mean restore savings were 102.21–148.26 ms for
short records and 90.68–111.14 ms for long records. Query-time difference intervals
crossed zero for both inputs: **no resolved query-speed advantage is claimed**.
The implementation recovers shared payload memory and faster preparation/restore,
while retaining roughly comparable local query performance in these tests.

Cache builds include process launch, text parse, serialization, and atomic save.
They were measured once per input/arm, so these are descriptive observations:

| Input | Mapped build, s | LZ4 build, s | Mapped file, MB | LZ4 file, MB |
|---|---:|---:|---:|---:|
| 150-base records | 0.706 | 43.511 | 114.56 | 43.39 |
| 30-kb records | 0.408 | 47.752 | 100.06 | 39.26 |

The existing LZ4 writer selects high compression for these large nucleotide
inputs. Mapped mode avoids that compression work at the cost of larger artifacts;
this does not establish that compression should be removed or that these build
ratios generalize.

## Correctness and limitations

All 23 CTest checks passed, including new mapped-format/lifecycle tests, injected
write/sync/rename failures for both storage modes, installed-consumer access,
and complete FASTQ service semantics over a mapped backing file. Directed tests
reject every truncated prefix and every single-byte mutation of a small mapped
cache, plus CRC-valid duplicate IDs, oversized lengths/counts, and trailing data.
Both index modes, empty records, compressed/mapped conversions, and readers held
across atomic replacement are covered. Bounded save allocation is checked for
both modes. An AddressSanitizer/UndefinedBehaviorSanitizer/libFuzzer run completed
285,862 iterations in 31 seconds with exit status 0 and no findings. The first
sandboxed run completed 320,586 iterations but LeakSanitizer could not finalize
under ptrace; the successful repeat ran outside that sandbox with leak detection
enabled. Both logs are retained.

This remains single-host, warm-filesystem-cache, synthetic sequence-only timing
evidence. FASTQ correctness is tested, but FASTQ performance is not measured here.
The host is not exclusive and frequency is uncontrolled. PSS excludes the
controller and unmapped filesystem page-cache memory; it is not total-system
memory. No new HTSlib/pyfastx comparison or real research workflow is claimed.
The previous indexed-reader and service losses remain valid for their tested paths.
Windows mapped replacement and independent-host replication are untested.

The [roadmap](../ROADMAP.md) now prioritizes an actual repeated-record workflow,
then a narrow binding if its users need one. Mapper-engine replacement and
service expansion are deferred. This change restores a demonstrated architectural
benefit; it does not establish novelty or publication readiness.

## Reproduction

```sh
python3 benchmarks/mapped_storage_study.py \
  /tmp/traceon-mapped-check-20260922 build/NEW-mapped-storage-run
```

The runner requires the sealed inputs and equally warmed worker from
`build/mapped-cache-2026-09-22`, plus local HTSlib 1.24 and the specified current
build. Output must be a new directory. It records commands, build/core hashes,
source snapshots, expected digests, individual observations, paired analysis,
and file hashes. Full local evidence: `build/mapped-storage-2026-09-22`.

Compact [evidence archive](results/mapped-storage-2026-09-22.tar.gz) and
[SHA-256](results/mapped-storage-2026-09-22.tar.gz.sha256) preserve scripts,
source, raw observations, and validation logs. Regenerable caches and executable
binaries are omitted; their hashes are retained in the manifest. Large input
hashes and reference digests are in `inputs.json`.
