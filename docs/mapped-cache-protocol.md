# Local mapped-cache historical validation

This bounded follow-up tests v1.0.0 (`600734f`, resolved full revision in the
run metadata), the current uncommitted library, and HTSlib 1.24 on this host.
It follows the earliest-prototype rerun. This protocol was documented after
execution began; it is not a preregistration.

- Generate deterministic, single-line ACGT FASTA with 150-base and 30,000-base
  records, at approximately 10 and 100 decimal MB of sequence payload each.
  Python Random seed 20260922; independent SHA-256 expectations derived during
  generation. These are new synthetic inputs, not the missing original inputs.
- Build each cache and FAI once. Record parsing/indexing and save times separately;
  these single observations are descriptive, not repeatable build-time estimates.
- Compile the same worker at O3 against the historical and current libraries.
  Preserve v1 source, including its Robin Hood map header. Force-include cstdint
  for compiler compatibility; use current system zlib. No historical algorithm edits.
- Six repetitions per input, with reproducibly shuffled arm order. Single query
  worker pinned to CPU 8. Builds restricted to CPUs 8–11. No other benchmark runs
  concurrently; the workstation is not otherwise exclusive.
- Precompute 200,000 random IDs with replacement. For the two cache libraries,
  warm those lookups and time five passes (one million lookups). Consume lengths
  only, checking every length. This tests borrowed-view access, not payload bandwidth.
  Do not report HTSlib fetches as equivalent zero-copy lookup operations.
- In each query worker, separately time 10,000 complete sequences using the same
  length-framed ID/sequence SHA-256 operation. HTSlib's native returned buffer is
  hashed directly, with no additional sequence-string copy. Cache workers have
  just completed lookup warmup; the indexed arm has no equivalent prefetch pass.
  These are warm-host observations, not a controlled cold-cache comparison.
- Time each fresh restore/index open separately. OS file caches are not flushed.
  Query times exclude process launch, cache construction, restore, and ID generation.
- For resident memory, start one or four fresh independent executable processes,
  with one CPU each (8–11). Each verifies every sequence in source ID order and
  then waits. Sample Linux smaps_rollup RSS/PSS, including anonymous/file PSS;
  sum across workers. Repeat six times. No cache is inherited through fork.
  Parent/controller and unmapped OS page-cache memory are excluded, so PSS is
  resident process attribution, not whole-system or cgroup memory consumption.
- Stop on build failure, wrong record count, wrong length, or wrong digest.
  Save commands, raw observations, environment, hashes, and summary medians.

Scope: sequence-only FASTA, two sizes, one host. No claims about FASTQ quality,
wrapped FASTA, gzip access, safety of historical formats, pyfastx/BioPython ratios,
or mapping/alignment throughput follow from this experiment. Restoring v1's
unchecked binary reader to production is not part of this study.

Run from the repository root with the existing local dependencies:

```sh
python3 benchmarks/mapped_cache_study.py
```

The script uses the portable current build in
`/tmp/traceon-release-check-20260922` and HTSlib under
`build/research-inputs/samtools-1.24/htslib-1.24`. It exports historical files
using `git show`; it does not switch branches or modify production source.

Two exploratory follow-ups were added after inspecting the primary results:

1. Compile v1's original `traceon_driver.cpp` unchanged, including its original
   key generation and volatile length accumulator. Run one million lookups six
   times per cache; require zero misses and the expected aggregate length.
   This is the closest historical lookup rerun, still on newly generated inputs.
2. Compile a minimal worker variant that skips the million cache-only lookups
   and first issues the same 10,000 untimed query IDs to every arm. Then time
   the common SHA-256-consuming loop. This removes the unequal query warmup;
   cache warmup obtains views, while HTSlib warmup fetches native buffers.
   These controlled follow-up timings are the report's primary payload comparison.

The exact follow-up scripts and generated variant source are retained in the
artifact directory and evidence archive. Six repetitions per input/arm are
retained without dropping outliers. Exploratory bootstrap intervals describe
paired mean time differences within this single host; they are not evidence
of generality across machines or datasets.
