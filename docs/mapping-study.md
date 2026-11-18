# TracEon mapping integration study

Status: feasibility study, not submission-ready evidence.
The initial per-process-restore experiment below is a control; it does not test
TracEon's intended cache-once, query-from-many-processes deployment.
For that model, start with the [corrected shared-server experiment](#corrected-experiment-cache-once-query-from-many-processes).
Started 2026-09-21. The research sequence is: establish a suitable integration,
validate equivalent outputs, and measure complete workflows. Do not infer a speedup from cache
lookup microbenchmarks.

## First candidate: minimap2

Pinned upstream: [v2.30, 79c9cc186b95f50bd899f69b48eba995ced810c6](https://github.com/lh3/minimap2/tree/79c9cc186b95f50bd899f69b48eba995ced810c6).

Source inspection identifies three distinct operations:

- `mm_idx_get` queries a minimizer-to-position index. TracEon's sequence-ID index
  cannot replace it while preserving the mapping algorithm.
- `mm_idx_getseq` and `mm_idx_getseq_rev` retrieve numeric-ID intervals from
  packed reference bases already in RAM. Any TracEon adapter must preserve
  interval bounds, numeric base encoding, ambiguity handling, and reverse
  complements. ASCII views are not interchangeable with these buffers.
- `mm_bseq_read3` streams query sequences in batches. Replacing that stream
  with an in-memory ID cache is not automatically helpful. Repeated subset
  retrieval is a separate workload that needs evidence from a real workflow.

Sources: pinned [index.c](https://github.com/lh3/minimap2/blob/79c9cc186b95f50bd899f69b48eba995ced810c6/index.c),
[bseq.c](https://github.com/lh3/minimap2/blob/79c9cc186b95f50bd899f69b48eba995ced810c6/bseq.c),
and [alignment callers](https://github.com/lh3/minimap2/blob/79c9cc186b95f50bd899f69b48eba995ced810c6/align.c).

Proceed with an internal reference adapter only if profiling shows a meaningful
retrieval cost and there is a credible way to improve it without excessive memory
or conversion overhead. If not, select a demonstrated repeated-read retrieval
workflow rather than forcing a replacement into minimap2.

## Baseline reproduction

Requires Linux, Python 3, git, a C compiler, make, and zlib development files.
The runner uses only the Python standard library. It neither installs tools nor
downloads data. Use a new output directory for every run.

```sh
git clone --branch v2.30 --depth 1 https://github.com/lh3/minimap2.git /tmp/traceon-minimap2-pilot
git -C /tmp/traceon-minimap2-pilot rev-parse HEAD
# Must print 79c9cc186b95f50bd899f69b48eba995ced810c6.
make -C /tmp/traceon-minimap2-pilot -j4
curl -fL --retry 2 https://github.com/lh3/minimap2/releases/download/v2.10/cookbook-data.tgz -o /tmp/traceon-minimap2-cookbook.tgz
mkdir -p /tmp/traceon-mapping-data
tar -xzf /tmp/traceon-minimap2-cookbook.tgz -C /tmp/traceon-mapping-data ecoli_ref.fa ecoli_p6_25x_canu.fa
python benchmarks/mapping_baseline.py \
  --minimap2 /tmp/traceon-minimap2-pilot/minimap2 \
  --reference /tmp/traceon-mapping-data/ecoli_ref.fa \
  --reads /tmp/traceon-mapping-data/ecoli_p6_25x_canu.fa \
  --output build/mapping-pilot-2026-09-21 \
  --preset map-pb --threads 1 --repetitions 5 --expected-reads 13124 \
  --source-revision 79c9cc186b95f50bd899f69b48eba995ced810c6 \
  --dataset-source https://github.com/lh3/minimap2/releases/download/v2.10/cookbook-data.tgz
python tests/test_mapping_baseline.py
```

Dataset choice follows the upstream [PacBio mapping example](https://github.com/lh3/minimap2/blob/79c9cc186b95f50bd899f69b48eba995ced810c6/cookbook.md#map-pb).
This is an example dataset, not a representative multi-organism evaluation.

Outputs include all commands, exit codes, binary/input hashes, CPU and memory
information, index construction, one warmup, five fresh-process mapping samples,
per-child peak RSS, and every SAM/stderr file. Linux `wait4` avoids accidentally
reporting the accumulated peak RSS of earlier child processes.

Every mapping must emit the independently counted 13,124 primary records and
match the warmup's canonical SAM digest. Canonicalization ignores only program
headers, record order, and optional-tag order; it retains alignment coordinates,
CIGAR, flags, sequence, quality, reference headers, tag values, and duplicates.
This checks reproducibility, not biological accuracy against ground truth.
Canonicalization holds the SAM in memory outside the timed process; this small
pilot runner is not intended for whole-human-genome SAM files.

OS page cache is uncontrolled; these are warmed runs, not cold-disk measurements.
Index construction has one sample. Do not use it for a precise amortization claim.
No privileged cache flushing or CPU-governor changes are made.

## Local baseline results — 2026-09-21

Intel Core Ultra 5 125H, Linux x86-64, GCC 16.2.1, upstream default `-g -O2`
build, minimap2 2.30-r1287, `map-pb`, one mapping thread, no CPU pinning.
Reference: 4,639,675 bases. Input: 13,124 upstream PacBio example reads.

| Measured run | Mapping wall time (s) | Peak child RSS (MiB) |
|---|---:|---:|
| 1 | 41.076 | 804.07 |
| 2 | 40.224 | 800.97 |
| 3 | 40.685 | 815.77 |
| 4 | 40.022 | 801.96 |
| 5 | 39.961 | 804.77 |

Median mapping time: **40.224 s**. Separate one-time index build: **0.256 s**.
The warmup took 41.116 s and is excluded from the median. All six mapping
outputs match: 14,496 total SAM records, 13,124 primary records, of which
12,800 are mapped. Canonical alignment SHA-256:
`4afd749160a1fe40f4789a04d9aac773ffc7fd8b8bcfdcb1c3c2f88ab2895af5`.

SHA-256 provenance:

| Artifact | SHA-256 |
|---|---|
| Upstream cookbook archive | `2cb3626a3c1cf1d26102aa3a98d0f73c7eac10f40428e4f6075ab9e2f179ac31` |
| `ecoli_ref.fa` | `419863da8fee3c819df3c15cc827a537f2f3bfd1f0577f6f416d50748572bff8` |
| `ecoli_p6_25x_canu.fa` | `7ae0fb74eb707799b5afa822e60b5425f6b1f134c48360cb6fdb7bc1304f5872` |
| Baseline minimap2 binary | `7de89ed252b22891f674d18d0d9e67dd8f49a660572a918c988c05ddd9716661` |

Raw commands, per-run measurements, SAM files, stderr, and host details remain
in `build/mapping-pilot-2026-09-21/`. This ignored build directory is a local
artifact, not a public archive. Preserve/archive it before citing these results
in a manuscript. These measurements contain **no TracEon implementation or
comparison**, so they establish neither a speedup nor research adoption.

## Profiling decision

A separate build with `CFLAGS='-g -O2 -pg -Wall -Wc++-compat'` produced identical
canonical SAM output. GNU gprof's 10 ms samples attributed 0.24 s (0.62% of
sampled self CPU time) to `mm_idx_getseq`, versus 27.55 s (70.68%) to
`ksw_extd2_sse41`. The 3.65 s (9.36%) in `mm_idx_get` belongs to the minimizer
index, not sequence-ID retrieval. The profile is diagnostic: instrumentation,
sampling, inlining, and pipeline threads preclude treating these percentages as
precise wall-time savings. The profiling run is excluded from baseline timings.

Decision: **do not build an internal minimap2 reference adapter on this evidence**.
The potentially replaceable work is too small in this workload. This does not
rule out other datasets or a different retrieval workload.

Profile reproduction, after the baseline:

```sh
git clone --local /tmp/traceon-minimap2-pilot /tmp/traceon-minimap2-profile
make -C /tmp/traceon-minimap2-profile -j4 'CFLAGS=-g -O2 -pg -Wall -Wc++-compat'
mkdir -p build/mapping-profile-reproduction
cd build/mapping-profile-reproduction
/tmp/traceon-minimap2-profile/minimap2 -a -x map-pb -t 1 \
  ../mapping-pilot-2026-09-21/reference.mmi \
  /tmp/traceon-mapping-data/ecoli_p6_25x_canu.fa > profile.stdout 2> profile.stderr
gprof /tmp/traceon-minimap2-profile/minimap2 gmon.out > gprof.txt
```

The original profile, instrumented binary hash, matching SAM digest, and timing
are retained in `build/mapping-profile-2026-09-21/`.

## Initial control: per-process restore before remapping

The pilot integrates TracEon as an external sequence-access step that supplies
ordinary FASTA to an unchanged minimap2. It does not replace minimap2's index.
`traceon_subset` restores a cache, validates requested IDs, borrows sequence
views, and writes the requested records in order. It deliberately supports
FASTA only: silently dropping FASTQ qualities would invalidate the comparison.
Header descriptions are not retained, and the mapper is run without `-y`.

Use 1,000 reads selected with seed 42, preserving source order. Rotate the order
of three arms across repetitions, after one warmup of each:

1. `seqtk subseq` scans the complete FASTA, materializes the subset, then minimap2
   maps it.
2. TracEon restores the full cache, retrieves IDs, materializes the same subset,
   then the same minimap2 maps it.
3. Reuse a subset previously materialized by seqtk and run minimap2 only. This
   is an essential baseline when the requested subset does not change.

Cache construction and saved-subset construction are reported separately.
Steady-state TracEon timings include full cache restoration in each process;
there is no resident service. Sequence and SAM validation happens outside timed
child processes, so sums exclude validation/orchestration overhead. This is a
component-sum workflow measurement, not a stopwatch around the entire Python
script. Peak RSS is the maximum of sequential extraction and mapping children.

```sh
git clone --branch v1.5 --depth 1 https://github.com/lh3/seqtk.git /tmp/traceon-seqtk-pilot
# Pinned seqtk commit: 94e707082d39b0a038f234df676e32d9802c0dc7.
make -C /tmp/traceon-seqtk-pilot -j4
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target traceon_subset traceon_make_cache -j4
python tests/test_subset_mapping_pilot.py build/traceon_make_cache build/traceon_subset
python benchmarks/subset_mapping_pilot.py \
  --minimap2 /tmp/traceon-minimap2-pilot/minimap2 \
  --seqtk /tmp/traceon-seqtk-pilot/seqtk \
  --cache-builder build/traceon_make_cache --extractor build/traceon_subset \
  --index build/mapping-pilot-2026-09-21/reference.mmi \
  --reads /tmp/traceon-mapping-data/ecoli_p6_25x_canu.fa \
  --output build/subset-mapping-pilot-2026-09-21 --count 1000 --repetitions 5
```

The random subset is a feasibility fixture, not an established biological
analysis. A later study needs a motivated selection rule (for example, reads
recruited to particular contigs) and evidence that changing subsets or reuse
across analyses warrants the cache. A fixed subset does not justify repeatedly
extracting it if saving it once is sufficient.

### Subset integration results — 2026-09-21

All 18 mapping outputs (three arms, one warmup and five measured runs each)
match: 1,098 SAM records, 1,000 primary reads, 976 mapped primary reads.
All extracted ID/sequence pairs match the independent FASTA reader in order.
Canonical SAM SHA-256:
`1358ec25464a7d2505115d55afd041754fc88c79996b0e1067abdd1631b29478`.
Selected-ID SHA-256:
`a15cbdcc3c636f5c7f1f8fc8c91adfcd8ad5c4cb741049c4f96fa05d19a4c5de`.

| Arm | Five workflow samples (s) | Median (s) | Median extraction (s) | Max child RSS (MiB) |
|---|---|---:|---:|---:|
| seqtk + minimap2 | 2.099850, 2.164830, 2.186779, 2.221135, 2.174009 | 2.174 | 0.041 | 466.31 |
| TracEon + minimap2 | 2.102644, 2.196028, 2.244985, 2.414009, 2.291141 | 2.245 | 0.067 | 466.75 |
| Saved subset + minimap2 | 2.073456, 2.009665, 2.156060, 2.159870, 2.264916 | 2.156 | 0 | 463.68 |

One-time preparation: TracEon cache construction **29.908 s**; seqtk saved-subset
construction **0.0366 s**. The cache occupies 45,149,428 bytes. These preparation
costs are additional to the steady-state table. The TracEon executable is built
with the repository's Release configuration; seqtk uses its default `-g -O2`.

**Limited conclusion: no demonstrated benefit for this per-process-restore adapter.**
Do not claim a statistically significant slowdown from five close, noisy samples,
but there is no measured speedup that could repay the additional preparation cost.
Saving a fixed subset is the simpler practical solution here. Both this result
and the internal profiling result are specific to this dataset and workload.
In particular, this result must not be generalized to a shared resident cache.

Raw per-phase commands/timings, binary and input hashes, cache, selected IDs,
FASTA outputs, SAM files, and validation summary are retained in
`build/subset-mapping-pilot-2026-09-21/`. This is a mutable development run, not
archived release evidence or evidence of external research adoption.

TracEon already implements a resident server for **one loaded dataset queried by
many independent processes**. The initial adapter did not exercise that design.
The corrected experiment below uses the existing service, without adding a new
server or changing the wire protocol.

## Corrected experiment: cache once, query from many processes

`shared_mapping_pilot.py` builds one cache and starts **one** `remote_bench`
server bound to loopback on an ephemeral port. That same server stays alive
through all process counts and repetitions. Remote clients receive only the
host, port, and ID list; they neither restore the cache nor load the source to
enumerate keys. Each client holds one connection for its 1,000 GET requests.
The client collects the returned sequence payloads and emits FASTA for minimap2.

Compare 1, 2, and 4 independent clients, five measured repetitions after one
warmup per arm/process count. Each process and repetition gets a different
seeded subset, preserving source order; each arm receives exactly the same
subset for a given trial. Rotate arm order across repetitions:

- `seqtk`: independent source scans by each process.
- `local-restore`: independent complete cache restores, retained as a control.
- `shared-server`: independent clients querying the one resident server.

All extraction processes in a group run concurrently, followed by concurrent
minimap2 processes. Measure each stage using one common elapsed interval, not
the sum of concurrent process durations. The reported workflow time is the sum
of these two stage intervals, including process-launch/wait overhead but
excluding oracle validation between stages. This is a staged batch workflow,
not an overlapped extraction/mapping pipeline or an internal mapper adapter.

Sample simultaneous aggregate PSS from Linux `/proc/*/smaps_rollup` approximately
every 5 ms. PSS apportions shared pages instead of counting each process's shared
libraries in full. Include the resident server in the shared arm at both stages;
exclude the Python controller. The server remains alive but idle during control
arms and is excluded from their reported memory. Sampling can miss short peaks,
and sequential reads of `/proc` are not an atomic system snapshot. Report
retrieval and mapping memory separately: minimap2 still owns its native index
and alignment allocations, which TracEon does not share. PSS does not include
unmapped filesystem page-cache or all kernel/socket memory; it is not total
machine memory usage.

Both preparation costs are separate: cache construction once, and server
restore/startup once. Neither is charged per client or per query. A standalone
test removes the source and cache files after server startup, then verifies two
independent clients can still retrieve correct records. Missing IDs fail the
run, and all returned sequences and downstream alignments must match across
arms before results are accepted.

```sh
cmake --build build --target traceon_subset remote_bench traceon_make_cache -j4
python tests/test_shared_mapping_pilot.py build/traceon_make_cache build/traceon_subset build/remote_bench
python benchmarks/shared_mapping_pilot.py \
  --minimap2 /tmp/traceon-minimap2-pilot/minimap2 \
  --seqtk /tmp/traceon-seqtk-pilot/seqtk \
  --cache-builder build/traceon_make_cache --extractor build/traceon_subset \
  --server build/remote_bench \
  --index build/mapping-pilot-2026-09-21/reference.mmi \
  --reads /tmp/traceon-mapping-data/ecoli_p6_25x_canu.fa \
  --output build/shared-mapping-pilot-2026-09-21 \
  --processes 1 2 4 --count 1000 --repetitions 5
```

This experiment tests the actual deployment model. It remains a workload-specific
pilot: loopback networking, one modest FASTA dataset, sequence-only queries,
an existing minimap2 index shared as an input file, and uncontrolled warm OS page
cache. A fixed, reusable read subset remains a different case where saving it
once is a valid simpler baseline. Larger-scale claims also need an established
indexed-access baseline and representative research workloads.

### Shared-server results — 2026-09-21

The same server served **42,000 sequence requests from 42 independent client
processes**, including warmups, across the entire study. Its log contains exactly
one restore, its PID remained unchanged, and shutdown exited successfully.
All 126 extraction outputs matched the independent sequence oracle. The 126
mapping outputs passed per-job equivalence checks across arms (different jobs
intentionally have different ID subsets).

| Processes | Retrieval backend | Median concurrent retrieval (ms) | Max sampled aggregate retrieval PSS (MiB) | Median staged workflow (s) |
|---:|---|---:|---:|---:|
| 1 | seqtk scan | 41.57 | 0.50 | 2.157 |
| 1 | Independent cache restore | 63.29 | 155.42 | 2.135 |
| 1 | Shared TracEon server | 32.07 | 122.63 | 2.096 |
| 2 | seqtk scan | 47.57 | 0.93 | 2.475 |
| 2 | Independent cache restore | 71.05 | 288.84 | 2.592 |
| 2 | Shared TracEon server | 33.36 | 132.03 | 2.448 |
| 4 | seqtk scan | 77.34 | 1.78 | 4.293 |
| 4 | Independent cache restore | 84.21 | 514.73 | 4.448 |
| 4 | Shared TracEon server | 40.71 | 149.47 | 3.269 |

Each process retrieves 1,000 records. Times are medians of five measured trials;
memory is the largest sampled aggregate PSS among those trials, including the
server for the shared arm. It is not a sum of unrelated individual RSS peaks.

Preparation happened once: **29.926 s** to construct the cache, then **0.0577 s**
to restore/start the server. Construction is common to both TracEon deployments;
seqtk requires neither step. The initial server PSS was 113.19 MiB. Cache SHA-256
matched the earlier build:
`39d58f0343c16dac50e3827595cdd764790548f3a6a5521a3f181b3182354197`.

At four processes, the shared design roughly halves retrieval wall time and
reduces retrieval-phase process memory by about 71% relative to separate
restores. Seqtk retains a much smaller memory footprint; the shared cache is
not a universal memory improvement over streaming.

Do **not** attribute the entire staged-workflow difference to TracEon: unchanged
minimap2 mapping-stage times varied substantially. For four processes, mapping
alone ranged from 2.81 to 4.50 s with seqtk input, 2.88 to 4.67 s after local
restore, and 2.93 to 4.54 s after shared retrieval. These differences greatly
exceed the retrieval-time differences. Scheduling on this heterogeneous CPU
was uncontrolled. Stable end-to-end speedup and a break-even reuse count need
a controlled follow-up; they are not established by this pilot.

**Corrected conclusion:** the cache-once, multi-process deployment shows useful
retrieval scaling on this workload. The previous cold adapter cannot justify
rejecting that design. The research question should evaluate shared data reuse,
total server/client resources, and setup amortization, rather than treating
TracEon as a replacement for the mapper's seeding or alignment algorithms.

Raw metadata, all commands and phase measurements, selected IDs, returned FASTA,
SAM outputs, and the server lifecycle log are retained in
`build/shared-mapping-pilot-2026-09-21/`; `summary.json` contains every measured
sample and aggregate. This remains local, mutable feasibility evidence.

## Evidence required for the integration study

1. Name the actual research analysis and identify where repeated retrieval occurs.
   A tool's popularity alone does not establish a suitable workload.
2. Preserve the mapper's parameters, seeding, chaining, alignment, names, record
   order where required, qualities, pairing, and ambiguous-base semantics.
3. Compare equivalent complete outputs to upstream. Add fixtures for reverse
   strands, ambiguous bases, absent IDs, empty inputs, and long records, plus
   representative real data. Explain any normalization of nondeterministic output.
4. Measure baseline versus integrated runs in alternating order, with repeated
   samples. Include both preparation and steady-state work, output writing,
   memory, and repeated-workflow break-even. Separate cold-start experiments
   from warm-cache experiments and identify all excluded costs.
5. Publish the integration patch, exact source revisions, dataset provenance and
   checksums, environment, raw samples, output validation, and failure cases.
   Avoid performance claims until equivalent-output checks pass.

