# Repeated record-access claim check — 2026-09-22

**Supported claim:** TracEon lets independent processes retrieve whole records
from one resident immutable dataset, preserving identifier, sequence, FASTQ
quality, request order and duplicate requests. Batched access reduced CLI-level
retrieval time against samtools indexed access in these tests. Setup and resident
memory remain material costs; the results do not establish a general advantage
over indexed files or other resident services.

This follow-up tests the record service directly. It does not replace a mapper
engine or measure downstream mapping. The separate
[resident mapper control](shared-minimap-index-pilot.md) found no consistent
TracEon lookup-backend advantage.

The subsequent [persistent-reader control](resident-record-study.md) removes
repeated reader initialization and expands to five inputs and 50 conditions.
It retains batched random short-record gains but finds scalar, long-record and
some ordered-access disadvantages. Use both studies when interpreting the claim.

## Design

Two public minimap2 cookbook datasets were used: 10,000 simulated 100-base
single-end E. coli FASTQ reads, and 13,124 corrected PacBio E. coli FASTA reads.
The portable Release TracEon candidate was compared with samtools 1.24
`fqidx`/`faidx`, using persistent file indexes and warm OS caches. Source IDs
were checked for uniqueness; duplicate *requests* were tested separately.

Each dataset was built into a fresh TracEon cache and a fresh samtools index.
One server remained alive across the entire dataset study. Independent clients
used snapshot-pinned scalar requests or batches of up to 64 records. All arms
received the same ID lists and wrote their output to local files. Indexed
extraction disabled output wrapping. No client loaded the full TracEon cache.

The matrix covered 1 and 4 simultaneous client processes, 10 and 1,000 requested
records per process, and four patterns:

- **Changing:** independently sampled IDs change between workers and repetitions.
- **Ordered:** the changing subsets are emitted in source-file order.
- **Repeated:** the same fixed subset is reused across workers and repetitions.
- **Duplicates:** requests sample with replacement from the first 32 IDs.

There were six measured repetitions per cell, covering all six arm-order
permutations, plus one warmup. Patterns/cells ran in a fixed order. Server and
workers shared CPU affinity 8–11; validation and PSS sampling used CPU 16.
The host was not exclusive and CPU frequency was not controlled. The datasets
were scanned for validation before timing, so these are warm-cache experiments.
No other experiment was deliberately run concurrently.

All output records were compared against the source, including order, duplicates
and FASTQ qualities. Timing includes process launch, samtools index loading or
TracEon connection/handshake, request handling and output writing. It is not a
measurement of isolated hash lookup speed. A resident HTSlib-based service was
not compared; gains over per-invocation index loading need not be unique to
TracEon.

## Results

All **672 cohorts / 1,680 retrieval processes / 848,400 delivered records** passed
exact-output checks. This includes 576 measured cohorts / 727,200 delivered
records; the remainder were warmups. Each of the 32 test conditions had a lower
batched median time than the indexed CLI, and every one of the 192 measured
indexed-versus-batched pairs favored batching. Scalar access had a higher median
time than indexed access in 11 conditions.

Four workers, 1,000 requested records per worker; median cohort times:

| Dataset and request pattern | Indexed ms | Scalar ms | Batch-64 ms | Projected setup break-even cohorts |
|---|---:|---:|---:|---:|
| Short, changing | 52.69 | 54.76 | 25.79 | 2.7 |
| Short, ordered | 45.27 | 59.71 | 27.72 | 4.2 |
| Short, repeated | 56.91 | 56.82 | 25.15 | 2.5 |
| Short, duplicates | 37.51 | 52.78 | 24.56 | 5.9 |
| Long, changing | 140.85 | 80.68 | 79.87 | 832 |
| Long, ordered | 126.75 | 83.16 | 77.66 | 933 |
| Long, repeated | 141.17 | 83.04 | 74.61 | 772 |
| Long, duplicates | 122.80 | 85.23 | 78.42 | 1,201 |

For changing subsets at four workers, median sampled peak aggregate PSS was
7.98 MiB indexed versus 8.19 MiB batched for short reads, and 11.47 versus
133.92 MiB for long reads. Batched TracEon used more aggregate PSS in **31 of 32**
conditions. The service shares its cache among clients, but that does not imply
less memory than an indexed file. PSS includes server plus clients for TracEon
and clients only for samtools; the controller and idle server in the indexed arm
are excluded. PSS sampling at 5 ms can miss brief peaks and does not account for
all machine-wide filesystem page-cache memory. Server allocations retained from
earlier cells influence later resident-memory observations.

## Preparation and reuse

| Dataset | Cache build s | Server restore/readiness s | Indexed-file preparation s | Extra TracEon setup s |
|---|---:|---:|---:|---:|
| Short | 0.0813 | 0.0411 | 0.0464 | 0.0759 |
| Long | 49.3772 | 0.2649 | 0.3284 | 49.3137 |

Projected break-even is extra setup divided by the **paired mean** cohort-time
saving; the displayed timing columns are medians. Across all cells it ranged
from 2.5–8.5 short-read cohorts and 772–3,380 long-read cohorts. Fractional
cohorts are estimates, not completed jobs. These projections assume comparable
future jobs and an uninterrupted server lifetime. They were not validated by
running thousands of cohorts, and setup was measured only once per dataset.
Small-query gains include avoiding repeated CLI/index initialization; they do
not establish faster per-record algorithms.

## Interpretation

Present TracEon as shared, repeated whole-record access with explicit snapshot
identity, query semantics and bounded batching. The experiment supports lower
retrieval time for the tested warm, local CLI workflows. It also identifies
reasons to prefer indexed files: little reuse, memory constraints, or workloads
where setup and communication dominate. Small example datasets, one host,
fixed cell ordering, and a limited comparator set prevent broad scaling claims.
No cold-cache, network-distance, memory-pressure or downstream speedup claim is
supported by this follow-up. Missing-ID/error behavior remains covered by the
existing protocol/CLI tests rather than this successful-retrieval matrix.

The next evidence needed for broader claims is a larger independent dataset,
equivalently resident indexed-access alternatives, and an external research
workflow. Those are evidence gaps, not prerequisites for using the existing
record service on an appropriate measured workload.

## Reproduce and inspect

```sh
python benchmarks/record_reuse_study.py \
  --traceon /tmp/traceon-release-check-20260922/traceon \
  --samtools build/research-inputs/samtools-1.24/samtools \
  --reads build/research-inputs/cookbook/ecoli_mason_1.fq --fastq \
  --cpus 8 9 10 11 --controller-cpu 16 \
  --output build/record-reuse-short-new
```

For long reads use `ecoli_p6_25x_canu.fa`, omit `--fastq`, and choose a new output
directory. Select available CPUs on another host. The benchmark needs permission
to open loopback sockets. Both studies ran sequentially after that permission
was granted; the initial sandbox-blocked startup is not a measurement.

Local raw outputs and complete reports are in
`build/record-reuse-{short,long}-2026-09-22-v2/`. The
[evidence archive](results/record-reuse-study-2026-09-22.tar.gz) contains source
snapshots, input/tool hashes, preparation measurements, all validated cohort
metrics and output hashes, request lists, logs, and complete per-cell summaries.
Large retrieved payloads and cache files remain local and can be regenerated.

Evidence archive SHA-256:
`62be8468c58d1be81c3115dafa63b2724fa3d1fce2d38a45c073e9dae8ea4068`.
