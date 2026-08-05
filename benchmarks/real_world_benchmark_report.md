# TracEon Real-World Benchmark Report

Regenerated fresh from two harnesses:
- `benchmarks/real_comparative_benchmark.py` — TracEon vs SeqTk/SeqKit/BioPython/PyFastX on real NCBI genomes + realistic FASTQ files.
- `benchmarks/integrated_pipeline_benchmark.py` — total wall-clock time for a simulated repeated-use pipeline (parse once, reuse many times) vs re-doing the work from scratch each time with baseline tools.

Environment: `traceon_bench` micromamba env (python 3.13.11, biopython 1.85, pyfastx 2.3.1, psutil 7.0.0, seqkit 2.12.0, seqtk 1.5). Driver: `build/traceon_driver` (x86-64, `CMAKE_BUILD_TYPE=Release`).

All 14 datasets were present from the prior download/generation session: 11 real NCBI genomes in `real_data/` (gzip-validated) and 3 synthetic-but-realistic FASTQ files in `benchmarks/` (sequential `@r0..@rN` headers, which is the ID scheme the `lookup` subcommand and the PyFastX shim both generate — see script docstrings). No datasets had to be (re)generated; both scripts ran as-is with zero code changes.

Single-run wall-clock numbers on a shared workstation; expect run-to-run noise of roughly ±20-50% on the smaller rows and several-fold on PyFastX's long-read lookup (see note in §1).

## 1. Comparative benchmark (real data)

All times in seconds unless noted. FASTA rows are real NCBI genomes; FASTQ rows are realistic synthetic reads. Lookup throughput (`ops/s`) requires the `prefix+index` ID scheme, so it is only measured on the FASTQ rows. Peak RSS (MB) measured per subprocess with psutil.

| Dataset | Size (MB) | Records | SeqTk* | SeqKit* | BioPy | **TracEon** | Restore | TracEon ops/s | PyFastX ops/s | Mem TracEon (MB) | Mem PyFastX (MB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| phix174.fasta.gz | 0.00 | 1 | - | 0.03 | 0.13 | **0.01** | 0.00 | n/a | n/a | 1 | - |
| lambda_phage.fasta.gz | 0.01 | 1 | - | 0.02 | 0.14 | **0.01** | 0.00 | n/a | n/a | 2 | - |
| ecoli_quick.fasta.gz | 1.32 | 1 | - | 0.05 | 0.16 | **0.03** | 0.00 | n/a | n/a | 15 | - |
| yeast.fasta.gz | 3.67 | 17 | - | 0.08 | 0.20 | **0.07** | 0.01 | n/a | n/a | 25 | - |
| human_chr21.fasta.gz | 10.60 | 1 | - | 0.29 | 0.45 | **0.23** | 0.03 | n/a | n/a | 85 | - |
| mouse_sample.fasta.gz | 35.18 | 1 | - | 0.72 | 1.02 | **0.67** | 0.08 | n/a | n/a | 235 | - |
| arabidopsis.fasta.gz | 35.75 | 7 | - | **0.53** | 0.95 | 0.61 | 0.08 | n/a | n/a | 235 | - |
| fruit_fly.fasta.gz | 42.19 | 1,870 | - | **0.61** | 1.15 | 0.71 | 0.09 | n/a | n/a | 281 | - |
| zebrafish_sample.fasta.gz | 48.09 | 3 | - | **0.74** | 1.31 | 0.84 | 0.10 | n/a | n/a | 308 | - |
| human_chr1.fasta.gz | 62.38 | 1 | - | 1.58 | 1.88 | **1.31** | 0.16 | n/a | n/a | 480 | - |
| human_chr1_3.fasta.gz | 181.96 | 3 | - | 3.20 | 5.04 | **2.76** | 0.43 | n/a | n/a | 1,647 | - |
| bench_wgs_100mb.fastq | 89.32 | 299,593 | 0.09 | **0.06** | 2.38 | 0.10 | 0.09 | 13,886,083 | 15,886,913 | 146 | 11 |
| bench_wgs_100mb.fastq.gz | 14.98 | 299,593 | 0.23 | **0.20** | 2.66 | 0.32 | 0.09 | 14,931,288 | 15,548,164 | 164 | 10 |
| bench_Nanopore_500MB.fastq | 499.63 | 8,730 | 2.06 | **0.14** | 3.09 | 0.41 | 0.30 | 48,991,766 | 207,205 | 522 | 12 |

\* SeqTk (`fqchk`) and SeqKit (`stats`) compute per-file statistics (quality/base-composition profiles for fqchk), which is more work than TracEon's parse-only `load`; the two columns are not strictly apples-to-apples with the other three.

**Takeaways (honest, including negatives):**
- **TracEon beats BioPython on every dataset**: 1.4x (human_chr1) to 24x (WGS FASTQ) faster than a plain Python parse.
- **vs SeqKit the parse picture is mixed — TracEon does NOT win everywhere.** TracEon wins on large single-contig FASTA (human_chr1 1.31 vs 1.58, human_chr1_3 2.76 vs 3.20) and tiny files, but SeqKit is faster on multi-contig FASTA (arabidopsis 0.53 vs 0.61, fruit_fly 0.61 vs 0.71, zebrafish 0.74 vs 0.84) and **2-3x faster on the FASTQ rows** (WGS 0.06 vs 0.10, WGS.gz 0.20 vs 0.32, Nanopore **0.14 vs 0.41**). SeqKit's Rust core and streaming design handle record-dense and long-read FASTQ better than TracEon's parse-and-buffer-everything path. This is a real negative result: for one-shot parsing of FASTQ, SeqKit is the faster tool.
- **Lookup throughput: tied on short reads, massively faster on long reads.** On 100bp-class WGS reads TracEon (13.9M–14.9M ops/s) is essentially tied with PyFastX (15.5M–15.9M), with PyFastX marginally ahead on the plain file. On the Nanopore-shaped file (8,730 records, avg ~57KB) TracEon sustains **49.0M ops/s vs PyFastX's 207K ops/s — ~236x faster**. PyFastX's index cost scales with the length of the record being fetched; TracEon's hash-map lookup is size-agnostic. Note this number is noisy: the prior session measured 560K ops/s for PyFastX on the same file, still a 58x gap — either way it is a decisive TracEon win on long reads.
- **Restore from binary cache** is 3-10x faster than a fresh parse (0.43s vs 2.76s on human_chr1_3; 0.09s vs 0.32s on WGS.gz) — the intended fast path once a cache exists.
- **Memory is the real trade-off.** TracEon's peak RSS runs 5-10x the compressed file size on FASTA (1,647 MB for the 182 MB human_chr1_3 = +805%; 480 MB for 62 MB human_chr1) because the cache holds decompressed sequence in an arena. On the FASTQ rows the overhead is smaller (+64% on WGS, +5% on Nanopore — per-record overhead amortizes against long reads). PyFastX stays at 10-12 MB by mmapping and lazy-parsing. TracEon trades memory for its lock-free, zero-copy read path; the numbers above quantify that cost.

## 2. Integrated pipeline benchmark

Simulates realistic repeated-use workflows (M=5 simulated pipeline runs) rather than one-shot ops.

### A. Repeated full-file reload, no lookups (FASTA)

| Dataset | Scenario | TracEon total | BioPython total | Result |
|---|---|---:|---:|---|
| yeast.fasta.gz | 1 parse+save + 5 restores | 2.88s | 1.18s | **0.4x — TracEon slower** |
| human_chr21.fasta.gz | 1 parse+save + 5 restores | 9.71s | 2.24s | **0.2x — TracEon slower** |
| human_chr1.fasta.gz | 1 parse+save + 5 restores | 58.80s | 9.11s | **0.2x — TracEon slower** |

Confirmed negative result, reproduced from the prior session: `saveBinary()` picks LZ4HC compression for large DNA/RNA payloads (`SmartStrategy::selectCompressionStrategy`), which is expensive up front and dominates the total. For a pure "reload the same file repeatedly, never query it" workload there are no ID lookups to amortize that cost against, so plain re-parsing with BioPython wins at M=5 (3-5x). The binary cache's value is in what you *do* with it afterward (fast restores + O(1) lookups), not in reload speed alone.

### B. Repeated load+lookup pipeline (FASTQ, 50,000 lookups × 5 runs)

`bench_wgs_100mb.fastq.gz`, 299,593 records:

| Approach | Total time (5 runs) | vs TracEon |
|---|---:|---|
| **TracEon** (cache once, restore+lookup ×5) | **8.16s** | — |
| PyFastX (index once, lookup ×5) | 180.87s | **22.2x slower** |
| BioPython (full re-parse ×5, no lookup index) | 6.20s | 0.8x — BioPython slightly faster in this no-lookup-needed comparison* |

\* BioPython's number only covers re-parsing, since it has no queryable index — it cannot answer the 50,000 lookups per run at all. It is included as the "cost of having no cache and no index" floor, not as a fair like-for-like against TracEon's full parse+cache+lookup pipeline.

**Bottom line:** TracEon's integrated advantage is workload-dependent and large when the task involves repeated ID lookups against cached data (~22x over PyFastX here on short reads; the per-lookup gap is far larger on long reads, §1). But it is *not* a universal win: for reload-only, no-query workflows the LZ4HC save cost isn't recouped within 5 repetitions (0.2-0.4x), and for one-shot parsing SeqKit beats TracEon on FASTQ. The strongest honest case for TracEon remains: parse once, cache, then serve O(1) lookups from a lock-free immutable arena — especially with long-read or large-record datasets where indexing competitors degrade.

## Reproducing

```bash
cmake --build build --target traceon_driver -j
cd benchmarks
micromamba run -n traceon_bench python3 real_comparative_benchmark.py
micromamba run -n traceon_bench python3 integrated_pipeline_benchmark.py
```

All 14 datasets must exist: 11 real genomes in `real_data/` and the 3 sequential-ID FASTQ files in `benchmarks/` (299,593 reads ~89 MB for `bench_wgs_100mb.fastq`/`.gz`, 8,730 long reads ~500 MB for `bench_Nanopore_500MB.fastq`). The scripts clean up their own sidecar/cache files (`.traceon_bench_tmp`, `.fxi`, `.fai`); source files are never deleted.
