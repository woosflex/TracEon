# Repeated subset remapping with one resident dataset

TracEon supplies complete records by ID to independent mapping jobs. Each worker
retrieves its selected reads from the same resident snapshot and runs an unchanged
minimap2. This supports repeated inspection/remapping of selected read cohorts;
it does not replace minimap2's seed index, chaining, or alignment algorithms.

## Using the installed service

Build/install as described in the README, then prepare one cache and start one
server. For example, the public minimap2 SIRV Iso-Seq example supplies FASTQ reads
and a matching reference suitable for the `splice:hq` preset:

```bash
traceon build SIRV_iso-seq.fq sirv.traceon
traceon serve sirv.traceon --port 9876
```

Wait for `TRACEON_SERVER READY`, and record its `snapshot=` digest from the trusted
launch. The server hashes/restores only at startup and stays running across jobs.
Each worker needs the endpoint, expected digest, and an ID list; it does not need
the full reads or cache file. In separate worker processes, run:

```bash
traceon query --port 9876 --expect-snapshot "$SNAPSHOT_SHA256" \
  --fastq --ids selected.ids > selected.fastq &&
minimap2 -ax splice:hq -t 1 SIRV_E2.fa selected.fastq > selected.sam
```

Give each worker distinct output paths. Check the query exit status before
mapping: the CLI can have emitted earlier records before a later ID fails. FASTQ
qualities survive retrieval; descriptions and original plus-line annotations do
not. Plain `query` produces FASTA and is suitable for FASTA source datasets such
as the cookbook's corrected E. coli PacBio reads (`map-pb`). This example is
single-end; it does not infer paired-read matching from filenames or IDs.

## Reproducible integration and comparison

`benchmarks/indexed_mapping_study.py` runs the complete workflow with the installed
CLI, then compares it with samtools faidx/fqidx retrieval. It builds each cache,
FAI, and mapping index once; one server spans all workers/repetitions. It changes
subsets deterministically, preserves requested order, checks every returned
sequence/quality against the source, and compares canonical SAM alignments
between arms. Only program-provenance headers, alignment-record order, and
optional-tag order are excluded from SAM comparison.

Example (choose CPU numbers available on your machine):

```bash
python3 benchmarks/indexed_mapping_study.py \
  --traceon /path/to/install/bin/traceon --samtools /path/to/samtools \
  --minimap2 /path/to/minimap2 --reference SIRV_E2.fa --reads SIRV_iso-seq.fq \
  --fastq --preset splice:hq --cpus 8 9 10 11 --controller-cpu 16 \
  --processes 1 2 4 --counts 10 50 --repetitions 6 --batch-size 64 \
  --output new-study-directory
```

The output includes exact commands, input/binary/harness hashes, tool versions,
CPU topology/affinity, setup costs, raw per-stage measurements, and validated SAM
summaries. Process wall time includes startup; it is not per-request tail latency.
The same CPU pool is available to all three arms, including the server in the
shared arms. Six repetitions cover all six arm-order permutations; without
`--batch-size`, two arms require an even number of repetitions.
The controller can be pinned outside that pool. Affinity does not reserve
cores or fix CPU frequency, so short timing differences need repeated examination.

PSS samples are simultaneous process memory, with the server counted for the
shared arm. The idle server remains running while the indexed arm is measured
but is excluded from that arm's deployment total. Controller memory and unallocated
filesystem page-cache pages are excluded. Inputs are read for validation before
measurement, so this is a warm-page-cache comparison, not cold disk I/O. Report
those boundaries rather than treating process PSS as total machine memory.

## Public inputs and baseline

The [minimap2 cookbook archive](https://github.com/lh3/minimap2/releases/download/v2.10/cookbook-data.tgz)
contains the E. coli and SIRV examples. The SIRV FASTQ has 53 records; the E. coli
corrected-read FASTA has 13,124. The archive also contains Mason-simulated short
reads, which must be labelled simulated and are distinct from experimental data.

The checked integration uses minimap2 2.30 at commit
`79c9cc186b95f50bd899f69b48eba995ced810c6` and samtools/HTSlib 1.24, built from their
public sources. [samtools faidx](https://www.htslib.org/doc/samtools-faidx.html) and
[fqidx](https://www.htslib.org/doc/samtools-fqidx.html) provide the indexed baseline;
`-n 0` emits unwrapped records and `--fai-idx` keeps study indexes beside results.
The samtools documentation cautions that fqidx indexes for millions of short
reads can become expensive. The current example datasets do not establish
behavior at that scale.

The successful SIRV smoke run is in `build/installed-sirv-smoke-2026-09-21/`.
It establishes integration correctness, including qualities and matching SAM
alignments, not a general performance advantage or independent biological finding.
