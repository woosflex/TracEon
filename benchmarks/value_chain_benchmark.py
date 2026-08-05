#!/usr/bin/env python3
"""
value_chain_benchmark.py -- Complete apples-to-apples value-chain comparison
for TracEon (parse + index + cache + lookup), fixing the unfair comparison in
real_world_benchmark_report.md (TracEon 'load' = parse+index+buffer vs
SeqKit 'stats' = streaming parse + stats computation).

Three measurements:

  A. FAIR PARSE COMPARISON (parse-only vs parse+index)
     - seqkit seq  <file>  : streaming parse only (no stats, no index).
                             Output discarded to /dev/null; the cost of
                             record formatting is still included.
     - traceon_driver load : parse + full hash-index build + arena buffering.
     - seqkit faidx <file> : FASTA index build (closest competitor feature to
                             TracEon's in-memory index). v2.12.0 rejects
                             gzipped input, so this runs on one-time
                             decompressed copies (decompress time reported
                             separately as the SeqKit path's preprocessing).
     The delta (TracEon load - seqkit seq) is an UPPER BOUND on TracEon's
     index-build + arena-copy cost, stated with caveats: the two tools use
     different decompressors (zlib-ng vs flate2) and seqkit seq spends CPU
     formatting output lines even when discarded.

  B. END-TO-END REPEATED-USE WORKFLOW (the real use case)
     5 simulated pipeline runs of "load genome, answer N lookups":
       FASTQ (bench_wgs_100mb.fastq.gz short reads, bench_Nanopore_500MB.fastq
       long reads), N = 50,000 lookups/run:
         - TracEon:    parse+save once, then per run: restore + 50k lookups
                       (lookup subcommand restores in-process, no double count)
         - SeqKit:     per run: seqkit grep -f ids.txt (scan-based lookup --
                       SeqKit cannot random-access FASTQ), plus the
                       seqkit seq parse-only floor
         - PyFastX:    index once, per run: fresh process, 50k lookups
                       (fresh process = cold per-record cache, the honest
                       pipeline-restart regime; in-process warm lookups are
                       far faster -- noted in the report)
         - BioPython:  per run: full re-parse (no index, no cache -- floor)
       FASTA (yeast.fasta.gz, 17 contigs), N = 1,000 lookups/run (realistic
       for genome FASTA: query a handful of contigs; seqkit faidx's per-ID
       cost is linear ~0.87 ms/ID so 50k would take 43 s/run -- measured and
       reported as a scaling note):
         - TracEon:    load+save once, per run: restore + 1k lookups
         - SeqKit:     one-time decompress + faidx index, per run: batch
                       lookup via `seqkit faidx -l ids.txt`
         - PyFastX:    index once, per run: fresh process, 1k lookups of
                       real accession IDs
         - BioPython:  per run: full re-parse (floor)

Methodology (honest):
  - Every measurement is a fresh subprocess, wall-clock, single run. Shared
    workstation; expect ±20-50% noise on small rows. OS page cache is NOT
    purged between runs (realistic warm-disk pipelines).
  - random.seed(42) for every generated ID list.
  - Sidecar/cache files (.traceon, .fxi, .fai, .gzi) and the .valuechain_tmp/
    decompressed copies are cleaned up; real_data/ sources are never touched.

Usage:  micromamba run -n traceon_bench python3 value_chain_benchmark.py
"""

import os
import sys
import time
import random
import shutil
import subprocess
import gzip

# === CONFIGURATION ===
TRACEON_BINARY = "../build/traceon_driver"
SEQKIT = "seqkit"
REAL_DATA_DIR = "../real_data"
TMP_DIR = ".valuechain_tmp"

M_RUNS = 5
LOOKUPS_FASTQ = 50_000          # task spec: 50k lookups/run on FASTQ
LOOKUPS_FASTA = 1_000           # scaled for FASTA (see docstring)

FASTA_DATASETS = [
    "phix174.fasta.gz", "lambda_phage.fasta.gz", "ecoli_quick.fasta.gz",
    "yeast.fasta.gz", "human_chr21.fasta.gz", "mouse_sample.fasta.gz",
    "arabidopsis.fasta.gz", "fruit_fly.fasta.gz", "zebrafish_sample.fasta.gz",
    "human_chr1.fasta.gz", "human_chr1_3.fasta.gz",
]
FASTQ_DATASETS = [
    "bench_wgs_100mb.fastq", "bench_wgs_100mb.fastq.gz",
    "bench_Nanopore_500MB.fastq",
]

BP_PARSE_SCRIPT = (
    "import sys, gzip; from Bio import SeqIO; "
    "open_func = gzip.open if sys.argv[1].endswith('.gz') else open; "
    "n = sum(1 for _ in SeqIO.parse(open_func(sys.argv[1], 'rt'), sys.argv[2]))"
)

PYFASTX_BUILD_SCRIPT = """
import sys, pyfastx
fname = sys.argv[1]
db = pyfastx.Fastq(fname) if 'fastq' in fname else pyfastx.Fasta(fname)
"""

PYFASTX_LOOKUP_PREFIX_SCRIPT = """
import sys, random, pyfastx
fname, n, prefix, m = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
db = pyfastx.Fastq(fname) if 'fastq' in fname else pyfastx.Fasta(fname)
random.seed(42)
for _ in range(n):
    try: _ = str(db[f"{prefix}{random.randint(0, m-1)}"])
    except: pass
"""

PYFASTX_LOOKUP_REALIDS_SCRIPT = """
import sys, random, pyfastx
fname, n = sys.argv[1], int(sys.argv[2])
db = pyfastx.Fasta(fname)
ids = list(db.keys())
random.seed(42)
for _ in range(n):
    try: _ = str(db[random.choice(ids)])
    except: pass
"""


# === UTILITIES ===
def run_timed(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL):
    """Fresh subprocess, wall-clock. Returns (seconds, stdout_text).
    Wall-clock is used uniformly for every tool (seqkit, traceon_driver,
    python), so process-spawn overhead affects all rows the same way."""
    t0 = time.time()
    r = subprocess.run(cmd, stdout=stdout, stderr=stderr, text=True)
    return time.time() - t0, (r.stdout or "")


def count_records_and_avglen(fname, fmt):
    opener = gzip.open if fname.endswith(".gz") else open
    n = 0
    total_len = 0
    cur_len = 0
    is_fastq = "fastq" in fmt
    with opener(fname, "rt") as f:
        if is_fastq:
            for i, line in enumerate(f):
                if i % 4 == 0:
                    n += 1
                elif i % 4 == 1:
                    total_len += len(line.strip())
        else:
            for line in f:
                if line.startswith(">"):
                    n += 1
                    if cur_len:
                        total_len += cur_len
                        cur_len = 0
                else:
                    cur_len += len(line.strip())
            if cur_len:
                total_len += cur_len
    return n, (total_len / n) if n else 1


def extract_fasta_ids(plain_path):
    ids = []
    with open(plain_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def decompress_copy(gz_path, dst):
    t0 = time.time()
    with gzip.open(gz_path, "rb") as src, open(dst, "wb") as out:
        shutil.copyfileobj(src, out, length=1 << 20)
    return time.time() - t0


def gen_ids_file(path, ids, n, seed=42):
    rnd = random.Random(seed)
    with open(path, "w") as f:
        for _ in range(n):
            f.write(rnd.choice(ids) + "\n")


def cleanup():
    if os.path.isdir(TMP_DIR):
        shutil.rmtree(TMP_DIR, ignore_errors=True)
    for d in [REAL_DATA_DIR, "."]:
        for suf in [".fai", ".gzi", ".fxi", ".traceon_bench_tmp"]:
            for f in os.listdir(d):
                if f.endswith(suf):
                    try:
                        os.remove(os.path.join(d, f))
                    except OSError:
                        pass


# === SECTION A: fair parse comparison ===
def section_a():
    print("\n" + "=" * 118)
    print("A. FAIR PARSE COMPARISON: seqkit seq (parse-only) vs TracEon load (parse+index) vs seqkit faidx (FASTA index build)")
    print("=" * 118)
    hdr = (f"{'DATASET':<26} | {'SIZE':>8} | {'RECS':>9} | "
           f"{'SEQKIT SEQ':>10} | {'TRACEON LOAD':>12} | {'DELTA(est idx)':>13} | {'idx frac':>8} | "
           f"{'SEQKIT FAIDX':>12} | {'DECOMP(1x)':>10}")
    print(hdr)
    print("-" * 118)

    rows = []
    for fname in FASTA_DATASETS:
        path = os.path.join(REAL_DATA_DIR, fname)
        if not os.path.exists(path):
            continue
        size_mb = os.path.getsize(path) / (1024 * 1024)
        n, _ = count_records_and_avglen(path, "fasta")

        t_seq, _ = run_timed([SEQKIT, "seq", path])                       # parse-only
        t_load, _ = run_timed([TRACEON_BINARY, "load", path])             # parse+index

        plain = os.path.join(TMP_DIR, fname.replace(".gz", ""))
        t_dec = decompress_copy(path, plain)
        t_faidx, _ = run_timed([SEQKIT, "faidx", plain])                  # index build

        delta = max(t_load - t_seq, 0.0)
        frac = (delta / t_load * 100) if t_load > 0 else 0.0
        rows.append((fname, size_mb, n, t_seq, t_load, delta, frac, t_faidx, t_dec))
        print(f"{fname:<26} | {size_mb:>8.2f} | {n:>9,} | {t_seq:>10.3f} | {t_load:>12.3f} | "
              f"{delta:>13.3f} | {frac:>7.0f}% | {t_faidx:>12.3f} | {t_dec:>10.3f}")

    for fname in FASTQ_DATASETS:
        if not os.path.exists(fname):
            continue
        size_mb = os.path.getsize(fname) / (1024 * 1024)
        n, _ = count_records_and_avglen(fname, "fastq")

        t_seq, _ = run_timed([SEQKIT, "seq", fname])
        t_load, _ = run_timed([TRACEON_BINARY, "load", fname])

        delta = max(t_load - t_seq, 0.0)
        frac = (delta / t_load * 100) if t_load > 0 else 0.0
        rows.append((fname, size_mb, n, t_seq, t_load, delta, frac, None, None))
        print(f"{fname:<26} | {size_mb:>8.2f} | {n:>9,} | {t_seq:>10.3f} | {t_load:>12.3f} | "
              f"{delta:>13.3f} | {frac:>7.0f}% | {'n/a (FASTQ)':>12} | {'-':>10}")
    print("-" * 118)
    return rows


# === SECTION B: end-to-end repeated-use workflow ===
def fastq_e2e(fname, fmt):
    print("\n" + "-" * 118)
    n_seqs, _ = count_records_and_avglen(fname, fmt)
    print(f"B. END-TO-END (FASTQ): {fname} — {n_seqs:,} records, "
          f"{LOOKUPS_FASTQ:,} lookups/run x {M_RUNS} runs")
    print("-" * 118)
    bin_file = os.path.join(TMP_DIR, os.path.basename(fname) + ".traceon")
    clean_fmt = "fastq"
    prefix = "r"

    # --- TracEon: parse+save once, then M x (restore + lookup) in-process.
    # The `lookup` subcommand restores the cache internally, so each run is
    # a single process (restore + N lookups) -- no double-counted restores.
    t_parse, _ = run_timed([TRACEON_BINARY, "load", fname])
    t_save, _ = run_timed([TRACEON_BINARY, "save", fname, bin_file])
    t_runs = []
    for _ in range(M_RUNS):
        t_runs.append(run_timed([TRACEON_BINARY, "lookup", bin_file,
                                 str(LOOKUPS_FASTQ), prefix, str(n_seqs)])[0])
    tra_total = t_parse + t_save + sum(t_runs)

    # --- SeqKit: scan-based lookup (grep) per run + parse-only floor ---
    ids_all = [f"{prefix}{i}" for i in range(n_seqs)]
    ids_file = os.path.join(TMP_DIR, os.path.basename(fname) + ".ids")
    gen_ids_file(ids_file, ids_all, LOOKUPS_FASTQ)
    seqkit_grep = 0.0
    for _ in range(M_RUNS):
        seqkit_grep += run_timed([SEQKIT, "grep", "-f", ids_file, fname])[0]
    seqkit_seq_floor = 0.0
    for _ in range(M_RUNS):
        seqkit_seq_floor += run_timed([SEQKIT, "seq", fname])[0]

    # --- PyFastX: index once, then M x (fresh process, 50k lookups) ---
    t_pfx_build = run_timed([sys.executable, "-c", PYFASTX_BUILD_SCRIPT, fname])[0]
    pfx_runs = []
    for _ in range(M_RUNS):
        pfx_runs.append(run_timed([sys.executable, "-c", PYFASTX_LOOKUP_PREFIX_SCRIPT,
                                   fname, str(LOOKUPS_FASTQ), prefix, str(n_seqs)])[0])
    pfx_total = t_pfx_build + sum(pfx_runs)

    # --- BioPython: full re-parse per run (floor, no lookup capability) ---
    bio_total = 0.0
    for _ in range(M_RUNS):
        bio_total += run_timed([sys.executable, "-c", BP_PARSE_SCRIPT, fname, clean_fmt])[0]

    print(f"  TracEon   (parse {t_parse:.3f}s + save {t_save:.3f}s + "
          f"{M_RUNS}x[restore+{LOOKUPS_FASTQ:,} lookup {sum(t_runs):.3f}s]): "
          f"{tra_total:>8.2f}s total")
    print(f"  SeqKit    (scan-based lookup via grep, {M_RUNS}x full-file scan): "
          f"{seqkit_grep:>8.2f}s total")
    print(f"  SeqKit    (parse-only floor via seq, {M_RUNS}x): "
          f"{seqkit_seq_floor:>8.2f}s total")
    print(f"  PyFastX   (index {t_pfx_build:.2f}s once + {M_RUNS}x fresh-proc lookup "
          f"{sum(pfx_runs):.2f}s): {pfx_total:>8.2f}s total")
    print(f"  BioPython (full re-parse x{M_RUNS}, no lookup): "
          f"{bio_total:>8.2f}s total")
    print(f"  Ratios vs TracEon: SeqKit-grep {seqkit_grep/tra_total:.1f}x | "
          f"SeqKit-seq {seqkit_seq_floor/tra_total:.1f}x | "
          f"PyFastX {pfx_total/tra_total:.1f}x | "
          f"BioPython {bio_total/tra_total:.1f}x")
    return tra_total, seqkit_grep, seqkit_seq_floor, pfx_total, bio_total


def fasta_e2e(fname, fmt):
    print("\n" + "-" * 118)
    path = os.path.join(REAL_DATA_DIR, fname)
    n_seqs, _ = count_records_and_avglen(path, fmt)
    print(f"B. END-TO-END (FASTA): {fname} — {n_seqs:,} contigs, "
          f"{LOOKUPS_FASTA:,} lookups/run x {M_RUNS} runs (scaled down from 50k: "
          f"seqkit faidx lookup cost is linear ~0.87ms/ID, see report)")
    print("-" * 118)
    bin_file = os.path.join(TMP_DIR, fname + ".traceon")
    plain = os.path.join(TMP_DIR, fname.replace(".gz", ""))

    # --- TracEon (lookup subcommand restores in-process; wall-clock) ---
    t_parse, _ = run_timed([TRACEON_BINARY, "load", path])
    t_save, _ = run_timed([TRACEON_BINARY, "save", path, bin_file])
    t_runs = []
    for _ in range(M_RUNS):
        t_runs.append(run_timed([TRACEON_BINARY, "lookup", bin_file,
                                 str(LOOKUPS_FASTA), "", str(n_seqs)])[0])
    tra_total = t_parse + t_save + sum(t_runs)

    # --- SeqKit faidx: one-time decompress + index, then M x batch lookup ---
    t_dec = decompress_copy(path, plain)
    t_idx = run_timed([SEQKIT, "faidx", plain])[0]
    ids = extract_fasta_ids(plain)
    ids_file = os.path.join(TMP_DIR, fname + ".ids")
    gen_ids_file(ids_file, ids, LOOKUPS_FASTA)
    faidx_runs = []
    for _ in range(M_RUNS):
        faidx_runs.append(run_timed([SEQKIT, "faidx", plain, "-l", ids_file])[0])
    faidx_total = t_dec + t_idx + sum(faidx_runs)

    # --- PyFastX (real accession IDs, fresh process per run) ---
    t_pfx_build = run_timed([sys.executable, "-c", PYFASTX_BUILD_SCRIPT, path])[0]
    pfx_runs = []
    for _ in range(M_RUNS):
        pfx_runs.append(run_timed([sys.executable, "-c", PYFASTX_LOOKUP_REALIDS_SCRIPT,
                                   path, str(LOOKUPS_FASTA)])[0])
    pfx_total = t_pfx_build + sum(pfx_runs)

    # --- BioPython floor ---
    bio_total = 0.0
    for _ in range(M_RUNS):
        bio_total += run_timed([sys.executable, "-c", BP_PARSE_SCRIPT, path, "fasta"])[0]

    print(f"  TracEon   (parse {t_parse:.3f}s + save {t_save:.3f}s + "
          f"{M_RUNS}x[restore+{LOOKUPS_FASTA:,} lookup {sum(t_runs):.3f}s]): "
          f"{tra_total:>8.2f}s total")
    print(f"  SeqKit    (faidx: decompress {t_dec:.2f}s + index {t_idx:.3f}s + "
          f"{M_RUNS}x{LOOKUPS_FASTA:,}-ID lookup {sum(faidx_runs):.2f}s): "
          f"{faidx_total:>8.2f}s total")
    print(f"  PyFastX   (index {t_pfx_build:.2f}s once + {M_RUNS}x fresh-proc lookup "
          f"{sum(pfx_runs):.2f}s): {pfx_total:>8.2f}s total")
    print(f"  BioPython (full re-parse x{M_RUNS}, no lookup): "
          f"{bio_total:>8.2f}s total")
    print(f"  Ratios vs TracEon: SeqKit-faidx {faidx_total/tra_total:.1f}x | "
          f"PyFastX {pfx_total/tra_total:.1f}x | BioPython {bio_total/tra_total:.1f}x")
    return tra_total, faidx_total, pfx_total, bio_total


def faidx_scaling_probe(plain_path, n_list=(1000, 10000)):
    """seqkit faidx batch-lookup cost scales ~linearly with the number of
    patterns (~0.87 ms/ID on yeast). Documented here so the report's
    extrapolation from 1k to 50k lookups is reproducible."""
    ids = extract_fasta_ids(plain_path)
    print("\n" + "-" * 118)
    print("B.2 seqkit faidx batch-lookup scaling (yeast, plain FASTA, 17 contigs)")
    print("-" * 118)
    for n in n_list:
        ids_file = os.path.join(TMP_DIR, f"scaling_{n}.ids")
        gen_ids_file(ids_file, ids, n)
        t, _ = run_timed([SEQKIT, "faidx", plain_path, "-l", ids_file])
        print(f"  {n:>6,} patterns: {t:.3f}s ({n / t:,.0f} IDs/s)")


def section_b():
    print("\n" + "=" * 118)
    print("B. END-TO-END REPEATED-USE WORKFLOW (load + N lookups, x5 runs)")
    print("=" * 118)
    fastq_e2e("bench_wgs_100mb.fastq.gz", "fastq.gz")
    fastq_e2e("bench_Nanopore_500MB.fastq", "fastq")
    fasta_e2e("yeast.fasta.gz", "fasta.gz")
    faidx_scaling_probe(os.path.join(TMP_DIR, "yeast.fasta"))


def main():
    if not os.path.exists(TRACEON_BINARY):
        print("Compile driver first! (expected ../build/traceon_driver)")
        sys.exit(1)
    if not os.path.isdir(TMP_DIR):
        os.makedirs(TMP_DIR)
    try:
        section_a()
        section_b()
    finally:
        cleanup()


if __name__ == "__main__":
    main()
