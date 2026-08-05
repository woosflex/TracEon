#!/usr/bin/env python3
"""Generate analysis plots from TracEon real-world benchmark data (2026-08-05)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Data from benchmarks/real_world_benchmark_report.md (verified 2026-08-05)
# FASTA rows: (name, size_mb, records, seqkit, biopy, traceon, restore, mem_traceon_mb)
fasta = [
    ("phix174", 0.00, 1, 0.03, 0.13, 0.01, 0.00, 1),
    ("lambda", 0.01, 1, 0.02, 0.14, 0.01, 0.00, 2),
    ("ecoli", 1.32, 1, 0.05, 0.16, 0.03, 0.00, 15),
    ("yeast", 3.67, 17, 0.08, 0.20, 0.07, 0.01, 25),
    ("chr21", 10.60, 1, 0.29, 0.45, 0.23, 0.03, 85),
    ("mouse", 35.18, 1, 0.72, 1.02, 0.67, 0.08, 235),
    ("arabid", 35.75, 7, 0.53, 0.95, 0.61, 0.08, 235),
    ("fruitfly", 42.19, 1870, 0.61, 1.15, 0.71, 0.09, 281),
    ("zebraf", 48.09, 3, 0.74, 1.31, 0.84, 0.10, 308),
    ("chr1", 62.38, 1, 1.58, 1.88, 1.31, 0.16, 480),
    ("chr1_3", 181.96, 3, 3.20, 5.04, 2.76, 0.43, 1647),
]

# FASTQ rows: (name, size_mb, records, seqkit, biopy, traceon, restore, traceon_ops, pyfastx_ops, mem_t, mem_p)
fastq = [
    ("WGS 100MB", 89.32, 299593, 0.06, 2.38, 0.10, 0.09, 13_886_083, 15_886_913, 146, 11),
    ("WGS .gz", 14.98, 299593, 0.20, 2.66, 0.32, 0.09, 14_931_288, 15_548_164, 164, 10),
    ("Nanopore 500MB", 499.63, 8730, 0.14, 3.09, 0.41, 0.30, 48_991_766, 207_205, 522, 12),
]

names = [r[0] for r in fasta]
sizes = [r[1] for r in fasta]
seqkit_f = [r[3] for r in fasta]
biopy_f = [r[4] for r in fasta]
traceon_f = [r[5] for r in fasta]
restore_f = [r[6] for r in fasta]
mem_t = [r[7] for r in fasta]

fq_names = [r[0] for r in fastq]
fq_seqkit = [r[3] for r in fastq]
fq_biopy = [r[4] for r in fastq]
fq_traceon = [r[5] for r in fastq]
fq_t_ops = [r[7] for r in fastq]
fq_p_ops = [r[8] for r in fastq]
fq_mem_t = [r[9] for r in fastq]
fq_mem_p = [r[10] for r in fastq]

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})

# ── Plot 1: FASTA parse time (log scale) ──────────────────────────────
fig, ax = plt.subplots(figsize=(11, 5.5))
x = np.arange(len(names))
w = 0.27
ax.bar(x - w, biopy_f, w, label="BioPython", color="#d62728")
ax.bar(x, seqkit_f, w, label="SeqKit", color="#ff7f0e")
ax.bar(x + w, traceon_f, w, label="TracEon", color="#1f77b4")
ax.set_xticks(x); ax.set_xticklabels(names, rotation=45, ha="right")
ax.set_ylabel("Parse time (s)")
ax.set_title("FASTA parse time by tool (real NCBI genomes, smaller = better)")
ax.legend()
fig.tight_layout()
fig.savefig("plot1_fasta_parse.png", dpi=130)
plt.close(fig)

# ── Plot 2: FASTQ parse time + speedup ────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5))
x = np.arange(len(fq_names))
w = 0.27
ax.bar(x - w, fq_biopy, w, label="BioPython", color="#d62728")
ax.bar(x, fq_seqkit, w, label="SeqKit", color="#ff7f0e")
ax.bar(x + w, fq_traceon, w, label="TracEon", color="#1f77b4")
ax.set_xticks(x); ax.set_xticklabels(fq_names)
ax.set_ylabel("Parse time (s)")
ax.set_title("FASTQ parse time (realistic reads, smaller = better)")
ax.legend()
fig.tight_layout()
fig.savefig("plot2_fastq_parse.png", dpi=130)
plt.close(fig)

# ── Plot 3: Lookup throughput short vs long reads ─────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
labels = ["WGS 100MB\n(short reads)", "WGS .gz\n(short reads)", "Nanopore 500MB\n(long reads)"]
t_ops = [v/1e6 for v in fq_t_ops]
p_ops = [v/1e6 for v in fq_p_ops]
x = np.arange(3); w = 0.35
bars1 = ax.bar(x - w/2, t_ops, w, label="TracEon", color="#1f77b4")
bars2 = ax.bar(x + w/2, p_ops, w, label="PyFastX", color="#2ca02c")
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("Million ops/s")
ax.set_title("Lookup throughput (higher = better)")
for b in list(bars1)+list(bars2):
    ax.annotate(f"{b.get_height():.1f}", (b.get_x()+b.get_width()/2, b.get_height()),
                ha="center", va="bottom", fontsize=8)
ax.legend()
fig.tight_layout()
fig.savefig("plot3_lookup_throughput.png", dpi=130)
plt.close(fig)

# ── Plot 4: Memory usage ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5.5))
x = np.arange(len(names))
ax.bar(x, mem_t, color="#1f77b4", label="TracEon peak RSS")
ax.bar(x, [s*6 for s in sizes], width=0.4, color="#999999", alpha=0.6, label="~6x compressed size (reference)")
ax.set_xticks(x); ax.set_xticklabels(names, rotation=45, ha="right")
ax.set_ylabel("Peak RSS (MB)")
ax.set_title("TracEon memory use on FASTA (arena holds decompressed data)")
ax.legend()
fig.tight_layout()
fig.savefig("plot4_memory.png", dpi=130)
plt.close(fig)

# ── Plot 5: Restore vs parse speedup ──────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))
speedups = [traceon_f[i]/restore_f[i] if restore_f[i] > 0 else 0 for i in range(len(names))]
ax.bar(names, speedups, color="#1f77b4")
ax.set_ylabel("Parse time / restore time (x)")
ax.set_title("Binary cache restore speedup (parse ÷ restore, higher = better)")
ax.tick_params(axis="x", rotation=45)
fig.tight_layout()
fig.savefig("plot5_restore_speedup.png", dpi=130)
plt.close(fig)

print("Plots written:")
import os
for f in ["plot1_fasta_parse.png","plot2_fastq_parse.png","plot3_lookup_throughput.png","plot4_memory.png","plot5_restore_speedup.png"]:
    print(" ", f, os.path.getsize(f)//1024, "KB")
