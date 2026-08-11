#!/usr/bin/env python3
"""Plots for the COMPLETE TracEon value-chain benchmark story (2026-08-05).
Data from benchmarks/real_world_benchmark_report.md §3 (verified)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})

# ── §3.1 Fair parse: SeqKit seq (parse-only) vs TracEon load (parse+index) ──
# (name, seqkit_seq, traceon_load)
parse_fasta = [
    ("phix174", 0.012, 0.004), ("lambda", 0.012, 0.004), ("ecoli", 0.040, 0.026),
    ("yeast", 0.065, 0.065), ("chr21", 0.290, 0.224), ("mouse", 0.774, 0.696),
    ("arabid", 0.561, 0.609), ("fruitfly", 0.640, 0.711), ("zebraf", 0.843, 0.828),
    ("chr1", 1.776, 1.262), ("chr1_3", 3.219, 2.661),
]
parse_fastq = [
    ("WGS 100MB", 0.061, 0.093), ("WGS .gz", 0.195, 0.318), ("Nanopore 500MB", 0.155, 0.319),
]

names = [r[0] for r in parse_fasta]
s_seq = [r[1] for r in parse_fasta]
t_load = [r[2] for r in parse_fasta]

fq_names = [r[0] for r in parse_fastq]
fq_seq = [r[1] for r in parse_fastq]
fq_load = [r[2] for r in parse_fastq]

# ── Plot A: FASTA fair parse ──────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 5.5))
x = np.arange(len(names)); w = 0.36
ax.bar(x - w/2, s_seq, w, label="SeqKit seq (parse-only)", color="#ff7f0e")
ax.bar(x + w/2, t_load, w, label="TracEon load (parse+index+arena)", color="#1f77b4")
ax.set_xticks(x); ax.set_xticklabels(names, rotation=45, ha="right")
ax.set_ylabel("Time (s)")
ax.set_title("FAIR parse: TracEon load (does MORE work: index+arena) vs SeqKit parse-only — 8/11 FASTA wins")
ax.legend()
fig.tight_layout(); fig.savefig("plotA_fair_parse_fasta.png", dpi=130); plt.close(fig)

# ── Plot B: FASTQ fair parse ──────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(3); w = 0.36
ax.bar(x - w/2, fq_seq, w, label="SeqKit seq (parse-only)", color="#ff7f0e")
ax.bar(x + w/2, fq_load, w, label="TracEon load (parse+index+arena)", color="#1f77b4")
for i, v in enumerate(fq_seq):
    ax.annotate(f"{v:.3f}", (i - w/2, v), ha="center", va="bottom", fontsize=8)
for i, v in enumerate(fq_load):
    ax.annotate(f"{v:.3f}", (i + w/2, v), ha="center", va="bottom", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(fq_names)
ax.set_ylabel("Time (s)")
ax.set_title("FASTQ: honest negative — SeqKit wins 1.5-2.1x (index insert = 35-52% of load)")
ax.legend()
fig.tight_layout(); fig.savefig("plotB_fair_parse_fastq.png", dpi=130); plt.close(fig)

# ── Plot C: End-to-end 5-run workflow (log scale) ─────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# WGS.gz
tools = ["SeqKit\nscan", "BioPython\nre-parse", "TracEon\ncache+lookup", "PyFastX\nindex+lookup"]
vals = [1.06, 6.56, 7.92, 182.28]
colors = ["#ff7f0e", "#d62728", "#1f77b4", "#2ca02c"]
ax = axes[0]
bars = ax.bar(tools, vals, color=colors)
ax.set_yscale("log")
for b, v in zip(bars, vals):
    ax.annotate(f"{v:.1f}s", (b.get_x()+b.get_width()/2, b.get_height()), ha="center", va="bottom", fontsize=8)
ax.set_title("WGS.gz short reads\n(5 runs × 50K lookups)")
ax.set_ylabel("Total time (s, log)")

# Yeast FASTA
tools = ["BioPython\nre-parse", "TracEon\ncache+lookup", "SeqKit\nfaidx", "PyFastX\nindex+lookup"]
vals = [1.10, 2.89, 4.63, 16.65]
colors = ["#d62728", "#1f77b4", "#ff7f0e", "#2ca02c"]
ax = axes[1]
bars = ax.bar(tools, vals, color=colors)
ax.set_yscale("log")
for b, v in zip(bars, vals):
    ax.annotate(f"{v:.2f}s", (b.get_x()+b.get_width()/2, b.get_height()), ha="center", va="bottom", fontsize=8)
ax.set_title("Yeast FASTA\n(5 runs × 1K lookups)")
ax.set_ylabel("Total time (s, log)")

# Nanopore
tools = ["SeqKit\nscan", "BioPython\nre-parse", "PyFastX\nindex+lookup", "TracEon\ncache+lookup"]
vals = [0.83, 6.61, 2.81, 45.32]
colors = ["#ff7f0e", "#d62728", "#2ca02c", "#1f77b4"]
ax = axes[2]
bars = ax.bar(tools, vals, color=colors)
ax.set_yscale("log")
for b, v in zip(bars, vals):
    ax.annotate(f"{v:.1f}s", (b.get_x()+b.get_width()/2, b.get_height()), ha="center", va="bottom", fontsize=8)
ax.set_title("Nanopore 500MB\n(5 runs × 50K lookups, LZ4HC save 43.3s)")
ax.set_ylabel("Total time (s, log)")

fig.suptitle("End-to-end repeated-use workflow — honest: LZ4HC save cost shows in Nanopore", fontsize=11)
fig.tight_layout(); fig.savefig("plotC_e2e_workflow.png", dpi=130); plt.close(fig)

# ── Plot D: Lookup per-query latency (server-style, no save cost) ─────
fig, ax = plt.subplots(figsize=(8, 5))
labels = ["TracEon\nin-memory\n(50K lookups)", "SeqKit scan\nper run\n(full file)", "PyFastX cold\ngz random\naccess", "seqkit faidx\nper-ID\n(linear)"]
vals = [0.004, 0.21, 36.4, 0.9]  # seconds for 50K lookups / per-run equivalent
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]
bars = ax.bar(labels, vals, color=colors)
ax.set_yscale("log")
for b, v in zip(bars, vals):
    ax.annotate(f"{v:.3f}s" if v < 1 else f"{v:.1f}s", (b.get_x()+b.get_width()/2, b.get_height()), ha="center", va="bottom", fontsize=8)
ax.set_ylabel("Time for 50K lookups (s, log)")
ax.set_title("Per-query latency: server-style deployment (load once, no save/restore)")
fig.tight_layout(); fig.savefig("plotD_latency.png", dpi=130); plt.close(fig)

print("Plots written:")
import os
for f in ["plotA_fair_parse_fasta.png","plotB_fair_parse_fastq.png","plotC_e2e_workflow.png","plotD_latency.png"]:
    print(" ", f, os.path.getsize(f)//1024, "KB")
