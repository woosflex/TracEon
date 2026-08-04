#!/usr/bin/env python3
"""
Quick benchmark for TracEon v1.4.0 — generates synthetic data and measures core operations.
Designed to complete in <5 minutes with meaningful performance metrics.
"""

import os
import sys
import time
import subprocess
import random
import gzip
from pathlib import Path

TRACEON_BINARY = "../build/traceon_driver"
SIZES_MB = [10, 100]  # Quick benchmark: just 10MB and 100MB
SCENARIOS = [
    ("WGS", "fastq", 150, False),
    ("WGS_GZ", "fastq", 150, True),
    ("RefGenome", "fasta", 1_000_000, False),
    ("RefGenome_GZ", "fasta", 1_000_000, True),
]

def generate_dataset(fname, fmt, target_mb, avg_len, compress=False):
    """Generate synthetic FASTA or FASTQ file."""
    bytes_per_rec = (avg_len * 2 + 50) if "fastq" in fmt else (avg_len + 20)
    num_seqs = max(1, int((target_mb * 1024 * 1024) / bytes_per_rec))

    buf_size = 10_000_000
    buf = "".join(random.choices("ACGT", k=buf_size))

    print(f"  Generating {fname} ({num_seqs:,} records, ~{avg_len}bp)...", end=" ", flush=True)
    t0 = time.time()

    if compress:
        f = gzip.open(fname, "wt")
    else:
        f = open(fname, "w")

    with f:
        for i in range(num_seqs):
            s_len = avg_len
            if s_len > buf_size:
                seq = "".join(random.choices("ACGT", k=s_len))
            else:
                start = random.randint(0, buf_size - s_len)
                seq = buf[start:start+s_len]

            if "fastq" in fmt:
                f.write(f"@r{i}\n{seq}\n+\n{'I'*s_len}\n")
            else:
                f.write(f">c{i}\n{seq}\n")

    gen_time = time.time() - t0
    file_size = os.path.getsize(fname) / (1024 * 1024)
    print(f"✓ ({gen_time:.2f}s, {file_size:.1f} MB)")
    return num_seqs

def run_benchmark(label, fmt, size_mb, avg_len, compress):
    """Run benchmark for a single scenario."""
    if compress:
        fname = f"bench_{label}_{size_mb}MB.{fmt}.gz"
    else:
        fname = f"bench_{label}_{size_mb}MB.{fmt}"

    bin_file = fname.replace(".gz", "").replace(".fastq", ".bin").replace(".fasta", ".bin")

    n_seqs = generate_dataset(fname, fmt, size_mb, avg_len, compress)
    file_size = os.path.getsize(fname) / (1024 * 1024)

    results = {}

    # Load
    print(f"  Loading {fname}...", end=" ", flush=True)
    t0 = time.time()
    ret = subprocess.run(
        [TRACEON_BINARY, "load", fname],
        capture_output=True, text=True
    )
    load_time = time.time() - t0
    print(f"✓ ({load_time:.3f}s)")
    results['load_time'] = load_time
    results['records'] = n_seqs
    results['file_size'] = file_size

    # Save
    print(f"  Saving binary cache...", end=" ", flush=True)
    t0 = time.time()
    ret = subprocess.run(
        [TRACEON_BINARY, "save", fname, bin_file],
        capture_output=True, text=True
    )
    save_time = time.time() - t0
    print(f"✓ ({save_time:.3f}s)")
    results['save_time'] = save_time

    # Check cache size
    if os.path.exists(bin_file):
        cache_size = os.path.getsize(bin_file) / (1024 * 1024)
        results['cache_size'] = cache_size
        results['compression_ratio'] = file_size / cache_size if cache_size > 0 else 0

        # Restore
        print(f"  Restoring from cache...", end=" ", flush=True)
        t0 = time.time()
        ret = subprocess.run(
            [TRACEON_BINARY, "restore", bin_file],
            capture_output=True, text=True
        )
        restore_time = time.time() - t0
        print(f"✓ ({restore_time:.3f}s)")
        results['restore_time'] = restore_time
    else:
        print(f"  ⚠ Cache file not created")
        results['cache_size'] = 0
        results['restore_time'] = 0

    # Cleanup
    if os.path.exists(fname):
        os.remove(fname)
    if os.path.exists(bin_file):
        os.remove(bin_file)

    return results

def main():
    if not os.path.exists(TRACEON_BINARY):
        print("ERROR: TracEon binary not found. Run: cmake --build ../build -j")
        sys.exit(1)

    print("\n" + "="*80)
    print("TracEon v1.4.0 Quick Benchmark")
    print("="*80 + "\n")

    all_results = []

    for size_mb in SIZES_MB:
        print(f"Testing {size_mb} MB datasets:\n")

        for label, fmt, avg_len, compress in SCENARIOS:
            full_label = f"{label}_{size_mb}MB"
            print(f"▶ {full_label}")

            try:
                results = run_benchmark(label, fmt, size_mb, avg_len, compress)
                results['scenario'] = full_label
                all_results.append(results)
            except Exception as e:
                print(f"  ✗ Failed: {e}")

            print()

    # Summary table
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80 + "\n")

    print(f"{'Scenario':<30} | {'Load (s)':<8} | {'Save (s)':<8} | {'Restore (s)':<10} | {'Cache (MB)':<10} | {'Ratio':<6}")
    print("-" * 100)

    for r in all_results:
        scenario = r['scenario']
        load = f"{r['load_time']:.3f}"
        save = f"{r['save_time']:.3f}"
        restore = f"{r['restore_time']:.3f}"
        cache = f"{r['cache_size']:.1f}" if r['cache_size'] > 0 else "—"
        if r['cache_size'] > 0 and 'compression_ratio' in r:
            ratio = f"{r['compression_ratio']:.1f}x"
        else:
            ratio = "—"
        print(f"{scenario:<30} | {load:<8} | {save:<8} | {restore:<10} | {cache:<10} | {ratio:<6}")

    print("\n✅ Quick benchmark complete!")
    print("\nNote: This is a quick benchmark on 2 sizes × 4 scenarios.")
    print("For comprehensive results, see benchmark_runner.py (full matrix: 4 sizes × 8 scenarios)")

if __name__ == "__main__":
    main()
