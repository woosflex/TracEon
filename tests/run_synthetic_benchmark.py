#!/usr/bin/env python3
"""
TracEon Synthetic GZIP Benchmark Suite
========================================
Generates synthetic FASTA.gz files at multiple sizes, runs the
chrono-instrumented benchmark binary 3 iterations per size, and
produces a JSON results file plus human-readable report.

Usage:
    python run_synthetic_benchmark.py [--sizes 10,50,100,500] [--iterations 3]

Output:
    benchmarks/v1.1.0_results.json   (structured timing data)
    stdout                            (human-readable report)
"""

import argparse
import gzip
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path
from statistics import mean, stdev

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK_BIN = PROJECT_ROOT / "tests" / "gzip_bench_test"
RESULTS_FILE  = PROJECT_ROOT / "benchmarks" / "v1.1.0_results.json"
TEMP_DIR      = Path("/tmp/traceon_bench")

# ── Helpers ────────────────────────────────────────────────────────────────

BASES = "ACGT"


def generate_fasta_gz(path: Path, target_mb: int, seq_len: int = 5000):
    """Generate a FASTA.gz file of approximately *target_mb* MB compressed size.

    Uses calibration to estimate compression ratio, then writes the
    file in one pass from scratch.
    """
    target_bytes = target_mb * 1024 * 1024
    seq = "".join(random.choices(BASES, k=seq_len))

    # ── Calibration ────────────────────────────────────────────────────────
    calib_path = path.with_suffix(".calib.tmp")
    N = 500
    uncomp = 0
    with gzip.open(str(calib_path), "wt", compresslevel=6) as f:
        for i in range(N):
            r = f">{i}\n{seq}\n"
            f.write(r)
            uncomp += len(r)
    comp = calib_path.stat().st_size
    calib_path.unlink()
    ratio = comp / uncomp if uncomp > 0 else 0.3

    # ── Generate full file ─────────────────────────────────────────────────
    bytes_per_rec = seq_len + 12
    target_recs = int(target_bytes / (bytes_per_rec * ratio)) + 2000

    # Write in medium batches to balance throughput vs memory
    BATCH = 5000
    with gzip.open(str(path), "wt", compresslevel=6) as f:
        for start in range(0, target_recs, BATCH):
            end = min(start + BATCH, target_recs)
            lines = [f">{i}\n{seq}\n" for i in range(start, end)]
            f.write("".join(lines))

    actual_mb = path.stat().st_size / (1024 * 1024)
    return actual_mb


def run_benchmark(file_path: Path) -> dict:
    """Execute the chrono-instrumented C++ benchmark once.

    Returns parsed JSON dict with timing fields.
    """
    if not BENCHMARK_BIN.exists():
        print(f"ERROR: benchmark binary not found at {BENCHMARK_BIN}", file=sys.stderr)
        sys.exit(1)

    try:
        result = subprocess.run(
            [str(BENCHMARK_BIN), str(file_path)],
            capture_output=True, text=True, timeout=300,
        )
        if result.returncode != 0:
            print(f"WARN: benchmark exited {result.returncode} for {file_path.name}")
            print(f"  stderr: {result.stderr.strip()}")
            return None

        # Parse the JSON output line
        for line in result.stdout.strip().splitlines():
            line = line.strip()
            if line.startswith("{"):
                return json.loads(line)
        print(f"WARN: no JSON output from benchmark for {file_path.name}")
        return None
    except subprocess.TimeoutExpired:
        print(f"ERROR: benchmark timed out on {file_path.name}")
        return None
    except json.JSONDecodeError as e:
        print(f"ERROR: JSON parse failure for {file_path.name}: {e}")
        print(f"  stdout: {result.stdout[:500]}")
        return None


def compute_stats(values):
    """Return (mean, std) tuple; std=0 if fewer than 2 samples."""
    if len(values) < 2:
        return (values[0] if values else 0.0), 0.0
    return mean(values), stdev(values)


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="TracEon Synthetic GZIP Benchmark Suite"
    )
    parser.add_argument(
        "--sizes", default="10,50,100,500",
        help="Comma-separated target compressed sizes in MB (default: 10,50,100,500)"
    )
    parser.add_argument(
        "--iterations", type=int, default=3,
        help="Number of iterations per file size (default: 3)"
    )
    parser.add_argument(
        "--seq-len", type=int, default=5000,
        help="Sequence length in bases (default: 5000)"
    )
    parser.add_argument(
        "--keep-files", action="store_true",
        help="Do not clean up generated FASTA.gz files"
    )
    parser.add_argument(
        "--build-only", action="store_true",
        help="Only compile the benchmark binary, do not run"
    )
    args = parser.parse_args()

    sizes = [int(s.strip()) for s in args.sizes.split(",")]
    iterations = args.iterations
    seq_len = args.seq_len

    # ── Build check ────────────────────────────────────────────────────────
    if not BENCHMARK_BIN.exists():
        print(f"Benchmark binary not found, compiling...")
        build_cmd = [
            "g++", "-std=c++20", "-O2", "-DNDEBUG",
            "-I", str(PROJECT_ROOT / "include"),
            "-I", str(PROJECT_ROOT / "third_party" / "lz4" / "lib"),
            "-I", str(PROJECT_ROOT / "third_party"),
            "-o", str(BENCHMARK_BIN),
            str(PROJECT_ROOT / "tests" / "gzip_benchmark_test.cpp"),
            "-L", str(PROJECT_ROOT / "build"),
            "-L", str(PROJECT_ROOT / "build" / "third_party" / "lz4" / "build" / "cmake"),
            "-ltraceon_core", "-lz", "-llz4",
        ]
        ret = subprocess.run(build_cmd, capture_output=True, text=True)
        if ret.returncode != 0:
            print(f"Build failed:\n{ret.stderr}", file=sys.stderr)
            sys.exit(1)
        print("Build successful.\n")

    if args.build_only:
        print("Benchmark binary ready at:", BENCHMARK_BIN)
        return

    # ── Prepare temp directory ─────────────────────────────────────────────
    TEMP_DIR.mkdir(parents=True, exist_ok=True)

    # ── Results accumulator ────────────────────────────────────────────────
    results = {"sizes": {}, "metadata": {}}
    metadata = {
        "tool": "TracEon",
        "version": "1.1.0",
        "benchmark_binary": str(BENCHMARK_BIN),
        "iterations": iterations,
        "seq_len": seq_len,
        "date": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": os.uname().nodename,
        "platform": f"{os.uname().sysname} {os.uname().machine}",
    }
    results["metadata"] = metadata

    generated_files = []

    print("=" * 72)
    print("  TracEon v1.1.0 — Synthetic GZIP Benchmark Suite")
    print("=" * 72)
    print(f"  Sizes:      {', '.join(f'{s}MB' for s in sizes)}")
    print(f"  Iterations: {iterations} per size")
    print(f"  Seq length: {seq_len} bp")
    print(f"  Binary:     {BENCHMARK_BIN}")
    print("=" * 72)

    overall_ok = True

    for target_mb in sizes:
        print(f"\n{'─' * 72}")
        print(f"  FILE SIZE: {target_mb} MB (compressed target)")
        print(f"{'─' * 72}")

        fname = f"bench_fasta_{target_mb}mb.fasta.gz"
        fpath = TEMP_DIR / fname

        # ── Generate ───────────────────────────────────────────────────────
        print(f"  Generating {fname} ...", end=" ", flush=True)
        t0 = time.time()
        try:
            actual_mb = generate_fasta_gz(fpath, target_mb, seq_len)
        except Exception as e:
            print(f" FAILED: {e}")
            overall_ok = False
            continue
        gen_time = time.time() - t0
        print(f"done ({actual_mb:.1f} MB compressed, {gen_time:.1f}s)")

        generated_files.append(fpath)

        # Container for raw iteration samples
        raw = {
            "decompress_s": [],
            "parse_s": [],
            "total_s": [],
            "records": [],
            "raw_bytes": [],
            "compressed_bytes": fpath.stat().st_size,
        }

        # ── Benchmark loop ─────────────────────────────────────────────────
        failed_iters = 0
        for i in range(iterations):
            print(f"    Iteration {i+1}/{iterations} ...", end=" ", flush=True)
            result = run_benchmark(fpath)
            if result is None:
                print("FAILED")
                failed_iters += 1
                continue

            raw["decompress_s"].append(result["decompress_s"])
            raw["parse_s"].append(result["parse_s"])
            raw["total_s"].append(result["total_s"])
            raw["records"].append(result["records"])
            raw["raw_bytes"].append(result["raw_bytes"])

            print(f"decomp={result['decompress_s']:.4f}s  "
                  f"parse={result['parse_s']:.4f}s  "
                  f"total={result['total_s']:.4f}s  "
                  f"recs={result['records']}")

        if failed_iters > 0:
            print(f"  WARNING: {failed_iters}/{iterations} iterations failed")
            if failed_iters == iterations:
                overall_ok = False
                continue

        # ── Compute summary statistics ─────────────────────────────────────
        def summarize(values):
            """Return dict with mean, std, min, max."""
            m, s = compute_stats(values)
            return {"mean": round(m, 6), "std": round(s, 6),
                    "min": round(min(values), 6), "max": round(max(values), 6)}

        summary = {
            "decompress_s": summarize(raw["decompress_s"]),
            "parse_s":      summarize(raw["parse_s"]),
            "total_s":      summarize(raw["total_s"]),
            "records":      summarize(raw["records"]),
            "raw_bytes":    summarize(raw["raw_bytes"]),
            "compressed_bytes": fpath.stat().st_size,
            "compressed_mb": round(actual_mb, 2),
        }
        # Store both raw samples and summary
        raw["compressed_bytes"] = fpath.stat().st_size
        results["sizes"][str(target_mb)] = {
            "samples": raw,
            "summary": summary,
        }

        # ── Print summary line ─────────────────────────────────────────────
        d = summary["decompress_s"]
        p = summary["parse_s"]
        t = summary["total_s"]
        print(f"\n  -> {target_mb}MB:  decompress {d['mean']:.4f} ± {d['std']:.4f}s  "
              f"parse {p['mean']:.4f} ± {p['std']:.4f}s  "
              f"total {t['mean']:.4f} ± {t['std']:.4f}s  "
              f"({int(summary['records']['mean'])} records)")

    # ── Cleanup generated files ────────────────────────────────────────────
    if not args.keep_files:
        for f in generated_files:
            try:
                f.unlink()
            except Exception:
                pass

    # ── Write results JSON ─────────────────────────────────────────────────
    RESULTS_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(str(RESULTS_FILE), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved to: {RESULTS_FILE}")

    # ── Final Summary Table ────────────────────────────────────────────────
    print("\n")
    print("=" * 104)
    print("  COMPONENT TIMING BREAKDOWN  (mean ± std, seconds)")
    print("=" * 104)
    header = f"  {'Size':>8} | {'Compressed':>10} | {'Raw Bytes':>10} | "
    header += f"{'Records':>8} | {'Decompress':>14} | {'Parse':>14} | {'Total':>14}"
    print(header)
    print("  " + "-" * 100)

    for target_mb in sizes:
        key = str(target_mb)
        if key not in results["sizes"]:
            continue
        s = results["sizes"][key]["summary"]
        comp_mb = s["compressed_mb"]
        raw_bytes_mb = s["raw_bytes"]["mean"] / (1024 * 1024)
        recs = int(s["records"]["mean"])
        d = s["decompress_s"]
        p = s["parse_s"]
        t = s["total_s"]
        print(f"  {target_mb:>4}MB | {comp_mb:>8.1f}MB | {raw_bytes_mb:>8.1f}MB | "
              f"{recs:>8} | "
              f"{d['mean']:.4f}±{d['std']:.4f} | "
              f"{p['mean']:.4f}±{p['std']:.4f} | "
              f"{t['mean']:.4f}±{t['std']:.4f}")

    print("=" * 104)

    # ── Throughput table ──────────────────────────────────────────────────
    print("\n")
    print("=" * 72)
    print("  THROUGHPUT (MB/s raw bytes processed)")
    print("=" * 72)
    print(f"  {'Size':>8} | {'Decompress':>14} | {'Parse':>14} | {'Total':>14}")
    print("  " + "-" * 56)

    for target_mb in sizes:
        key = str(target_mb)
        if key not in results["sizes"]:
            continue
        s = results["sizes"][key]["summary"]
        raw_mb = s["raw_bytes"]["mean"] / (1024 * 1024)
        d = s["decompress_s"]
        p = s["parse_s"]
        t = s["total_s"]

        def mbps(mbytes, time_dict):
            t_mean = time_dict["mean"]
            return mbytes / t_mean if t_mean > 0 else 0

        d_mbps = mbps(raw_mb, d)
        p_mbps = mbps(raw_mb, p)
        t_mbps = mbps(raw_mb, t)
        print(f"  {target_mb:>4}MB | {d_mbps:>8.1f} MB/s     | "
              f"{p_mbps:>8.1f} MB/s     | {t_mbps:>8.1f} MB/s")

    print("=" * 72)

    # ── Check for gnuplot/matplotlib ──────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        print("\n  Generating comparison plot (benchmark_comparison.png) ...", end=" ", flush=True)

        fig, axes = plt.subplots(1, 3, figsize=(16, 5))

        sizes_plot = []
        decompress_means = []
        parse_means = []
        total_means = []
        decompress_stds = []
        parse_stds = []
        total_stds = []

        for target_mb in sizes:
            key = str(target_mb)
            if key not in results["sizes"]:
                continue
            s = results["sizes"][key]["summary"]
            sizes_plot.append(target_mb)
            decompress_means.append(s["decompress_s"]["mean"])
            parse_means.append(s["parse_s"]["mean"])
            total_means.append(s["total_s"]["mean"])
            decompress_stds.append(s["decompress_s"]["std"])
            parse_stds.append(s["parse_s"]["std"])
            total_stds.append(s["total_s"]["std"])

        # Plot 1: Component breakdown (stacked bars)
        x = range(len(sizes_plot))
        width = 0.6
        ax = axes[0]
        ax.bar(x, decompress_means, width, yerr=decompress_stds, capsize=4,
               label="Decompress", color="#2ecc71", alpha=0.85)
        ax.bar(x, parse_means, width, yerr=parse_stds, capsize=4,
               bottom=decompress_means, label="Parse", color="#3498db", alpha=0.85)
        ax.set_xlabel("Compressed File Size (MB)")
        ax.set_ylabel("Wall Time (s)")
        ax.set_title("Load Time Breakdown")
        ax.set_xticks(list(x))
        ax.set_xticklabels([f"{s}MB" for s in sizes_plot])
        ax.legend()
        ax.grid(axis="y", alpha=0.3)

        # Plot 2: Total load time with error bars
        ax = axes[1]
        ax.errorbar(sizes_plot, total_means, yerr=total_stds,
                     fmt="-o", color="#e74c3c", capsize=5, capthick=2, linewidth=2)
        ax.set_xlabel("Compressed File Size (MB)")
        ax.set_ylabel("Total Load Time (s)")
        ax.set_title("Total Load Time vs File Size")
        ax.grid(True, alpha=0.3)

        # Plot 3: Throughput comparison
        ax = axes[2]
        raw_mbs = [results["sizes"][str(s)]["summary"]["raw_bytes"]["mean"] / (1024*1024)
                   for s in sizes_plot if str(s) in results["sizes"]]
        decomp_mbps = [rm / d_m if d_m > 0 else 0
                       for rm, d_m in zip(raw_mbs, decompress_means)]
        parse_mbps  = [rm / p_m if p_m > 0 else 0
                       for rm, p_m in zip(raw_mbs, parse_means)]
        total_mbps  = [rm / t_m if t_m > 0 else 0
                       for rm, t_m in zip(raw_mbs, total_means)]

        ax.plot(sizes_plot, decomp_mbps, "-s", label="Decompress", color="#2ecc71", linewidth=2)
        ax.plot(sizes_plot, parse_mbps, "-^", label="Parse", color="#3498db", linewidth=2)
        ax.plot(sizes_plot, total_mbps, "-o", label="Total", color="#e74c3c", linewidth=2)
        ax.set_xlabel("Compressed File Size (MB)")
        ax.set_ylabel("Throughput (MB/s)")
        ax.set_title("Processing Throughput")
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plot_path = PROJECT_ROOT / "benchmarks" / "benchmark_comparison.png"
        plt.savefig(str(plot_path), dpi=150, bbox_inches="tight")
        plt.close()
        print(f"done -> {plot_path}")

    except ImportError:
        print("\n  [SKIP] matplotlib not available — skipping plot generation")
    except Exception as e:
        print(f"\n  [WARN] plot generation failed: {e}")

    print("\n" + "=" * 72)
    print("  BENCHMARK COMPLETE")
    print("=" * 72)

    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
