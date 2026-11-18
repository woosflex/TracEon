#!/usr/bin/env python3
"""Profile an existing validated study's resident cache without client restores.

Request times exclude handshake/warmup and include CRC plus releasing the owned
result. Native per-process windows start independently; cohort wall includes
startup and timing-file output. Server CPU is sampled from Linux /proc ticks.
"""
import argparse
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess

from mapping_baseline import positive, run_measured, sha256
from shared_mapping_pilot import run_group, start_server


def quantile(values, fraction):
    ordered = sorted(values)
    return ordered[int(fraction * (len(ordered) - 1))]


def fields(text, prefix):
    lines = [line for line in text.splitlines() if line.startswith(prefix + "\t")]
    if len(lines) != 1:
        raise ValueError(f"expected one {prefix} result")
    return dict(field.split("=", 1) for field in lines[0].split("\t")[1:])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("study", "traceon", "benchmark", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--cpus", type=int, nargs="+", required=True)
    parser.add_argument("--controller-cpu", type=int, required=True)
    parser.add_argument("--processes", type=positive, nargs="+", default=[1, 2, 4])
    parser.add_argument("--requests", type=positive, default=20000)
    parser.add_argument("--repetitions", type=positive, default=6)
    parser.add_argument("--batch-size", type=positive, default=1)
    args = parser.parse_args()
    if args.batch_size > 256:
        parser.error("batch size must be 1–256")
    allowed = os.sched_getaffinity(0)
    if not set(args.cpus) <= allowed or args.controller_cpu not in allowed or args.controller_cpu in args.cpus:
        parser.error("worker CPUs and a distinct controller CPU must be available")
    os.sched_setaffinity(0, {args.controller_cpu})
    study, binary, benchmark = (p.resolve(strict=True) for p in (args.study, args.traceon, args.benchmark))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    prior = json.loads((study / "summary.json").read_text())
    inputs = json.loads((study / "manifest.json").read_text())
    fastq = inputs["arguments"]["fastq"]
    cache, fai = study / "reads.traceon", study / "reads.fai"
    if sha256(cache) != prior["cache_sha256"] or sha256(fai) != prior["fai_sha256"]:
        raise ValueError("the validated study cache/index changed")
    entries = [line.split("\t") for line in fai.read_text().splitlines()]
    if len(entries) != inputs["source_reads"]:
        raise ValueError("index record count mismatch")
    ids = output / "ids.txt"
    ids.write_text("".join(row[0] + "\n" for row in entries))
    lengths = [int(row[1]) * (2 if fastq else 1) for row in entries]
    manifest = {"arguments": vars(args) | {"study": str(study), "traceon": str(binary), "benchmark": str(benchmark), "output": str(output)},
                "platform": platform.platform(), "binary_sha256": sha256(binary), "benchmark_sha256": sha256(benchmark),
                "source_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                "working_tree": subprocess.check_output(["git", "status", "--short"], text=True),
                "harness_sha256": {p.name: sha256(p) for p in (Path(__file__), Path(__file__).with_name("shared_mapping_pilot.py"),
                    Path(__file__).with_name("mapping_baseline.py"), Path(__file__).with_name("subset_mapping_pilot.py"))},
                "benchmark_source_sha256": sha256(Path(__file__).with_name("remote_bench.cpp")),
                "validated_study_manifest_sha256": sha256(study / "manifest.json"),
                "cache_sha256": sha256(cache), "snapshot": prior["server_setup"]["snapshot"], "fastq": fastq,
                "payload_bytes": {"min": min(lengths), "median": statistics.median(lengths), "max": max(lengths)},
                "latency_unit": "request (one batch when batch_size > 1)",
                "records_per_request": args.batch_size, "warmup_requests_per_process": 1000, "server_cpu_tick_seconds": 1 / os.sysconf("SC_CLK_TCK"),
                "policy": "same affinity pool for server/clients; independent native timing windows; warm reused cache; no exclusive cores/frequency control"}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    prefix = ["taskset", "-c", ",".join(map(str, args.cpus))]
    server, port, setup = start_server(binary, cache, output, supported_cli=True, cpus=args.cpus)
    samples = []
    try:
        if setup["snapshot"] != manifest["snapshot"]:
            raise ValueError("restored snapshot identity changed")
        for processes in args.processes:
            for repetition in range(args.repetitions):
                label = f"p{processes}-r{repetition}"
                commands, files = [], []
                for worker in range(processes):
                    raw = output / f"{label}-w{worker}-latency-us.txt"
                    files.append(raw)
                    command = prefix + [str(benchmark), "--mode", "remote", "--ids", str(ids),
                        "--expect-snapshot", setup["snapshot"], "--port", str(port), "--threads", "1",
                        "--count", str(args.requests), "--warmup", "1000", "--seed", str(42 + repetition * 10 + worker),
                        "--latencies", str(raw), "--batch-size", str(args.batch_size)]
                    if fastq:
                        command.append("--fastq")
                    commands.append(command)
                group = run_group(commands, output, label, server.pid)
                latencies, native, payload_bytes = [], [], 0
                for worker, raw in enumerate(files):
                    values = [float(value) for value in raw.read_text().splitlines()]
                    if len(values) != args.requests or any(not math.isfinite(v) or v < 0 for v in values):
                        raise ValueError("invalid request timing sample count/value")
                    latencies.extend(values)
                    stdout = (output / f"{label}-{worker}.stdout").read_text()
                    row = fields(stdout, "TRACEON_BENCH")
                    if int(row["count"]) != args.requests or int(row["threads"]) != 1:
                        raise ValueError("native benchmark request/thread count mismatch")
                    request = fields(stdout, "TRACEON_REQUESTS")
                    if int(request["batch_size"]) != args.batch_size or request["latency_unit"] != "request":
                        raise ValueError("native request units do not match the profile")
                    native.append(row)
                    payload_bytes += int(fields(stdout, "TRACEON_PAYLOAD")["bytes"])
                samples.append({"processes": processes, "repetition": repetition, "group": group,
                    "native": native, "requests": len(latencies), "records": len(latencies) * args.batch_size,
                    "mean_payload_bytes": payload_bytes / len(latencies),
                    "mean_record_bytes": payload_bytes / (len(latencies) * args.batch_size),
                    "median_us": statistics.median(latencies), "p95_us": quantile(latencies, .95), "p99_us": quantile(latencies, .99),
                    "native_windows_max_seconds": max(float(row["wall_s"]) for row in native)})
    finally:
        server.terminate()
        try:
            server.wait(timeout=10)
        except subprocess.TimeoutExpired:
            server.kill()
            server.wait()
        server.stdout.close()
    if server.returncode:
        raise RuntimeError("server failed to exit cleanly")
    codecs = []
    for size in (0, 100, 200, 1000, 10000, 100000):
        command = prefix + [str(benchmark), "--mode", "codec", "--bytes", str(size), "--count", "20000", "--warmup", "1000",
                            "--latencies", str(output / f"codec-{size}-us.tsv")]
        measurement = run_measured(command, output, f"codec-{size}")
        codecs.append({"measurement": measurement,
                       "medians": fields((output / f"codec-{size}.stdout").read_text(), "TRACEON_CODEC")})
    aggregate = []
    for processes in args.processes:
        selected = [s for s in samples if s["processes"] == processes]
        aggregate.append({"processes": processes, "repetitions": len(selected),
            "requests_per_repetition": processes * args.requests,
            **{f"median_{metric}": statistics.median(s[metric] for s in selected) for metric in ("median_us", "p95_us", "p99_us", "mean_payload_bytes")},
            "median_server_cpu_seconds_per_group_including_warmup": statistics.median(s["group"]["server_cpu_seconds"] for s in selected),
            "peak_aggregate_pss_kib": max(s["group"]["sampled_peak_aggregate_pss_kib"] for s in selected)})
    report = {"status": "validated request profile", "batch_size": args.batch_size, "setup": setup, "server_exit": server.returncode,
              "samples": samples, "aggregate": aggregate, "codecs": codecs}
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(aggregate, indent=2))


if __name__ == "__main__":
    main()
