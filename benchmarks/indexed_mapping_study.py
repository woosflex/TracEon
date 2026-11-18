#!/usr/bin/env python3
"""Installed, snapshot-pinned TracEon versus samtools indexed read retrieval.

One resident server serves independent subset-remapping workers. All query/map
processes share the same explicit CPU budget. This is a single-end workflow;
qualities are preserved for FASTQ. No mapper algorithm is replaced.
"""
import argparse
import json
import itertools
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess

from mapping_baseline import positive, run_measured, sam_summary, sha256
from shared_mapping_pilot import run_group, start_server
from subset_mapping_pilot import fasta_records


def sequence_records(path, fastq):
    if not fastq:
        for key, sequence in fasta_records(path):
            yield key, sequence, ""
        return
    with Path(path).open() as stream:
        while header := stream.readline():
            sequence, plus, quality = (stream.readline() for _ in range(3))
            if not header.startswith("@") or not plus.startswith("+") or not quality:
                raise ValueError(f"expected four-line FASTQ in {path}")
            names = header[1:].split()
            sequence, quality = sequence.rstrip("\r\n"), quality.rstrip("\r\n")
            if not names or len(sequence) != len(quality):
                raise ValueError(f"invalid FASTQ ID or payload lengths in {path}")
            yield names[0], sequence, quality


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("traceon", "samtools", "minimap2", "reference", "reads", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--fastq", action="store_true")
    parser.add_argument("--preset", required=True)
    parser.add_argument("--cpus", type=int, nargs="+", required=True,
                        help="same Linux CPU affinity pool for all arms, including server")
    parser.add_argument("--controller-cpu", type=int, help="keep validation/resource sampling outside the worker CPU pool")
    parser.add_argument("--processes", type=positive, nargs="+", default=[1, 2, 4])
    parser.add_argument("--counts", type=positive, nargs="+", default=[1000])
    parser.add_argument("--repetitions", type=positive, default=6)
    parser.add_argument("--batch-size", type=positive, help="add a batched arm (2–256 IDs per request)")
    args = parser.parse_args()
    if args.batch_size is not None and not 2 <= args.batch_size <= 256:
        parser.error("batch size must be 2–256")
    arms = ("indexed", "shared", "batched") if args.batch_size else ("indexed", "shared")
    orders = list(itertools.permutations(arms))
    if args.repetitions % len(orders):
        parser.error(f"repetitions must be a multiple of {len(orders)} for balanced arm order")
    allowed = os.sched_getaffinity(0)
    if not set(args.cpus) <= allowed or len(set(args.cpus)) != len(args.cpus):
        parser.error("CPUs must be distinct and available in the current affinity")
    if args.controller_cpu is not None:
        if args.controller_cpu not in allowed or args.controller_cpu in args.cpus:
            parser.error("controller CPU must be available and outside the worker pool")
        os.sched_setaffinity(0, {args.controller_cpu})
    paths = {k: v.resolve(strict=True) for k, v in vars(args).items() if isinstance(v, Path) and k != "output"}
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    records = {}
    for key, sequence, quality in sequence_records(paths["reads"], args.fastq):
        if key in records:
            raise ValueError("duplicate IDs are not supported by the indexed baseline")
        records[key] = (sequence, quality)
    if not records or max(args.counts) > len(records):
        raise ValueError("not enough distinct source reads for the requested subsets")
    names = list(records)
    prefix = ["taskset", "-c", ",".join(map(str, args.cpus))]
    manifest = {"platform": platform.platform(), "arguments": vars(args) | {k: str(v) for k, v in paths.items()} | {"output": str(output)},
                "arm_orders": orders, "source_reads": len(records), "inputs": {k: {"path": str(v), "sha256": sha256(v)} for k, v in paths.items()},
                "harness_sha256": {p.name: sha256(p) for p in (Path(__file__), Path(__file__).with_name("shared_mapping_pilot.py"),
                                   Path(__file__).with_name("mapping_baseline.py"), Path(__file__).with_name("subset_mapping_pilot.py"))},
                "versions": {"minimap2": subprocess.check_output([str(paths["minimap2"]), "--version"], text=True).strip(),
                             "samtools": subprocess.check_output([str(paths["samtools"]), "--version"], text=True)},
                "cpu_topology": subprocess.check_output(["lscpu", "-e=CPU,CORE,SOCKET,ONLINE"], text=True),
                "cpu_info": Path("/proc/cpuinfo").read_text(), "memory_info": Path("/proc/meminfo").read_text(),
                "source_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                "working_tree": subprocess.check_output(["git", "status", "--short"], text=True),
                "cache_policy": "warm OS page cache: inputs scanned for validation; one unmeasured warmup per cohort; no cache dropping",
                "memory_scope": "5ms sampled aggregate process PSS; shared includes server, indexed excludes idle server; controller excluded",
                "cpu_policy": "identical affinity pool for every arm and server; not exclusive/reserved CPUs; frequency/governor unchanged"}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    cache, fai, index = output / "reads.traceon", output / "reads.fai", output / "reference.mmi"
    cache_build = run_measured(prefix + [str(paths["traceon"]), "build", str(paths["reads"]), str(cache)], output, "cache-build")
    operation = "fqidx" if args.fastq else "faidx"
    fai_build = run_measured(prefix + [str(paths["samtools"]), operation, "--fai-idx", str(fai), str(paths["reads"])], output, "fai-build")
    reference_build = run_measured(prefix + [str(paths["minimap2"]), "-x", args.preset, "-t", "1",
                                            "-d", str(index), str(paths["reference"])], output, "reference-index")
    server, port, setup = start_server(paths["traceon"], cache, output, supported_cli=True, cpus=args.cpus)
    samples = []
    try:
        if setup["records"] != len(records) or setup["cpu_affinity"] != sorted(args.cpus):
            raise ValueError("server dataset count or affinity mismatch")
        for processes in args.processes:
            for count in args.counts:
                for repetition in range(args.repetitions + 1):
                    jobs = []
                    for worker in range(processes):
                        ids = random.Random(42 + processes * 100000 + count * 100 + repetition * 10 + worker).sample(names, count)
                        path = output / f"ids-p{processes}-n{count}-r{repetition}-w{worker}.txt"
                        path.write_text("".join(key + "\n" for key in ids))
                        jobs.append((path, [(key, *records[key]) for key in ids]))
                    canonical = None
                    for arm in orders[(repetition - 1) % len(orders)]:
                        label = f"{arm}-p{processes}-n{count}-r{repetition}"
                        commands = []
                        for ids, _ in jobs:
                            if arm == "indexed":
                                command = [str(paths["samtools"]), operation, "--fai-idx", str(fai), "-n", "0", "-r", str(ids), str(paths["reads"])]
                            else:
                                command = [str(paths["traceon"]), "query", "--port", str(port),
                                           "--expect-snapshot", setup["snapshot"], "--ids", str(ids)]
                                if arm == "batched":
                                    command.extend(["--batch-size", str(args.batch_size)])
                                if args.fastq:
                                    command.append("--fastq")
                            commands.append(prefix + command)
                        resident = server.pid if arm != "indexed" else None
                        retrieval = run_group(commands, output, label + "-retrieve", resident)
                        for worker, (_, expected) in enumerate(jobs):
                            actual = list(sequence_records(output / f"{label}-retrieve-{worker}.stdout", args.fastq))
                            if actual != expected:
                                raise ValueError(f"{label}: payload/order mismatch for worker {worker}")
                        commands = [prefix + [str(paths["minimap2"]), "-a", "-x", args.preset, "-t", "1",
                                             str(index), str(output / f"{label}-retrieve-{worker}.stdout")] for worker in range(processes)]
                        mapping = run_group(commands, output, label + "-map", resident)
                        summaries = [sam_summary(output / f"{label}-map-{worker}.stdout") for worker in range(processes)]
                        if any(s["primary_records"] != count for s in summaries):
                            raise ValueError(f"{label}: primary read count mismatch")
                        if canonical is not None and summaries != canonical:
                            raise ValueError(f"{label}: SAM alignment/quality mismatch between arms")
                        canonical = summaries
                        if server.poll() is not None:
                            raise RuntimeError("server exited during workflow")
                        if repetition:
                            samples.append({"arm": arm, "processes": processes, "count": count, "repetition": repetition,
                                            "retrieval": retrieval, "mapping": mapping, "sam": summaries,
                                            "workflow_seconds": retrieval["wall_seconds"] + mapping["wall_seconds"]})
    finally:
        server.terminate()
        try:
            server.wait(timeout=10)
        except subprocess.TimeoutExpired:
            server.kill()
            server.wait()
        server.stdout.close()
    if server.returncode:
        raise RuntimeError(f"unclean server shutdown: {server.returncode}")
    aggregates = []
    for processes in args.processes:
        for count in args.counts:
            pair = {}
            for arm in arms:
                selected = [s for s in samples if (s["arm"], s["processes"], s["count"]) == (arm, processes, count)]
                row = {"arm": arm, "processes": processes, "count_per_worker": count, "samples": len(selected),
                       "median_retrieval_seconds": statistics.median(s["retrieval"]["wall_seconds"] for s in selected),
                       "median_workflow_seconds": statistics.median(s["workflow_seconds"] for s in selected),
                       "peak_retrieval_pss_kib": max(s["retrieval"]["sampled_peak_aggregate_pss_kib"] for s in selected)}
                row["queries_per_second"] = processes * count / row["median_retrieval_seconds"]
                pair[arm] = row
                aggregates.append(row)
            setup_difference = cache_build["wall_seconds"] + setup["startup_seconds"] - fai_build["wall_seconds"]
            for arm in arms[1:]:
                saving = pair["indexed"]["median_retrieval_seconds"] - pair[arm]["median_retrieval_seconds"]
                pair[arm]["reuse_break_even_groups_estimate"] = max(0, setup_difference) / saving if saving > 0 else None
    report = {"status": "validated installed shared-cache mapping workflow", "cache_build": cache_build,
              "fai_build": fai_build, "reference_build": reference_build, "server_setup": setup,
              "server_exit": server.returncode, "cache_sha256": sha256(cache), "fai_sha256": sha256(fai),
              "samples": samples, "aggregates": aggregates}
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(aggregates, indent=2))


if __name__ == "__main__":
    main()
