#!/usr/bin/env python3
"""One resident TracEon server, many independent querying processes (Linux).

Compare shared RPC, per-process restore, and seqtk scans on identical changing
FASTA subsets. Report concurrent stage wall time and sampled aggregate PSS.
"""
import argparse
from contextlib import ExitStack
import json
import os
from pathlib import Path
import platform
import random
import re
import selectors
import statistics
import subprocess
import threading
import time

from mapping_baseline import positive, run_measured, sam_summary, sha256
from subset_mapping_pilot import fasta_records


def pss_kib(pid):
    try:
        for line in Path(f"/proc/{pid}/smaps_rollup").read_text().splitlines():
            if line.startswith("Pss:"):
                return int(line.split()[1])
    except (FileNotFoundError, ProcessLookupError):
        return 0  # Process exited between sampling and reading /proc.
    raise RuntimeError(f"PSS unavailable for PID {pid}")


def process_cpu_seconds(pid):
    fields = Path(f"/proc/{pid}/stat").read_text().rsplit(") ", 1)[1].split()
    return (int(fields[11]) + int(fields[12])) / os.sysconf("SC_CLK_TCK")


def run_group(commands, output, label, server_pid=None):
    """Block in wait4 for timing; sample aggregate PSS independently at 5ms.

    Wall time includes process startup but not sampler shutdown or log writes.
    Sampling can miss brief memory peaks; it must not quantize elapsed time.
    """
    children, metrics, peaks, errors = [], [], [], []
    done = threading.Event()
    def sample():
        try:
            while not done.is_set():
                peaks.append((pss_kib(server_pid) if server_pid else 0) +
                             sum(pss_kib(child.pid) for child in children))
                done.wait(0.005)
        except BaseException as error:
            errors.append(error)
    with ExitStack() as files:
        server_cpu = process_cpu_seconds(server_pid) if server_pid else None
        started = time.perf_counter()
        monitor = None
        try:
            for i, command in enumerate(commands):
                stdout = files.enter_context((output / f"{label}-{i}.stdout").open("wb"))
                stderr = files.enter_context((output / f"{label}-{i}.stderr").open("wb"))
                children.append(subprocess.Popen(command, stdout=stdout, stderr=stderr))
            monitor = threading.Thread(target=sample)
            monitor.start()
            for child, command in zip(children, commands):
                _, status, usage = os.wait4(child.pid, 0)
                child.returncode = os.waitstatus_to_exitcode(status)
                metrics.append({"returncode": child.returncode, "user_seconds": usage.ru_utime,
                                "system_seconds": usage.ru_stime, "peak_rss_kib": usage.ru_maxrss})
                if child.returncode:
                    raise subprocess.CalledProcessError(child.returncode, command)
            elapsed = time.perf_counter() - started
            cpu_delta = process_cpu_seconds(server_pid) - server_cpu if server_pid else None
            done.set()
            monitor.join()
            if errors:
                raise errors[0]
            result = {"label": label, "commands": commands, "wall_seconds": elapsed,
                      "sampled_peak_aggregate_pss_kib": max(peaks, default=0), "memory_samples": len(peaks),
                      "includes_server": server_pid is not None, "server_cpu_seconds": cpu_delta,
                      "children": metrics}
            with (output / "measurements.jsonl").open("a") as stream:
                stream.write(json.dumps(result) + "\n")
            return result
        finally:
            done.set()
            if monitor:
                monitor.join()
            for child in children:
                if child.returncode is None:
                    child.kill()
                    child.wait()


def start_server(binary, cache, output, *, supported_cli=False, extra_args=(), cpus=None):
    started = time.perf_counter()
    with (output / "server.stderr").open("wb") as stderr:
        command = ([str(binary), "serve", str(cache), "--host", "127.0.0.1", "--port", "0"]
                   if supported_cli else [str(binary), "--mode=serve", "--host=127.0.0.1",
                                          "--port=0", f"--file={cache}"])
        if cpus is not None:
            command = ["taskset", "-c", ",".join(map(str, cpus))] + command
        child = subprocess.Popen(command + list(extra_args),
                                 stdout=subprocess.PIPE, stderr=stderr, text=True)
    try:
        with selectors.DefaultSelector() as selector:
            selector.register(child.stdout, selectors.EVENT_READ)
            if not selector.select(60):
                raise TimeoutError("server did not become ready within 60 seconds")
            ready = child.stdout.readline()
        (output / "server.stdout").write_text(ready)
        pattern = r"TRACEON_SERVER READY host=127\.0\.0\.1 port=(\d+) records=(\d+)"
        pattern += r" snapshot=([0-9a-f]{64}) protocol=2\n" if supported_cli else r"\n"
        match = re.fullmatch(pattern, ready)
        if not match or child.poll() is not None:
            raise RuntimeError("server failed; see server.stderr")
        return child, int(match[1]), {"startup_seconds": time.perf_counter() - started,
                                    "records": int(match[2]), "ready_pss_kib": pss_kib(child.pid),
                                    "snapshot": match[3] if supported_cli else None,
                                    "command": child.args, "pid": child.pid,
                                    "cpu_affinity": sorted(os.sched_getaffinity(child.pid)) if hasattr(os, "sched_getaffinity") else []}
    except BaseException:
        child.kill()
        child.wait()
        child.stdout.close()
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("minimap2", "seqtk", "cache-builder", "extractor", "server", "index", "reads", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--processes", type=positive, nargs="+", default=[1, 2, 4])
    parser.add_argument("--count", type=positive, default=1000)
    parser.add_argument("--repetitions", type=positive, default=5)
    args = parser.parse_args()
    if platform.system() != "Linux":
        parser.error("Linux /proc and wait4 are required")
    paths = {key: value.resolve(strict=True) for key, value in vars(args).items()
             if isinstance(value, Path) and key != "output"}
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    records = dict(fasta_records(paths["reads"]))
    names = [name for name, _ in fasta_records(paths["reads"])]
    if len(records) != len(names) or args.count > len(names):
        raise ValueError("source must contain unique IDs and enough records")
    manifest = {"platform": platform.platform(), "cpu_affinity": sorted(os.sched_getaffinity(0)),
                "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                "inputs": {k: {"path": str(v), "sha256": sha256(v)} for k, v in paths.items()},
                "cache_policy": "uncontrolled warm OS page cache; one server for the entire study",
                "memory": "sampled simultaneous PSS; server included in shared arm only; controller excluded",
                "source_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                "working_tree": subprocess.check_output(["git", "status", "--short"], text=True)}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    cache = output / "reads.traceon"
    build = run_measured([str(paths["cache_builder"]), str(paths["reads"]), str(cache)], output, "cache-build")
    server, port, setup = start_server(paths["server"], cache, output)
    samples = []
    try:
        if setup["records"] != len(names):
            raise ValueError("server record count differs from source")
        arms = ["seqtk", "local-restore", "shared-server"]
        for processes in args.processes:
            for repetition in range(args.repetitions + 1):
                jobs = []
                for client in range(processes):
                    seed = 42 + processes * 10000 + repetition * 100 + client
                    chosen = set(random.Random(seed).sample(names, args.count))
                    ids = [name for name in names if name in chosen]
                    path = output / f"ids-{processes}-{repetition}-{client}.txt"
                    path.write_text("".join(name + "\n" for name in ids))
                    jobs.append((path, [(name, records[name]) for name in ids]))
                expected_sam = None
                for arm in arms[repetition % 3:] + arms[:repetition % 3]:
                    label = f"{arm}-p{processes}-r{repetition}"
                    commands = []
                    for ids, _ in jobs:
                        if arm == "seqtk":
                            command = [str(paths["seqtk"]), "subseq", str(paths["reads"]), str(ids)]
                        elif arm == "local-restore":
                            command = [str(paths["extractor"]), str(cache), str(ids)]
                        else:
                            command = [str(paths["extractor"]), "--remote", "127.0.0.1", str(port), str(ids)]
                        commands.append(command)
                    resident = server.pid if arm == "shared-server" else None
                    extraction = run_group(commands, output, label + "-extract", resident)
                    for i, (_, expected) in enumerate(jobs):
                        if list(fasta_records(output / f"{label}-extract-{i}.stdout")) != expected:
                            raise ValueError(f"{label}: sequence mismatch for client {i}")
                    commands = [[str(paths["minimap2"]), "-a", "-x", "map-pb", "-t", "1",
                                 str(paths["index"]), str(output / f"{label}-extract-{i}.stdout")]
                                for i in range(processes)]
                    mapping = run_group(commands, output, label + "-map", resident)
                    sam = [sam_summary(output / f"{label}-map-{i}.stdout") for i in range(processes)]
                    if any(s["primary_records"] != args.count for s in sam):
                        raise ValueError(f"{label}: primary record count mismatch")
                    if expected_sam is not None and sam != expected_sam:
                        raise ValueError(f"{label}: alignment mismatch between arms")
                    expected_sam = sam
                    if server.poll() is not None:
                        raise RuntimeError("resident server exited during study")
                    if repetition:
                        samples.append({"arm": arm, "processes": processes, "repetition": repetition,
                                        "extraction": extraction, "mapping": mapping, "sam": sam,
                                        "workflow_seconds": extraction["wall_seconds"] + mapping["wall_seconds"]})
    finally:
        server.terminate()
        try:
            server.wait(timeout=10)
        except subprocess.TimeoutExpired:
            server.kill()
            server.wait()
        server.stdout.close()
    if server.returncode != 0:
        raise RuntimeError(f"server did not shut down cleanly: {server.returncode}")
    aggregates = []
    for processes in args.processes:
        for arm in arms:
            selected = [s for s in samples if s["processes"] == processes and s["arm"] == arm]
            aggregates.append({"processes": processes, "arm": arm,
                               "median_retrieval_seconds": statistics.median(s["extraction"]["wall_seconds"] for s in selected),
                               "median_workflow_seconds": statistics.median(s["workflow_seconds"] for s in selected),
                               "max_retrieval_pss_kib": max(s["extraction"]["sampled_peak_aggregate_pss_kib"] for s in selected),
                               "max_mapping_pss_kib": max(s["mapping"]["sampled_peak_aggregate_pss_kib"] for s in selected)})
    report = {"status": "validated shared-server pilot", "cache_build": build, "server_setup": setup,
              "server_exit": server.returncode, "cache_sha256": sha256(cache),
              "samples": samples, "aggregates": aggregates}
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(aggregates, indent=2))


if __name__ == "__main__":
    main()
