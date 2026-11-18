#!/usr/bin/env python3
"""Linux minimap2 baseline: index construction, repeated mapping, SAM checks.

No TracEon speedup is measured here. Inputs and binaries are supplied explicitly.
Each output directory must be new; raw outputs are retained for inspection.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sam_summary(path):
    """Ignore program provenance, record order and optional-tag order only.

    Preserve all alignment fields, reference headers, qualities, and tag values.
    Sorting retains duplicate records. Memory use is outside the measured child.
    """
    records, headers = [], []
    primary = mapped = 0
    with Path(path).open() as stream:
        for line in stream:
            line = line.rstrip("\n")
            if line.startswith("@"):
                if not line.startswith("@PG\t"):
                    headers.append(line)
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                raise ValueError(f"Malformed SAM record in {path}")
            flag = int(fields[1])
            if not flag & (0x100 | 0x800):
                primary += 1
                mapped += not bool(flag & 0x4)
            records.append("\t".join(fields[:11] + sorted(fields[11:])))
    if not records or not primary:
        raise ValueError(f"No primary SAM records in {path}")
    canonical = "\n".join(sorted(headers) + sorted(records)) + "\n"
    return {"canonical_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
            "records": len(records), "primary_records": primary,
            "mapped_primary_records": mapped}


def run_measured(command, output, label):
    """wait4 reports this child's RSS, not a cumulative maximum over prior runs."""
    started = time.perf_counter()
    with (output / f"{label}.stdout").open("wb") as stdout, \
            (output / f"{label}.stderr").open("wb") as stderr:
        with subprocess.Popen(command, stdout=stdout, stderr=stderr, cwd=output) as child:
            try:
                _, status, usage = os.wait4(child.pid, 0)
                child.returncode = os.waitstatus_to_exitcode(status)
            except BaseException:
                child.kill()
                child.wait()
                raise
    result = {"label": label, "command": command, "returncode": child.returncode,
              "wall_seconds": time.perf_counter() - started,
              "user_seconds": usage.ru_utime, "system_seconds": usage.ru_stime,
              "peak_rss_kib": usage.ru_maxrss}
    with (output / "measurements.jsonl").open("a") as stream:
        stream.write(json.dumps(result) + "\n")
    if child.returncode:
        raise subprocess.CalledProcessError(child.returncode, command)
    return result


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimap2", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preset", default="map-pb")
    parser.add_argument("--threads", type=positive, default=1)
    parser.add_argument("--repetitions", type=positive, default=5)
    parser.add_argument("--expected-reads", type=positive, required=True,
                        help="independent input count; single-end only")
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--dataset-source", required=True)
    args = parser.parse_args()
    if platform.system() != "Linux":
        parser.error("This pilot uses Linux wait4 RSS units")
    binary, reference, reads = [p.resolve(strict=True) for p in
                                (args.minimap2, args.reference, args.reads)]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    version = subprocess.check_output([str(binary), "--version"], text=True).strip()
    manifest = {"minimap2_version": version, "source_revision": args.source_revision,
                "dataset_source": args.dataset_source, "platform": platform.platform(),
                "cpu_affinity": sorted(os.sched_getaffinity(0)),
                "cpu_info": Path("/proc/cpuinfo").read_text(),
                "memory_info": Path("/proc/meminfo").read_text(),
                "cache_policy": "OS page cache uncontrolled; one warmup; fresh processes",
                "arguments": {k: str(v) if isinstance(v, Path) else v
                              for k, v in vars(args).items()},
                "inputs": {str(p): {"sha256": sha256(p), "bytes": p.stat().st_size}
                           for p in (binary, reference, reads)}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    index = output / "reference.mmi"
    build = run_measured([str(binary), "-x", args.preset, "-t", str(args.threads),
                          "-d", str(index), str(reference)], output, "index")
    command = [str(binary), "-a", "-x", args.preset, "-t", str(args.threads),
               str(index), str(reads)]
    samples, expected = [], None
    for i in range(args.repetitions + 1):
        label = "warmup" if i == 0 else f"mapping-{i}"
        result = run_measured(command, output, label)
        result["sam"] = sam_summary(output / f"{label}.stdout")
        if result["sam"]["primary_records"] != args.expected_reads:
            raise ValueError(f"{label}: primary record count differs from input count")
        if expected is not None and result["sam"] != expected:
            raise ValueError(f"{label}: SAM output differs from warmup")
        expected = result["sam"]
        if i:
            samples.append(result)
    report = {"status": "baseline-only; no TracEon comparison", "index": build,
              "index_sha256": sha256(index), "samples": samples,
              "median_mapping_seconds": statistics.median(x["wall_seconds"] for x in samples),
              "max_mapping_rss_kib": max(x["peak_rss_kib"] for x in samples),
              "sam": expected}
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "samples"}, indent=2))


if __name__ == "__main__":
    main()
