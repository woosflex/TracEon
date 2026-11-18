#!/usr/bin/env python3
"""Compare seqtk extraction, TracEon restore/extraction, and a saved read subset.

FASTA-only, single-end minimap2 map-pb pilot. All three arms materialize the same
ordered subset; no claim of a persistent in-process cache or mapper modification.
"""
import argparse
import json
from pathlib import Path
import platform
import random
import statistics
import subprocess

from mapping_baseline import positive, run_measured, sam_summary, sha256


def fasta_records(path):
    name, parts = None, []
    with Path(path).open() as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                header = line[1:].split()
                if not header:
                    raise ValueError("empty FASTA identifier")
                name, parts = header[0], []
            elif name is None:
                raise ValueError("sequence before FASTA header")
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("minimap2", "seqtk", "cache-builder", "extractor", "index", "reads", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--count", type=positive, default=1000)
    parser.add_argument("--repetitions", type=positive, default=5)
    args = parser.parse_args()
    if platform.system() != "Linux":
        parser.error("Linux is required for per-child RSS measurement")
    paths = {key: value.resolve(strict=True) for key, value in vars(args).items()
             if isinstance(value, Path) and key != "output"}
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    manifest = {"status": "feasibility pilot; mutable TracEon working tree",
                "platform": platform.platform(), "subset_seed": 42,
                "cache_policy": "warm/uncontrolled OS page cache; fresh processes",
                "inputs": {k: {"path": str(v), "sha256": sha256(v)} for k, v in paths.items()},
                "source_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                "working_tree": subprocess.check_output(["git", "status", "--short"], text=True),
                "repetitions": args.repetitions, "count": args.count}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    names = [name for name, _ in fasta_records(paths["reads"])]
    if len(names) != len(set(names)):
        raise ValueError("duplicate source IDs are not supported by this pilot")
    if args.count > len(names):
        raise ValueError("requested subset exceeds input record count")
    selected = set(random.Random(42).sample(names, args.count))
    expected = [(name, sequence) for name, sequence in fasta_records(paths["reads"])
                if name in selected]
    ids, saved = output / "ids.txt", output / "expected.fa"
    ids.write_text("".join(name + "\n" for name, _ in expected))
    saved.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in expected))
    cache = output / "reads.traceon"
    build = run_measured([str(paths["cache_builder"]), str(paths["reads"]), str(cache)],
                         output, "cache-build")
    saved_build = run_measured([str(paths["seqtk"]), "subseq", str(paths["reads"]), str(ids)],
                               output, "saved-subset-build")
    if list(fasta_records(output / "saved-subset-build.stdout")) != expected:
        raise ValueError("seqtk saved subset differs from independent FASTA oracle")
    saved = output / "saved-subset-build.stdout"
    samples, expected_sam = [], None
    arms = ["seqtk", "traceon", "saved-subset"]
    for repetition in range(args.repetitions + 1):
        order = arms[repetition % 3:] + arms[:repetition % 3]
        for arm in order:
            label = f"{arm}-{repetition}"
            extraction = None
            subset = saved
            if arm != "saved-subset":
                command = ([str(paths["seqtk"]), "subseq", str(paths["reads"]), str(ids)]
                           if arm == "seqtk" else [str(paths["extractor"]), str(cache), str(ids)])
                extraction = run_measured(command, output, label + "-extract")
                subset = output / f"{label}-extract.stdout"
                if list(fasta_records(subset)) != expected:
                    raise ValueError(f"{label}: extracted identifiers/order/sequences differ")
            mapping = run_measured([str(paths["minimap2"]), "-a", "-x", "map-pb", "-t", "1",
                                    str(paths["index"]), str(subset)], output, label + "-map")
            sam = sam_summary(output / f"{label}-map.stdout")
            if sam["primary_records"] != args.count:
                raise ValueError(f"{label}: primary count differs from requested subset")
            if expected_sam is not None and sam != expected_sam:
                raise ValueError(f"{label}: alignment differs across arms/repetitions")
            expected_sam = sam
            if repetition:
                samples.append({"arm": arm, "repetition": repetition, "extraction": extraction,
                                "mapping": mapping,
                                "total_seconds": mapping["wall_seconds"] +
                                    (extraction["wall_seconds"] if extraction else 0),
                                "peak_rss_kib": max(mapping["peak_rss_kib"],
                                    extraction["peak_rss_kib"] if extraction else 0)})
    medians = {arm: statistics.median(s["total_seconds"] for s in samples if s["arm"] == arm)
               for arm in arms}
    report = {"status": "output-validated FASTA subset-remapping pilot", "samples": samples,
              "median_workflow_seconds": medians, "cache_build": build,
              "saved_subset_build": saved_build, "sam": expected_sam,
              "ids_sha256": sha256(ids), "cache_sha256": sha256(cache)}
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "samples"}, indent=2))


if __name__ == "__main__":
    main()
