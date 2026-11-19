#!/usr/bin/env python3
"""Small integration checks for audit regressions; uses only the Python stdlib."""
import argparse
import ast
import gzip
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import select
import socket
import time
import types
import typing

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--thread-shim", type=Path)
    args = parser.parse_args()
    build = args.build_dir.resolve()
    if args.thread_shim: args.thread_shim = args.thread_shim.resolve()
    suffix = ".exe" if os.name == "nt" else ""
    if not (build / ("traceon_driver" + suffix)).exists() and (build / "Release" / ("traceon_driver" + suffix)).exists():
        build = build / "Release"
    driver = build / ("traceon_driver" + suffix)
    scaling = build / ("read_scaling" + suffix)
    remote = build / ("remote_bench" + suffix)

    def run(command, **kwargs):
        return subprocess.run([str(x) for x in command], capture_output=True,
                              text=True, timeout=30, **kwargs)

    spec = importlib.util.spec_from_file_location("value_chain", ROOT / "benchmarks/value_chain_benchmark.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    try:
        module.run_timed([sys.executable, "-c", "raise SystemExit(7)"])
    except subprocess.CalledProcessError as error:
        require(error.returncode == 7, "failed command status changed")
    else:
        raise AssertionError("benchmark swallowed a subprocess failure")
    with tempfile.TemporaryDirectory(prefix="traceon-audit-checks-") as directory:
        tmp = Path(directory)
        (tmp / "pyfastx.py").write_text(
            "class Fasta:\n"
            " def __init__(self, path): pass\n"
            " def keys(self): return ['r0']\n"
            " def __getitem__(self, key): raise KeyError(key)\n"
            "Fastq = Fasta\n")
        env = dict(os.environ, PYTHONPATH=str(tmp))
        for script in (module.PYFASTX_LOOKUP_PREFIX_SCRIPT, module.PYFASTX_LOOKUP_REALIDS_SCRIPT):
            result = run([sys.executable, "-c", script, "reads.fastq", "1", "r", "1"], env=env)
            require(result.returncode != 0 and "KeyError" in result.stderr,
                    "lookup script swallowed a missing record")
        previous_cwd = Path.cwd()
        os.chdir(tmp)
        try:
            sentinel = tmp / "existing.fai"
            sentinel.write_text("keep me")
            old_scratch = tmp / ".valuechain_tmp"
            old_scratch.mkdir()
            (old_scratch / "keep").write_text("keep me too")
            module.TRACEON_BINARY = str(driver)
            workspaces = []

            def section_a():
                workspaces.append(Path(module.TMP_DIR))
                (workspaces[-1] / "owned.fxi").write_text("temporary")

            module.section_a = section_a
            module.section_b = lambda: None
            module.main()
            require(sentinel.read_text() == "keep me", "pre-existing index deleted")
            require((old_scratch / "keep").exists(), "pre-existing scratch directory deleted")
            require(not workspaces[0].exists(), "owned workspace leaked")

            # Load the local validation routine without its optional download dependencies.
            tree = ast.parse((ROOT / "benchmarks/validate_real_data.py").read_text())
            routine = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                           and node.name == "run_traceon_test")
            namespace = dict(Path=Path, tempfile=tempfile, time=time, subprocess=subprocess,
                             Optional=typing.Optional, Dict=typing.Dict, Any=typing.Any,
                             DEFAULT_LOOKUPS=1, Colors=types.SimpleNamespace(CYAN="", ENDC=""),
                             get_file_size_mb=lambda path: path.stat().st_size / (1024 * 1024))
            for name in ("print_success", "print_info", "print_error"):
                namespace[name] = lambda *args: None
            exec(compile(ast.Module(body=[routine], type_ignores=[]), "validation", "exec"), namespace)
            validation_source = tmp / "validation.fasta"
            validation_source.write_text(">r\nACGT\n")
            validation_cache = validation_source.with_suffix(".traceon")
            validation_cache.write_text("pre-existing cache")
            result = namespace["run_traceon_test"](driver, validation_source, lookups=4)
            require(result is not None, "validation failed")
            require(validation_cache.read_text() == "pre-existing cache", "validation replaced a user cache")

            base_args = [scaling, "--records=1", "--keys=1", "--count=1", "--warmup=0", "--threads=1"]
            for bad in ("--reps=0", "--threads=1,0", "--threads=1,-1", "--threads=1,",
                        "--count=garbage", "--records=-1", "--records=0"):
                result = run(base_args + [bad])
                require(result.returncode > 0, f"bad benchmark option crashed or succeeded: {bad}")
            result = run(base_args + ["--reps=1"])
            require(result.returncode == 0, result.stderr)
            require("hit_rate=1.0000" in result.stdout, "valid benchmark lost reads")
            if remote.exists():
                for bad in ("--keys=0", "--threads=0", "--count=-1", "--port=65536"):
                    result = run([remote, "--mode=remote", "--file=unused.traceon", bad])
                    require(result.returncode > 0, f"invalid remote option accepted/crashed: {bad}")

            if remote.exists():
                source = tmp / "remote.fasta"
                source.write_text(">a\nACGT\n>b\nTGCA\n")
                binary = tmp / "remote.traceon"
                require(run([driver, "save", source, binary]).returncode == 0, "cache build failed")
                server = subprocess.Popen([str(remote), "--mode=serve", "--port=0", "--file=" + str(binary)],
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                try:
                    require(select.select([server.stdout], [], [], 5)[0], "server startup timed out")
                    ready = server.stdout.readline()
                    require("READY" in ready, "server failed to start: " + ready)
                    port = int(next(field.split("=", 1)[1] for field in ready.split() if field.startswith("port=")))
                    result = run([remote, "--mode=remote", "--port=" + str(port), "--file=" + str(binary),
                                  "--keys=2", "--count=17", "--threads=4", "--warmup=2"])
                    require(result.returncode == 0, result.stderr)
                    fields = dict(item.split("=", 1) for item in result.stdout.split() if "=" in item)
                    require(fields["count"] == "17", "rounded-up work count")
                    expected_ops = 17 / float(fields["wall_s"])
                    require(abs(float(fields["ops"]) - expected_ops) < 2,
                            "throughput does not use elapsed wall time")

                    # Sequential healthcheck-style connections must not retain one stack each.
                    def churn(n):
                        for _ in range(n):
                            with socket.create_connection(("127.0.0.1", port), timeout=5) as client:
                                client.shutdown(socket.SHUT_WR)
                                require(client.recv(1) == b"", "server did not close idle connection")

                    churn(16)
                    status = Path(f"/proc/{server.pid}/status")
                    def virtual_kib():
                        return int(next(line.split()[1] for line in status.read_text().splitlines()
                                        if line.startswith("VmSize:")))
                    before = virtual_kib() if args.thread_shim else None
                    churn(128)
                    if before is not None:
                        require(virtual_kib() - before < 256 * 1024, "completed worker stacks accumulate")
                finally:
                    server.terminate()
                    try: server.communicate(timeout=10)
                    except subprocess.TimeoutExpired:
                        server.kill()
                        server.communicate()
                        raise AssertionError("server stop hung")

            # Execute the actual CI seed generator functions, not copies of them.
            workflow = (ROOT / ".github/workflows/fuzz.yml").read_text()
            block = workflow.split("python3 - <<'EOF'", 1)[1].split("        EOF", 1)[0]
            import textwrap
            tree = ast.parse(textwrap.dedent(block))
            functions = ast.Module(body=[node for node in tree.body if isinstance(node, ast.FunctionDef)], type_ignores=[])
            namespace = {}
            exec(compile(functions, "ci-seed-generators", "exec"), namespace)
            seed = tmp / "seed.fastq"
            namespace["gen_fastq"](seed, 36000)
            result = run([driver, "load", seed])
            require(result.returncode == 0 and "Sequences: 36000" in result.stdout, result.stderr)

            if args.thread_shim and (os.cpu_count() or 1) > 1:
                env = dict(os.environ, LD_PRELOAD=str(args.thread_shim.resolve()))
                fasta = tmp / "large.fasta"
                fasta.write_text(">r\n" + "A" * (11 * 1024 * 1024) + "\n")
                # Stored GZIP members exceed the parallel-decompression threshold.
                concatenated = tmp / "parallel.fasta.gz"
                concatenated.write_bytes(gzip.compress(fasta.read_bytes(), compresslevel=0)
                                         + gzip.compress(b">s\nACGT\n", compresslevel=0))
                for source in (seed, fasta, concatenated):
                    result = run([driver, "load", source], env=env)
                    require(result.returncode == 1 and "Exception:" in result.stderr,
                            f"thread launch failure did not propagate for {source}: {result.returncode} {result.stderr}")
        finally:
            os.chdir(previous_cwd)
    print("audit integration checks passed")


if __name__ == "__main__":
    main()
