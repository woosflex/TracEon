#!/usr/bin/env python3
"""Sustained Linux service stress with raw resource samples and recovery checks."""
import argparse
from contextlib import ExitStack
import json
import os
from pathlib import Path
import platform
import socket
import struct
import subprocess
import threading
import time

from mapping_baseline import positive, run_measured, sha256
from shared_mapping_pilot import pss_kib, start_server


def frame(kind, body=b""):
    return struct.pack("<Q", len(body) + 1) + bytes([kind]) + body


def receive(sock):
    def exact(n):
        parts = bytearray()
        while len(parts) < n:
            data = sock.recv(n - len(parts))
            if not data:
                raise EOFError("peer closed")
            parts.extend(data)
        return bytes(parts)
    size, = struct.unpack("<Q", exact(8))
    if not 1 <= size <= 64 * 1024 * 1024:
        raise ValueError("invalid response length")
    return exact(size)


def resource_sample(pid):
    status = dict(line.split(":", 1) for line in Path(f"/proc/{pid}/status").read_text().splitlines())
    return {"server_pss_kib": pss_kib(pid), "server_rss_kib": int(status["VmRSS"].split()[0]),
            "threads": int(status["Threads"]), "fds": len(list(Path(f"/proc/{pid}/fd").iterdir())),
            "aggregate_pss_kib": pss_kib(pid) + pss_kib(os.getpid())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--traceon", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--connections", type=positive, default=4)
    parser.add_argument("--cycles", type=positive, default=3)
    parser.add_argument("--phase-seconds", type=positive, default=2)
    args = parser.parse_args()
    harness = {p.name: sha256(p) for p in (Path(__file__), Path(__file__).with_name("shared_mapping_pilot.py"),
                Path(__file__).with_name("mapping_baseline.py"), Path(__file__).with_name("subset_mapping_pilot.py"))}
    binary = args.traceon.resolve(strict=True)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    source, cache = output / "stress.fa", output / "stress.traceon"
    source.write_text(">small\nACGT\n>large\n" + "ACGT" * (2 * 1024 * 1024) + "\n")
    build = run_measured([str(binary), "build", str(source), str(cache)], output, "build")
    server, port, setup = start_server(binary, cache, output, supported_cli=True,
        extra_args=["--max-connections", str(args.connections), "--timeout-ms", "250"])
    stop = threading.Event()
    phase = ["startup"]
    samples, errors, phases, settled = [], [], [], []
    started = time.monotonic()

    def sample():
        try:
            with (output / "resources.jsonl").open("w") as stream:
                while not stop.is_set():
                    row = {"seconds": time.monotonic() - started, "phase": phase[0], **resource_sample(server.pid)}
                    samples.append(row)
                    stream.write(json.dumps(row) + "\n")
                    stop.wait(0.01)
        except BaseException as error:
            errors.append(repr(error))

    def connect(handshake=False):
        sock = socket.create_connection(("127.0.0.1", port), timeout=2)
        try:
            if handshake:
                sock.sendall(frame(1, b"TRO-RPC\x01"))
                response = receive(sock)
                if response[:6] != b"\x81\x00\x00\x00\x00\x00":
                    raise ValueError("bad HELLO response")
            return sock
        except BaseException:
            sock.close()
            raise

    def get_small(sock):
        sock.sendall(frame(2, b"small"))
        response = receive(sock)
        if response[:10] != b"\x81\x00\x04\x00\x00\x00ACGT":
            raise ValueError("incorrect admitted query")

    def probe():
        command = [str(binary), "query", "--port", str(port), "--timeout-ms", "1000",
                   "--expect-snapshot", setup["snapshot"], "small"]
        deadline = time.monotonic() + 3
        while True:
            result = subprocess.run(command, capture_output=True, timeout=3)
            if result.returncode == 0:
                if result.stdout != b">small\nACGT\n":
                    raise ValueError("recovery query mismatch")
                return
            if time.monotonic() >= deadline:
                raise RuntimeError(f"capacity did not recover: {result.stderr!r}")
            time.sleep(0.01)

    monitor = threading.Thread(target=sample)
    monitor.start()
    try:
        baseline = resource_sample(server.pid)
        for cycle in range(args.cycles):
            for name in ("churn", "overload", "idle", "slow-reader", "partial-frame", "disconnect"):
                phase[0] = f"c{cycle}-{name}"
                deadline = time.monotonic() + args.phase_seconds
                attempts = rejected = queries = 0
                while time.monotonic() < deadline:
                    if name in ("churn", "disconnect"):
                        try:
                            with connect(name == "churn") as sock:
                                attempts += 1
                                if name == "churn":
                                    get_small(sock)
                                    queries += 1
                                elif attempts % 3 == 0:
                                    # Oversized length is rejected before allocation.
                                    sock.sendall(struct.pack("<Q", 64 * 1024 * 1024 + 1))
                                    if receive(sock)[0] != 0x82:
                                        raise ValueError("oversized frame accepted")
                                else:
                                    # EOF during HELLO/header/payload exercises cleanup.
                                    sock.sendall(frame(1, b"TRO-RPC\x02")[:3 if attempts % 2 else 12])
                        except (EOFError, ConnectionResetError, BrokenPipeError):
                            rejected += 1  # Admission may reject rapid fault/churn traffic.
                        continue
                    with ExitStack() as stack:
                        sockets = [stack.enter_context(connect(name in ("overload", "slow-reader")))
                                   for _ in range(args.connections)]
                        attempts += len(sockets)
                        if name == "overload":
                            for _ in range(16):
                                for sock in sockets:
                                    get_small(sock)
                                    queries += 1
                                with connect() as excess:
                                    if excess.recv(1) != b"":
                                        raise ValueError("excess connection was not rejected")
                                rejected += 1
                                attempts += 1
                        else:
                            for sock in sockets:
                                if name == "slow-reader":
                                    sock.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 1024)
                                    sock.sendall(frame(2, b"large"))
                                elif name == "partial-frame":
                                    sock.sendall(struct.pack("<Q", 64 * 1024 * 1024) + b"\x01")
                            # Keep sockets open beyond the 250ms server deadline.
                            time.sleep(0.4)
                            probe()  # Capacity must recover even before these sockets close.
                    # Reap completed workers before starting the next fully occupied cohort.
                    time.sleep(0.002)
                phase[0] = f"c{cycle}-{name}-recovery"
                probe()
                time.sleep(0.4)
                current = resource_sample(server.pid)
                if current["threads"] != baseline["threads"] or current["fds"] != baseline["fds"]:
                    raise RuntimeError(f"resources did not recover after {name}: {current}")
                phases.append({"cycle": cycle, "phase": name, "attempts": attempts,
                               "rejected": rejected, "successful_queries": queries, "settled": current})
            settled.append(resource_sample(server.pid))
        if errors:
            raise RuntimeError(errors)
        peak_threads = max(s["threads"] for s in samples)
        if peak_threads > baseline["threads"] + args.connections:
            raise RuntimeError("worker count exceeded admission cap")
        # After a complete warmup cycle, resource retention must plateau.
        late_growth = max(s["server_pss_kib"] for s in settled[1:]) - settled[0]["server_pss_kib"] if len(settled) > 1 else 0
        if late_growth > 8192:
            raise RuntimeError(f"post-warmup retained PSS grew by {late_growth} KiB")
        report = {"status": "passed", "binary_sha256": sha256(binary), "build": build,
                  "harness_sha256": harness, "fixture_sha256": sha256(source),
                  "platform": platform.platform(), "cpu_affinity": sorted(os.sched_getaffinity(0)),
                  "setup": setup, "arguments": vars(args) | {"traceon": str(binary), "output": str(output)},
                  "baseline": baseline, "phases": phases, "settled_cycles": settled,
                  "samples": len(samples), "peak_threads": peak_threads,
                  "peak_server_pss_kib": max(s["server_pss_kib"] for s in samples),
                  "peak_aggregate_pss_kib": max(s["aggregate_pss_kib"] for s in samples),
                  "post_warmup_pss_growth_kib": late_growth,
                  "memory_scope": "10ms sampled PSS; aggregate includes server and Python load-generator/controller",
                  "source_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                  "working_tree": subprocess.check_output(["git", "status", "--short"], text=True)}
    finally:
        stop.set()
        monitor.join(timeout=3)
        server.terminate()
        try:
            server.wait(timeout=5)
        except subprocess.TimeoutExpired:
            server.kill()
            server.wait()
        server.stdout.close()
    if server.returncode:
        raise RuntimeError(f"unclean server shutdown: {server.returncode}")
    if errors:
        raise RuntimeError(errors)
    report["server_exit"] = server.returncode
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "samples", "peak_threads", "peak_server_pss_kib", "post_warmup_pss_growth_kib")}))


if __name__ == "__main__":
    main()
