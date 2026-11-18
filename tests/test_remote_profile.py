"""Benchmark correctness: no client cache, pinned FASTQ, raw timing samples."""
import math
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from shared_mapping_pilot import start_server

traceon, benchmark = [Path(p).resolve() for p in sys.argv[1:]]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    source, cache, ids, raw = (root / name for name in ("reads.fq", "reads.traceon", "ids.txt", "latencies.txt"))
    source.write_text("@a\nACGT\n+\nI@+?\n@empty\n\n+\n\n")
    subprocess.run([str(traceon), "build", str(source), str(cache)], check=True, capture_output=True)
    server, port, setup = start_server(traceon, cache, root, supported_cli=True)
    try:
        source.unlink()
        cache.unlink()
        ids.write_text("a\n")
        command = [str(benchmark), "--mode", "remote", "--port", str(port), "--ids", str(ids),
                   "--expect-snapshot", setup["snapshot"], "--count", "32", "--threads", "2", "--warmup", "3",
                   "--latencies", str(raw)]
        result = subprocess.run(command + ["--fastq"], check=True, capture_output=True, text=True, timeout=10)
        assert "TRACEON_PAYLOAD\tbytes=256\tfastq=1" in result.stdout
        times = [float(x) for x in raw.read_text().splitlines()]
        assert len(times) == 32 and all(math.isfinite(x) and x >= 0 for x in times)
        result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=10)
        assert "TRACEON_PAYLOAD\tbytes=128\tfastq=0" in result.stdout
        result = subprocess.run(command + ["--fastq", "--batch-size", "4"], check=True,
                                capture_output=True, text=True, timeout=10)
        assert "TRACEON_PAYLOAD\tbytes=1024\tfastq=1" in result.stdout
        assert "TRACEON_REQUESTS\tbatch_size=4\tlatency_unit=request" in result.stdout
        assert len(raw.read_text().splitlines()) == 32
        for size in ("0", "257"):
            assert subprocess.run(command + ["--batch-size", size], capture_output=True).returncode != 0
        ids.write_text("empty\n")
        result = subprocess.run(command + ["--fastq"], check=True, capture_output=True, text=True, timeout=10)
        assert "TRACEON_PAYLOAD\tbytes=0\tfastq=1" in result.stdout
        for extras in (["--expect-snapshot", "0" * 64], ["--file", str(cache)]):
            assert subprocess.run(command + extras, capture_output=True, timeout=10).returncode != 0
        ids.write_text("missing\n")
        assert subprocess.run(command, capture_output=True, timeout=10).returncode != 0
        ids.write_text("\n")
        assert subprocess.run(command, capture_output=True, timeout=10).returncode != 0
    finally:
        server.terminate()
        server.wait(timeout=10)
        server.stdout.close()
    assert server.returncode == 0
    for size in (0, 200, 10000):
        subprocess.run([str(benchmark), "--mode", "codec", "--bytes", str(size), "--count", "32",
                        "--warmup", "3", "--latencies", str(raw)], check=True, capture_output=True)
        lines = raw.read_text().splitlines()
        assert lines[0].startswith("encode_crc_us\t") and len(lines) == 33
        assert all(len(row.split()) == 5 for row in lines[1:])
        assert all(math.isfinite(float(x)) and float(x) >= 0 for row in lines[1:] for x in row.split())
    assert subprocess.run([str(benchmark), "--mode", "codec", "--bytes", str(64 * 1024 * 1024)],
                          capture_output=True).returncode != 0
print("remote profile checks passed")
