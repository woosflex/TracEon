"""Run with cache-builder, subset-extractor, and remote_bench executable paths."""
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from shared_mapping_pilot import run_group, start_server

builder, extractor, server_binary = [Path(p).resolve() for p in sys.argv[1:]]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    fasta, cache, ids = root / "reads.fa", root / "reads.traceon", root / "ids.txt"
    fasta.write_text(">a description\nacgTN\n>b\nTTAA\n")
    ids.write_text("b\na\n")
    subprocess.run([str(builder), str(fasta), str(cache)], check=True, capture_output=True)
    server, port, setup = start_server(server_binary, cache, root)
    try:
        assert setup["records"] == 2
        # Once loaded, clients do not need access to either source or cache file.
        cache.unlink()
        fasta.unlink()
        command = [str(extractor), "--remote", "127.0.0.1", str(port), str(ids)]
        result = run_group([command, command], root, "clients", server.pid)
        assert result["includes_server"] and len(result["children"]) == 2
        assert result["sampled_peak_aggregate_pss_kib"] >= setup["ready_pss_kib"]
        for client in range(2):
            assert (root / f"clients-{client}.stdout").read_text() == ">b\nTTAA\n>a\nacgTN\n"
        ids.write_text("a\nmissing\n")
        missing = subprocess.run(command, capture_output=True)
        assert missing.returncode != 0 and b"missing ID" in missing.stderr and not missing.stdout
        try:
            run_group([[sys.executable, "-c", "raise SystemExit(7)"]], root, "failure")
            raise AssertionError("failed process group accepted")
        except subprocess.CalledProcessError as error:
            assert error.returncode == 7
    finally:
        server.terminate()
        server.wait(timeout=10)
        server.stdout.close()
    assert server.returncode == 0
print("shared-server pilot checks passed")
