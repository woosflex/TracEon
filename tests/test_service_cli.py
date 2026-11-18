"""CLI integration and independent SHA-256 oracle for the served logical snapshot."""
from pathlib import Path
import hashlib
import selectors
import struct
import subprocess
import sys
import tempfile

binary = str(Path(sys.argv[1]).resolve())
for extension, format_byte, source_text, records in (
    ("fa", 0, ">a description\nACGTN\n>b\nTTAA\n", {"a": ("ACGTN", ""), "b": ("TTAA", "")}),
    ("fq", 3, "@a description\nACGTN\n+description\n!@+IJ\n@b\n\n+\n\n",
     {"a": ("ACGTN", "!@+IJ"), "b": ("", "")}),
):
    fastq = extension == "fq"
    digest = hashlib.sha256(b"TracEon snapshot v1\n" + bytes([format_byte]) + struct.pack("<Q", len(records)))
    for key, (seq, qual) in sorted(records.items()):
        for field in (key, seq, qual):
            raw = field.encode()
            digest.update(struct.pack("<Q", len(raw)))
            digest.update(raw)
    expected_id = digest.hexdigest()
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source, cache = root / f"reads.{extension}", root / "reads.traceon"
        source.write_text(source_text)
        subprocess.run([binary, "build", str(source), str(cache)] + sys.argv[2:], check=True, capture_output=True)
        server = subprocess.Popen([binary, "serve", str(cache), "--port", "0", "--max-connections", "2"],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            with selectors.DefaultSelector() as selector:
                selector.register(server.stdout, selectors.EVENT_READ)
                assert selector.select(10), "server readiness timed out"
            line = server.stdout.readline()
            assert line.startswith("TRACEON_SERVER READY "), (line, server.poll())
            metadata = dict(x.split("=", 1) for x in line.split()[2:])
            assert metadata["snapshot"] == expected_id, metadata
            port = metadata["port"]
            source.unlink()
            cache.unlink()
            query = [binary, "query", "--port", port]
            clients = [subprocess.Popen(query + ["--expect-snapshot", expected_id, "b", "a"],
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) for _ in range(2)]
            expected_fasta = "".join(f">{key}\n{records[key][0]}\n" for key in ("b", "a"))
            for client in clients:
                stdout, stderr = client.communicate(timeout=10)
                assert client.returncode == 0, stderr
                assert stdout == expected_fasta
            result = subprocess.run([binary, "info", "--port", port], capture_output=True, text=True, timeout=10)
            assert result.returncode == 0, result.stderr
            assert result.stdout == f"protocol=2 snapshot={expected_id} records=2 sequence=1 fastq={int(fastq)} batch=1\n"
            result = subprocess.run(query + ["--expect-snapshot", "0" * 64, "a"], capture_output=True, timeout=10)
            assert result.returncode == 1 and b"snapshot mismatch" in result.stderr and not result.stdout
            result = subprocess.run(query + ["--protocol", "1", "b", "a"], capture_output=True, text=True, timeout=10)
            assert result.returncode == 0 and result.stdout == expected_fasta, result.stderr
            result = subprocess.run(query + ["--ids", "-"], input="a\r\na\n", text=True,
                                    capture_output=True, timeout=10)
            assert result.returncode == 0, result.stderr
            assert result.stdout == ">a\nACGTN\n>a\nACGTN\n"
            for size in (2, 256):
                result = subprocess.run(query + ["--batch-size", str(size), "b", "a", "a"], capture_output=True, text=True, timeout=10)
                assert result.returncode == 0, result.stderr
                assert result.stdout == expected_fasta + ">a\nACGTN\n"
                result = subprocess.run(query + ["--batch-size", str(size), "--ids", "-"],
                                        input="a\na\n", capture_output=True, text=True, timeout=10)
                assert result.returncode == 0 and result.stdout == ">a\nACGTN\n>a\nACGTN\n", result.stderr
                result = subprocess.run(query + ["--batch-size", str(size), "a", "absent", "b"], capture_output=True, text=True, timeout=10)
                assert result.returncode == 3 and result.stdout == ">a\nACGTN\n"
                if fastq:
                    result = subprocess.run(query + ["--batch-size", str(size), "--fastq", "b", "a", "b"], capture_output=True, text=True, timeout=10)
                    assert result.returncode == 0 and result.stdout == "@b\n\n+\n\n@a\nACGTN\n+\n!@+IJ\n@b\n\n+\n\n", result.stderr
            result = subprocess.run(query + ["--fastq", "b", "a"], capture_output=True, text=True, timeout=10)
            if fastq:
                assert result.returncode == 0, result.stderr
                assert result.stdout == "@b\n\n+\n\n@a\nACGTN\n+\n!@+IJ\n"
            else:
                assert result.returncode == 2 and not result.stdout
            for mode in (["--fastq"], []) if fastq else ([],):
                missing = subprocess.run(query + mode + ["absent"], capture_output=True, timeout=10)
                assert missing.returncode == 3 and not missing.stdout
            for args in ([], ["--port", "0", "a"], ["--port", "65536", "a"],
                         ["--ids", "-", "a"], ["--timeout-ms", "-1", "a"],
                         ["--expect-snapshot", "", "a"], ["--expect-snapshot", "bad", "a"],
                         ["--protocol", "1", "--fastq", "a"], ["--protocol", "3", "a"],
                         ["--batch-size", "0", "a"], ["--batch-size", "257", "a"],
                         ["--protocol", "1", "--batch-size", "2", "a"],
                         ["--protocol", "1", "--expect-snapshot", expected_id, "a"]):
                assert subprocess.run(query + args, capture_output=True, timeout=10).returncode == 2
        finally:
            server.terminate()
            stdout, stderr = server.communicate(timeout=10)
        assert server.returncode == 0, stderr
print("service CLI and snapshot identity checks passed")
