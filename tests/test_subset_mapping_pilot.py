"""Run: python tests/test_subset_mapping_pilot.py build/traceon_make_cache build/traceon_subset"""
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from subset_mapping_pilot import fasta_records

builder, extractor = [str(Path(p).resolve()) for p in sys.argv[1:]]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    fasta, cache, ids = root / "reads.fa", root / "reads.traceon", root / "ids.txt"
    fasta.write_text(">a description\nacg\nTN\n>b\nTTAA\n")
    assert list(fasta_records(fasta)) == [("a", "acgTN"), ("b", "TTAA")]
    subprocess.run([builder, str(fasta), str(cache)], check=True, capture_output=True)
    ids.write_text("b\na\n")
    result = subprocess.run([extractor, str(cache), str(ids)], check=True, capture_output=True, text=True)
    assert result.stdout == ">b\nTTAA\n>a\nacgTN\n"
    for invalid in ("a\nmissing\n", "a\na\n", "", "a description\n"):
        ids.write_text(invalid)
        result = subprocess.run([extractor, str(cache), str(ids)], capture_output=True)
        assert result.returncode != 0 and not result.stdout
    fastq = root / "reads.fq"
    fastq.write_text("@a\nACGT\n+\nIIII\n")
    subprocess.run([builder, str(fastq), str(cache)], check=True, capture_output=True)
    ids.write_text("a\n")
    result = subprocess.run([extractor, str(cache), str(ids)], capture_output=True)
    assert result.returncode != 0 and b"FASTA only" in result.stderr and not result.stdout
print("subset adapter checks passed")
