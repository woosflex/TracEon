"""Standalone checks: python tests/test_mapping_baseline.py."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile

spec = importlib.util.spec_from_file_location(
    "mapping_baseline", Path(__file__).resolve().parents[1] / "benchmarks/mapping_baseline.py")
baseline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(baseline)

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    a, b = root / "a.sam", root / "b.sam"
    record = "q\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIIII"
    a.write_text("@SQ\tSN:ref\tLN:4\n@PG\tID:one\n" + record + "\tNM:i:0\tAS:i:4\n")
    b.write_text("@SQ\tSN:ref\tLN:4\n@PG\tID:two\n" + record + "\tAS:i:4\tNM:i:0\n")
    assert baseline.sam_summary(a) == baseline.sam_summary(b)
    second = record.replace("q\t", "r\t", 1) + "\n"
    a.write_text(a.read_text() + second)
    b.write_text("@SQ\tSN:ref\tLN:4\n" + second + record + "\tAS:i:4\tNM:i:0\n")
    assert baseline.sam_summary(a) == baseline.sam_summary(b)
    b.write_text(a.read_text() + second)
    assert baseline.sam_summary(a) != baseline.sam_summary(b)
    b.write_text(a.read_text().replace("4M", "3M1S"))
    assert baseline.sam_summary(a) != baseline.sam_summary(b)
    b.write_text(a.read_text().replace("IIII", "IIHI"))
    assert baseline.sam_summary(a) != baseline.sam_summary(b)
    b.write_text("@SQ\tSN:ref\tLN:4\n")
    try:
        baseline.sam_summary(b)
        raise AssertionError("empty alignment accepted")
    except ValueError:
        pass
    try:
        baseline.run_measured([sys.executable, "-c", "raise SystemExit(7)"], root, "failure")
        raise AssertionError("failed child accepted")
    except subprocess.CalledProcessError as error:
        assert error.returncode == 7
    measured = baseline.run_measured([sys.executable, "-c", "print('ok')"], root, "success")
    assert measured["returncode"] == 0 and measured["peak_rss_kib"] > 0
    assert (root / "success.stdout").read_text() == "ok\n"
print("mapping baseline checks passed")
