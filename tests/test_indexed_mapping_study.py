"""Independent FASTQ payload checks used by the indexed workflow harness."""
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from indexed_mapping_study import sequence_records

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "reads.fq"
    path.write_text("@a description\nAC\n+description\n@+\n@empty\n\n+\n\n")
    assert list(sequence_records(path, True)) == [("a", "AC", "@+"), ("empty", "", "")]
    for malformed in ("@a\nAC\n+\nI\n", "@a\nAC\n", "@\nAC\n+\nII\n",
                      "@a\nAC\nwrong\nII\n"):
        path.write_text(malformed)
        try:
            list(sequence_records(path, True))
        except ValueError:
            pass
        else:
            raise AssertionError("malformed FASTQ accepted")
print("indexed mapping payload checks passed")
