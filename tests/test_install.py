"""Install, relocate, compile and run a consumer using only the installed package."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

cmake, build, source, compiler, config = sys.argv[1:]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    original, relocated = root / "original", root / "relocated"
    subprocess.run([cmake, "--install", build, "--prefix", str(original),
                    "--component", "TracEon", "--config", config], check=True)
    original.rename(relocated)
    for manifest in relocated.rglob("*.cmake"):
        text = manifest.read_text()
        assert str(Path(build).resolve()) not in text, manifest
        assert str(Path(source).resolve().parents[1]) not in text, manifest
    consumer = root / "consumer"
    shutil.copytree(source, consumer)
    subprocess.run([cmake, "-S", str(consumer), "-B", str(root / "build"),
                    f"-DCMAKE_PREFIX_PATH={relocated}", f"-DCMAKE_CXX_COMPILER={compiler}",
                    "-DCMAKE_BUILD_TYPE=Release"], check=True)
    subprocess.run([cmake, "--build", str(root / "build"), "--config", "Release", "-j2"], check=True)
    executable = root / "build" / "check_consumer"
    if not executable.exists():
        executable = root / "build" / "Release" / "check_consumer.exe"
    subprocess.run([str(executable)], cwd=root, check=True, timeout=10)
    cli = relocated / "bin" / "traceon"
    if cli.exists():
        subprocess.run([str(cli), "--help"], check=True, capture_output=True, timeout=10)
print("relocated installed consumer passed")
