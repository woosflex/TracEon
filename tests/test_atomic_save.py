"""Linux fault injection proves pre-commit failures leave the old cache intact."""
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time

builder, driver, shim = [str(Path(p).resolve()) for p in sys.argv[1:4]]
storage = sys.argv[4:]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    old, new, cache = root / "old.fa", root / "new.fa", root / "cache.traceon"
    old.write_text(">old\nACGT\n")
    new.write_text(">new\nTTTT\n>second\nGGGG\n")
    subprocess.run([builder, str(old), str(cache)] + storage, check=True, capture_output=True)
    original = cache.read_bytes()
    cache.chmod(0o640)
    for fault in ("write", "file-sync", "rename"):
        env = dict(os.environ, LD_PRELOAD=shim, TRACEON_TEST_ATOMIC_FAULT=fault)
        result = subprocess.run([builder, str(new), str(cache)] + storage, env=env, capture_output=True, timeout=10)
        assert result.returncode != 0, fault
        assert cache.read_bytes() == original, fault
        assert not list(root.glob(".traceon-tmp-*")), fault
        subprocess.run([driver, "restore", str(cache)], check=True, capture_output=True)
    env = dict(os.environ, LD_PRELOAD=shim, TRACEON_TEST_ATOMIC_FAULT="stop-before-rename")
    child = subprocess.Popen([builder, str(new), str(cache)] + storage, env=env,
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        stopped = False
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            pid, status = os.waitpid(child.pid, os.WNOHANG | os.WUNTRACED)
            if pid:
                if os.WIFSTOPPED(status):
                    stopped = True
                else:
                    child.returncode = os.waitstatus_to_exitcode(status)
                break
            time.sleep(0.005)
        assert stopped, "writer never reached the synchronized pre-rename boundary"
        assert cache.read_bytes() == original
        child.kill()
        child.wait(timeout=10)
        assert child.returncode == -signal.SIGKILL and cache.read_bytes() == original
    finally:
        if child.returncode is None:
            child.kill()
            child.wait()
    # SIGKILL cannot run destructors; orphan temporary files are safe to remove.
    for temporary in root.glob(".traceon-tmp-*"):
        temporary.unlink()
    env = dict(os.environ, LD_PRELOAD=shim, TRACEON_TEST_ATOMIC_FAULT="directory-sync")
    result = subprocess.run([builder, str(new), str(cache)] + storage, env=env, capture_output=True, timeout=10)
    assert result.returncode != 0 and b"replacement completed" in result.stderr
    assert cache.read_bytes() != original  # Post-commit durability error, not a failed replacement.
    assert cache.stat().st_mode & 0o777 == 0o640
    result = subprocess.run([driver, "restore", str(cache)], check=True, capture_output=True)
    assert b"Sequences: 2" in result.stdout
    assert not list(root.glob(".traceon-tmp-*"))
    env = dict(os.environ, LD_PRELOAD=shim, TRACEON_TEST_ATOMIC_FAULT="directory-sync-einval")
    result = subprocess.run([builder, str(new), str(cache)] + storage, env=env, capture_output=True, timeout=10)
    assert result.returncode == 0
    assert cache.stat().st_mode & 0o777 == 0o640
    result = subprocess.run([driver, "restore", str(cache)], check=True, capture_output=True)
    assert b"Sequences: 2" in result.stdout
    link = root / "link.traceon"
    link.symlink_to(cache.name)
    before = cache.read_bytes()
    assert subprocess.run([builder, str(old), str(link)] + storage, capture_output=True).returncode != 0
    assert link.is_symlink() and cache.read_bytes() == before
    subprocess.run([builder, str(old), str(cache)] + storage, check=True, capture_output=True)
    assert cache.read_bytes() == original
print("atomic save failure, interruption, replacement, and durability checks passed")
