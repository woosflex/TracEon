"""Standalone packaging check without a compiler or external services."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import zipfile

repo = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    install = root / 'install'
    (install / 'bin').mkdir(parents=True)
    (install / 'include/TracEon').mkdir(parents=True)
    (install / 'bin/tool').write_bytes(b'test executable bytes')
    (install / 'include/TracEon/test.h').write_text('// test header\n')
    command = [sys.executable, str(repo / 'tools/package_release.py'), str(install)]
    for _ in range(2):
        subprocess.run(command + [str(root / 'artifacts')], cwd=repo, check=True, capture_output=True)
    archive = next(p for p in (root / 'artifacts').iterdir() if not p.name.endswith('.sha256'))
    assert archive.with_name(archive.name + '.sha256').read_text().split()[0] == hashlib.sha256(archive.read_bytes()).hexdigest()
    if archive.suffix == '.zip':
        with zipfile.ZipFile(archive) as source:
            files = {name: source.read(name) for name in source.namelist() if not name.endswith('/')}
    else:
        with tarfile.open(archive) as source:
            files = {p.name.removeprefix('./'): source.extractfile(p).read() for p in source if p.isfile()}
    manifest = json.loads(files['BUILD-MANIFEST.json'])
    assert 'BUILD-MANIFEST.json' not in manifest['files_sha256']
    assert set(manifest['files_sha256']) == {'bin/tool', 'include/TracEon/test.h'}
    for name, digest in manifest['files_sha256'].items():
        assert hashlib.sha256(files[name]).hexdigest() == digest
    assert subprocess.run(command + [str(install / 'nested')], cwd=repo, capture_output=True).returncode != 0
print('release packaging checks passed')
