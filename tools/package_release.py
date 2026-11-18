#!/usr/bin/env python3
"""Archive an installed package with file hashes and build-host provenance."""
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

source, output = (Path(p).resolve() for p in sys.argv[1:])
if not (source / 'bin').is_dir() or not (source / 'include' / 'TracEon').is_dir():
    raise SystemExit('expected a TracEon installation with tools and headers')
if output == source or source in output.parents:
    raise SystemExit('archive output must be outside the installation')
output.mkdir(parents=True, exist_ok=True)
name = f"traceon-{os.environ.get('PACKAGE_OS', platform.system())}-{os.environ.get('PACKAGE_ARCH', platform.machine())}"
files = {p.relative_to(source).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
         for p in sorted(source.rglob('*')) if p.is_file() and p.name != 'BUILD-MANIFEST.json'}
manifest = {'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            'working_tree': subprocess.check_output(['git', 'status', '--short'], text=True),
            'platform': platform.platform(), 'files_sha256': files,
            'compatibility': 'Built on the recorded host; system runtime dependencies are not bundled. Windows packages contain core tools/libraries, not the POSIX service.'}
(source / 'BUILD-MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
archive = Path(shutil.make_archive(str(output / name), 'zip' if platform.system() == 'Windows' else 'gztar', source))
(output / (archive.name + '.sha256')).write_text(hashlib.sha256(archive.read_bytes()).hexdigest() + '  ' + archive.name + '\n')
print(archive)
